#include "graph.h"

#include <algorithm>
#include <numeric>
#include <chrono>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <atomic>

#include <argp.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using std::vector;
using std::cout;
using std::endl;

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

enum Heuristic { min_max, min_product };

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Find a maximum clique in a graph in DIMACS format";
static char args_doc[] = "FILENAME1 FILENAME2";
static struct argp_option options[] = {
    {"verbose", 'v', 0, 0, "Verbose output"},
    {"connected", 'c', 0, 0, "Solve max common CONNECTED subgraph problem"},
    {"timeout", 't', "timeout", 0, "Specify a timeout (seconds)"},
    { 0 }
};

static struct {
    bool verbose;
    bool connected;
    char *filename1;
    char *filename2;
    int timeout;
    int arg_num;
} arguments;

static std::atomic<bool> abort_due_to_timeout;

void set_default_arguments() {
    arguments.verbose = false;
    arguments.connected = false;
    arguments.filename1 = NULL;
    arguments.filename2 = NULL;
    arguments.timeout = 0;
    arguments.arg_num = 0;
}

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    switch (key) {
        case 'v':
            arguments.verbose = true;
            break;
        case 'c':
            arguments.connected = true;
            break;
        case 't':
            arguments.timeout = std::stoi(arg);
            break;
        case ARGP_KEY_ARG:
            if (arguments.arg_num == 0) {
                arguments.filename1 = arg;
            } else if (arguments.arg_num == 1) {
                arguments.filename2 = arg;
            } else {
                argp_usage(state);
            }
            arguments.arg_num++;
            break;
        case ARGP_KEY_END:
            if (arguments.arg_num == 0)
                argp_usage(state);
            break;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

/*******************************************************************************
                                     Stats
*******************************************************************************/

unsigned long long nodes{ 0 };

/*******************************************************************************
                                 MCS functions
*******************************************************************************/

struct VtxPair {
    int v;
    int w;
    VtxPair(int v, int w): v(v), w(w) {}
};

using Iter = std::vector<int>::iterator;

struct Bidomain {
    Iter l_start;
    Iter r_start;
    Iter l_end;
    Iter r_end;
//    int l,        r;        // start indices of left and right sets
//    int left_len, right_len;
    bool is_adjacent;
    int X_count;
};

void show_current(const vector<VtxPair>& current)
{
//    cout << current.size() << "   ";
    for (unsigned int i=0; i<current.size(); i++) {
        cout << "  (" << current[i].v << " " << current[i].w << ")";
    }
    cout << std::endl;
}

void show(const vector<VtxPair>& current, const vector<Bidomain> &domains,
        const vector<int>& left, const vector<int>& right)
{
    cout << "Nodes: " << nodes << std::endl;
    cout << "Length of current assignment: " << current.size() << std::endl;
    cout << "Current assignment:";
    for (unsigned int i=0; i<current.size(); i++) {
        cout << "  (" << current[i].v << " -> " << current[i].w << ")";
    }
    cout << std::endl;
    for (unsigned int i=0; i<domains.size(); i++) {
        struct Bidomain bd = domains[i];
        cout << "Left  ";
        for (Iter it=bd.l_start; it!=bd.l_end; it++)
            cout << *it << " ";
        cout << std::endl;
        cout << "Right  ";
        for (Iter it=bd.r_start; it!=bd.r_end; it++)
            cout << *it << " ";
        cout << std::endl;
    }
    cout << "\n" << std::endl;
}

bool check_sol(const Graph & g0, const Graph & g1 , const vector<VtxPair> & solution) {
    return true;
    vector<bool> used_left(g0.n, false);
    vector<bool> used_right(g1.n, false);
    for (unsigned int i=0; i<solution.size(); i++) {
        struct VtxPair p0 = solution[i];
        if (used_left[p0.v] || used_right[p0.w])
            return false;
        used_left[p0.v] = true;
        used_right[p0.w] = true;
        if (g0.label[p0.v] != g1.label[p0.w])
            return false;
        for (unsigned int j=i+1; j<solution.size(); j++) {
            struct VtxPair p1 = solution[j];
            if (g0.adjmat[p0.v][p1.v] != g1.adjmat[p0.w][p1.w])
                return false;
        }
    }
    return true;
}

int find_and_remove_first_val(Bidomain & bd, vector<bool> & X) {
    for (Iter it=bd.l_start; it!=bd.l_end; it++) {
        int v = *it;
        if (!X[v]) {
            bd.l_end--;
            std::swap(*it, *bd.l_end);
            return v;
        }
    }
    return -1;
}

int select_bidomain(const vector<Bidomain>& domains, const vector<int> & left,
        int current_matching_size)
{
    for (unsigned int i=0; i<domains.size(); i++) {
        const Bidomain &bd = domains[i];
        if (bd.l_end - bd.l_start == bd.X_count)
            continue;
        if (arguments.connected && current_matching_size>0 && !bd.is_adjacent)
            continue;
        return i;
    }
    return -1;
}

// Returns iter to one-past-end of left part
Iter partition(Iter start, Iter end, const vector<unsigned int> & adjrow) {
    return std::partition(start, end,
            [&](const int elem){ return 0 != adjrow[elem]; });
}

vector<Bidomain> filter_domains(const vector<Bidomain> & d, vector<int> & left,
        vector<int> & right, const Graph & g0, const Graph & g1, int v, int w,
        vector<bool> & X)
{
    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    for (const Bidomain &old_bd : d) {
////        int l = old_bd.l;
////        int r = old_bd.r;
        Iter l_middle = partition(old_bd.l_start, old_bd.l_end, g0.adjmat[v]);
        Iter r_middle = partition(old_bd.r_start, old_bd.r_end, g1.adjmat[w]);
//        int left_len_noedge = old_bd.left_len - left_len;
//        int right_len_noedge = old_bd.right_len - right_len;
        if (l_middle != old_bd.l_end && r_middle != old_bd.r_end) {
            int X_count = 0;
            for (Iter it=l_middle; it!=old_bd.l_end; it++) {
                X_count += X[*it];
            }
            new_d.push_back({l_middle, r_middle, old_bd.l_end, old_bd.r_end,
                    old_bd.is_adjacent, X_count});
        }
        if (old_bd.l_start != l_middle && old_bd.r_start != r_middle) {
            int X_count = 0;
            for (Iter it=old_bd.l_start; it!=l_middle; it++) {
                X_count += X[*it];
            }
            new_d.push_back({old_bd.l_start, old_bd.r_start, l_middle, r_middle, true, X_count});
        }
    }
    return new_d;
}

// returns and Iter to the smallest value in arr that is >w.
// Assumption: such a value exists
// Assumption: arr contains no duplicates
// Assumption: arr does not contain INT_MAX
Iter iter_to_next_smallest(Iter start, Iter end, int w) {
    Iter retval;
    int smallest = INT_MAX;
    for (Iter it=start; it!=end; it++) {
        if (*it>w && *it<smallest) {
            smallest = *it;
            retval=it;
        }
    }
    return retval;
}

void solve(const Graph & g0, const Graph & g1,
        vector<VtxPair> & current, vector<Bidomain> & domains,
        vector<int> & left, vector<int> & right, unsigned int matching_size_goal,
        vector<bool> & X)
{
    if (abort_due_to_timeout)
        return;

    if (arguments.verbose) show(current, domains, left, right);
    nodes++;

    int bd_idx = select_bidomain(domains, left, current.size());
    if (bd_idx == -1) {
        bool is_maximal = true;
        if (arguments.connected && !current.empty()) {
            for (auto & bd : domains) {
                if (bd.X_count && bd.is_adjacent) {
                    is_maximal = false;
                    break;
                }
            }
        } else {
            for (auto & bd : domains) {
                if (bd.X_count) {
                    is_maximal = false;
                    break;
                }
            }
        }
        if (is_maximal) {
            show_current(current);
//            std::cout << 1 << std::endl;
        }
        return;
    }
    Bidomain &bd = domains[bd_idx];

    int v = find_and_remove_first_val(bd, X);

    // Try assigning v to each vertex w beginning at bd.r, in turn
    int w = -1;
    bd.r_end--;
    int num_r_vals = bd.r_end - bd.r_start;
    for (int i=0; i<=num_r_vals; i++) {
        Iter iter = iter_to_next_smallest(bd.r_start, bd.r_end+1, w);
        w = *iter;

        // swap w to the end of its colour class
        *iter = *bd.r_end;
        *bd.r_end = w;

        auto new_domains = filter_domains(domains, left, right, g0, g1, v, w, X);
        current.push_back(VtxPair(v, w));
        solve(g0, g1, current, new_domains, left, right, matching_size_goal, X);
        current.pop_back();
    }
    bd.l_end++;
    bd.r_end++;
    X[v] = true;
    ++bd.X_count;
    solve(g0, g1, current, domains, left, right, matching_size_goal, X);
    X[v] = false;
}

void mcs(const Graph & g0, const Graph & g1) {
    vector<int> left;  // the buffer of vertex indices for the left partitions
    vector<int> right;  // the buffer of vertex indices for the right partitions
    left.reserve(g0.n);
    right.reserve(g1.n);

    auto domains = vector<Bidomain> {};

    std::set<unsigned int> left_labels;
    std::set<unsigned int> right_labels;
    for (unsigned int label : g0.label) left_labels.insert(label);
    for (unsigned int label : g1.label) right_labels.insert(label);
    std::set<unsigned int> labels;  // labels that appear in both graphs
    std::set_intersection(std::begin(left_labels),
                          std::end(left_labels),
                          std::begin(right_labels),
                          std::end(right_labels),
                          std::inserter(labels, std::begin(labels)));

    // Create a bidomain for each label that appears in both graphs
    for (unsigned int label : labels) {
        Iter start_l = left.end();
        Iter start_r = right.end();

        for (int i=0; i<g0.n; i++)
            if (g0.label[i]==label)
                left.push_back(i);
        for (int i=0; i<g1.n; i++)
            if (g1.label[i]==label)
                right.push_back(i);

        domains.push_back({start_l, start_r, left.end(), right.end(), false, 0});
    }

    vector<VtxPair> current;
    vector<bool> X(g0.n);
    solve(g0, g1, current, domains, left, right, 1, X);
}

vector<int> calculate_degrees(const Graph & g) {
    vector<int> degree(g.n, 0);
    for (int v=0; v<g.n; v++) {
        for (int w=0; w<g.n; w++) {
            unsigned int mask = 0xFFFFu;
            if (g.adjmat[v][w] & mask) degree[v]++;
            if (g.adjmat[v][w] & ~mask) degree[v]++;  // inward edge, in directed case
        }
    }
    return degree;
}

int sum(const vector<int> & vec) {
    return std::accumulate(std::begin(vec), std::end(vec), 0);
}

int main(int argc, char** argv) {
    set_default_arguments();
    argp_parse(&argp, argc, argv, 0, 0, 0);

    struct Graph g0 = readGraph(arguments.filename1);
    struct Graph g1 = readGraph(arguments.filename2);

    std::thread timeout_thread;
    std::mutex timeout_mutex;
    std::condition_variable timeout_cv;
    abort_due_to_timeout.store(false);
    bool aborted = false;

    if (0 != arguments.timeout) {
        timeout_thread = std::thread([&] {
                auto abort_time = std::chrono::steady_clock::now() + std::chrono::seconds(arguments.timeout);
                {
                    /* Sleep until either we've reached the time limit,
                     * or we've finished all the work. */
                    std::unique_lock<std::mutex> guard(timeout_mutex);
                    while (! abort_due_to_timeout.load()) {
                        if (std::cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                            /* We've woken up, and it's due to a timeout. */
                            aborted = true;
                            break;
                        }
                    }
                }
                abort_due_to_timeout.store(true);
                });
    }

    auto start = std::chrono::steady_clock::now();

//    vector<int> g0_deg = calculate_degrees(g0);
//    vector<int> g1_deg = calculate_degrees(g1);
//
//    vector<int> vv0(g0.n);
//    std::iota(std::begin(vv0), std::end(vv0), 0);
//    bool g1_dense = sum(g1_deg) > g1.n*(g1.n-1);
//    std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
//        return g1_dense ? (g0_deg[a]<g0_deg[b]) : (g0_deg[a]>g0_deg[b]);
//    });
//    vector<int> vv1(g1.n);
//    std::iota(std::begin(vv1), std::end(vv1), 0);
//    bool g0_dense = sum(g0_deg) > g0.n*(g0.n-1);
//    std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) {
//        return g0_dense ? (g1_deg[a]<g1_deg[b]) : (g1_deg[a]>g1_deg[b]);
//    });
//
//    struct Graph g0_sorted = induced_subgraph(g0, vv0);
//    struct Graph g1_sorted = induced_subgraph(g1, vv1);

    mcs(g0, g1);

    auto stop = std::chrono::steady_clock::now();
    auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    /* Clean up the timeout thread */
    if (timeout_thread.joinable()) {
        {
            std::unique_lock<std::mutex> guard(timeout_mutex);
            abort_due_to_timeout.store(true);
            timeout_cv.notify_all();
        }
        timeout_thread.join();
    }

////    for (int i=0; i<g0.n; i++)
////        for (unsigned int j=0; j<solution.size(); j++)
////            if (solution[j].v == i)
////                cout << "(" << solution[j].v << " -> " << solution[j].w << ") ";
////    cout << std::endl;

    cout << "Nodes:                      " << nodes << endl;
    cout << "CPU time (ms):              " << time_elapsed << endl;
    if (aborted)
        cout << "TIMEOUT" << endl;
}

