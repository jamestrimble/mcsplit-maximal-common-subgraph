#include "graph_implicit.h"

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

static char doc[] = "Find maximal common connected induced subgraphs";
static char args_doc[] = "FILENAME1 FILENAME2";
static struct argp_option options[] = {
    {"verbose", 'v', 0, 0, "Verbose output"},
    {"quiet", 'q', 0, 0, "Quiet output; this runs faster"},
    {"timeout", 't', "timeout", 0, "Specify a timeout (seconds)"},
    { 0 }
};

static struct {
    bool verbose;
    bool quiet;
    char *filename1;
    char *filename2;
    int timeout;
    int arg_num;
} arguments;

static std::atomic<bool> abort_due_to_timeout;

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    switch (key) {
        case 'v':
            arguments.verbose = true;
            break;
        case 'q':
            arguments.quiet = true;
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
unsigned long long solution_count{ 0 };

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
    int X_count;
};

void show_current(const vector<VtxPair>& current)
{
    ++solution_count;
    if (arguments.quiet) {
        return;
    }
    for (unsigned int i=0; i<current.size(); i++) {
        cout << "  (" << current[i].v << " " << current[i].w << ")";
    }
    cout << std::endl;
}

void show(const vector<VtxPair>& current, const vector<Bidomain> &domains)
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
            if (g0.adjsets[p0.v].count(p1.v) != g1.adjsets[p0.w].count(p1.w))
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

int select_bidomain(const vector<Bidomain>& domains)
{
    for (unsigned int i=0; i<domains.size(); i++) {
        const Bidomain &bd = domains[i];
        if (bd.l_end - bd.l_start == bd.X_count)
            continue;
        return i;
    }
    return -1;
}

// Returns iter to one-past-end of left part
Iter partition(Iter start, Iter end, const std::unordered_set<int> & adjset) {
    return std::partition(start, end,
            [&](const int elem){ return adjset.count(elem); });
}

vector<Bidomain> filter_domains(const vector<Bidomain> & d, vector<int> & left,
        vector<int> & right, const Graph & g0, const Graph & g1, int v, int w,
        vector<bool> & X)
{
    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    for (const Bidomain &old_bd : d) {
        Iter l_middle = partition(old_bd.l_start, old_bd.l_end, g0.adjsets[v]);
        Iter r_middle = partition(old_bd.r_start, old_bd.r_end, g1.adjsets[w]);
        if (l_middle != old_bd.l_end && r_middle != old_bd.r_end) {
            int X_count = 0;
            for (Iter it=l_middle; it!=old_bd.l_end; it++) {
                X_count += X[*it];
            }
            new_d.push_back({l_middle, r_middle, old_bd.l_end, old_bd.r_end, X_count});
        }
        if (old_bd.l_start != l_middle && old_bd.r_start != r_middle) {
            int X_count = 0;
            for (Iter it=old_bd.l_start; it!=l_middle; it++) {
                X_count += X[*it];
            }
            new_d.push_back({old_bd.l_start, old_bd.r_start, l_middle, r_middle, X_count});
        }
    }
    return new_d;
}

void add_adjacent_label_classes(vector<Bidomain> & new_d,
        vector<int> & left, vector<int> & right, const Graph & g0, const Graph & g1,
        vector<bool> & X)
{
    if (left.empty() || right.empty()) {
        return;
    }
    std::sort(left.begin(), left.end(), [&](int a, int b)
            { return g0.label[a] < g0.label[b]; });
    std::sort(right.begin(), right.end(), [&](int a, int b)
            { return g1.label[a] < g1.label[b]; });

    auto l = left.begin();
    auto r = right.begin();
    while (l!=left.end() && r!=right.end()) {
        unsigned int left_label = g0.label[*l];
        unsigned int right_label = g1.label[*r];
        if (left_label < right_label) {
            l++;
        } else if (left_label > right_label) {
            r++;
        } else {
            int X_count = 0;
            auto lmin = l;
            auto rmin = r;
            do {
                if (X[*l]) {
                    ++X_count;   // TODO: avoid branch?
                }
                l++;
            } while (l!=left.end() && g0.label[*l]==left_label);
            do { r++; } while (r!=right.end() && g1.label[*r]==left_label);
            new_d.push_back({lmin, rmin, l, r, X_count});
        }
    }
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

void set_all(const vector<int> & vals, vector<bool> & bools)
{
    for (int v : vals) {
        bools[v] = true;
    }
}

void unset_all(const vector<int> & vals, vector<bool> & bools)
{
    for (int v : vals) {
        bools[v] = false;
    }
}

// TODO: pass X by reference?
void solve(const Graph & g0, const Graph & g1,
        vector<VtxPair> & current, vector<Bidomain> & domains,
        vector<bool> & X, vector<bool> & D_G, vector<bool> & D_H)
{
    if (abort_due_to_timeout)
        return;

    if (arguments.verbose) show(current, domains);
    nodes++;

    if (domains.empty()) {
        // current is maximal
        show_current(current);
        return;
    }

    vector<int> left;
    vector<int> right;
    vector<int> added_to_X;

    for (auto & bd : domains) {
        int v = -1;
        bd.l_end--;
        int num_l_vals = bd.l_end - bd.l_start;
        for (int i=0; i<=num_l_vals; i++) {
            Iter iter = iter_to_next_smallest(bd.l_start, bd.l_end+1, v);
            v = *iter;

            if (X[v]) {
                continue;
            }

            // swap v to the end of its colour class
            *iter = *bd.l_end;
            *bd.l_end = v;

            left.clear();
            for (int u : g0.adjlists[v]) {
                if (D_G[u]) {
                    left.push_back(u);
                    D_G[u] = false;
                }
            }

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

                right.clear();
                for (int u : g1.adjlists[w]) {
                    if (D_H[u]) {
                        right.push_back(u);
                        D_H[u] = false;
                    }
                }

                auto new_domains = filter_domains(domains, left, right, g0, g1, v, w, X);
                add_adjacent_label_classes(new_domains, left, right, g0, g1, X);
                current.push_back(VtxPair(v, w));
                solve(g0, g1, current, new_domains, X, D_G, D_H);
                current.pop_back();
                set_all(right, D_H);
            }
            set_all(left, D_G);
            bd.r_end++;
            X[v] = true;
            added_to_X.push_back(v);
            ++bd.X_count;
        }
        bd.l_end++;
    }
    for (int v : added_to_X) {
        X[v] = false;
    }
}

void start_search(const Graph & g0, const Graph & g1,
        vector<Bidomain> & domains)
{
    vector<VtxPair> current;
    if (domains.empty()) {
        // an edge case
        show_current(current);
        return;
    }
    vector<bool> X(g0.n);
    vector<bool> D_G(g0.n, true);
    vector<bool> D_H(g1.n, true);
    for (Bidomain & bd : domains) {
        for (Iter it=bd.l_start; it!=bd.l_end; it++) {
            int v = *it;
            auto left = g0.adjlists[v];
            D_G[v] = false;
            unset_all(left, D_G);
            for (Iter it=bd.r_start; it!=bd.r_end; it++) {
                int w = *it;
                auto right = g1.adjlists[w];
                D_H[w] = false;
                unset_all(right, D_H);
                vector<Bidomain> domains;
                add_adjacent_label_classes(domains, left, right, g0, g1, X);
                current.push_back(VtxPair(v, w));
                solve(g0, g1, current, domains, X, D_G, D_H);
                current.pop_back();
                D_H[w] = true;
                set_all(right, D_H);
            }
            D_G[v] = true;
            set_all(left, D_G);
            X[v] = true;
        }
    }
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

        domains.push_back({start_l, start_r, left.end(), right.end(), 0});
    }

    start_search(g0, g1, domains);
}

int sum(const vector<int> & vec) {
    return std::accumulate(std::begin(vec), std::end(vec), 0);
}

int main(int argc, char** argv) {
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

    cout << "Solutions:                  " << solution_count << endl;
    cout << "Nodes:                      " << nodes << endl;
    cout << "CPU time (ms):              " << time_elapsed << endl;
    if (aborted)
        cout << "TIMEOUT" << endl;
}

