#include <limits.h>
#include <stdbool.h>

#include <unordered_set>
#include <vector>

struct Graph {
    int n;
    std::vector<std::unordered_set<int>> adjsets;
    std::vector<unsigned int> label;
    Graph(unsigned int n);
};

Graph readGraph(char* filename);

