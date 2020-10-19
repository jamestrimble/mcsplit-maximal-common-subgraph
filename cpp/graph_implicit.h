#include <limits.h>
#include <stdbool.h>

#include <vector>

struct Graph {
    int n;
    std::vector<std::vector<unsigned int>> adjmat;
    std::vector<std::vector<int>> adjlists;
    std::vector<unsigned int> label;
    Graph(unsigned int n);
};

Graph readGraph(char* filename);

