#include <limits.h>
#include <stdbool.h>

#include <vector>

struct Graph {
    int n;
    std::vector<std::vector<unsigned char>> adjmat;
    std::vector<unsigned int> label;
    Graph(unsigned int n);
};

Graph induced_subgraph(struct Graph& g, std::vector<int> vv);

Graph readGraph(char* filename);

