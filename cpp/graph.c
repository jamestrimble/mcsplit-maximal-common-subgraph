#include "graph.h"

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>

constexpr int BITS_PER_UNSIGNED_INT (CHAR_BIT * sizeof(unsigned int));

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

Graph::Graph(unsigned int n) {
    this->n = n;
    label = std::vector<unsigned int>(n, 0u);
    adjmat = {n, std::vector<unsigned int>(n, false)};
}

Graph induced_subgraph(struct Graph& g, std::vector<int> vv) {
    Graph subg(vv.size());
    for (int i=0; i<subg.n; i++)
        for (int j=0; j<subg.n; j++)
            subg.adjmat[i][j] = g.adjmat[vv[i]][vv[j]];

    for (int i=0; i<subg.n; i++)
        subg.label[i] = g.label[vv[i]];
    return subg;
}

void add_edge(Graph& g, int v, int w, bool directed=false, unsigned int val=1) {
    if (v != w) {
        if (directed) {
            g.adjmat[v][w] |= val;
            g.adjmat[w][v] |= (val<<16);
        } else {
            g.adjmat[v][w] = val;
            g.adjmat[w][v] = val;
        }
    } else {
        // To indicate that a vertex has a loop, we set the most
        // significant bit of its label to 1
        g.label[v] |= (1u << (BITS_PER_UNSIGNED_INT-1));
    }
}

struct Graph readGraph(char* filename) {
    FILE* f;
    
    if ((f=fopen(filename, "r"))==NULL)
        fail("Cannot open file");

    int n, m;
    fscanf(f, "%d%d", &n, &m);
    struct Graph g(n);
    for (int i=0; i<n; i++) {
        int label;
        fscanf(f, "%d", &label);
        g.label[i] |= label;
    }
    for (int i=0; i<m; i++) {
        int v, w;
        fscanf(f, "%d%d", &v, &w);
        add_edge(g, v, w);
    }

    return g;
}
