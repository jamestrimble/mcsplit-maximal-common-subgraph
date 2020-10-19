# Using a McSplit-like algorithm to solve *maximal* common subgraph

This is just a proof of concept.  The file format used is described
in the README of https://github.com/veluca93/parallel_enum/tree/bccliques .

To run:

```
python3 maximalcs.py g1.txt g2.txt
```

To generate:  The following will generate a random graph with 10
vertices, edge probability .5, and 3 distinct vertex labels.

```
python3 gen.py 10 .5 3
```

## C++ Versions

The `cpp` directory uses arrays of arrays for adjacency matrices, while
the `cpp-using-sets` directory uses vectors of sets.  The latter is
probably preferable for comparison with Versari's code, and the performance
penalty seems fairly small.

The `mcsp_implicit_b` variant avoids the `D_G` and `D_H` sets, so it's
a bit more implicit. It's maybe closer to Versari's `koch_implicit`
program from his batchelor's thesis, so it might be more interesting for
comparison with it.
