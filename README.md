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
