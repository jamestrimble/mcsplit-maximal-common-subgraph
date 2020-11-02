import sys


class Graph(object):
    def __init__(self, n, labels):
        self.n = n
        self.adjmat = [[0] * n for _ in range(n)]
        self.labels = labels

    def add_edge(self, v, w):
        self.adjmat[v][w] = 1
        self.adjmat[w][v] = 1


def reduce(F, H):
    F_labels = set(F.labels)
    H_labels = set(H.labels)
    G_vertices = []
    for label in F_labels & H_labels:
        F_vertices = [v for v in range(F.n) if F.labels[v] == label]
        H_vertices = [v for v in range(H.n) if H.labels[v] == label]
        for v in F_vertices:
            for w in H_vertices:
                G_vertices.append((v, w))
    # g is isomorphic to G, but has integers instead of vertex pairs as vertices
    g_adj = {v: set() for v in range(len(G_vertices))}
    for i, (u1, u2) in enumerate(G_vertices):
        for j, (w1, w2) in enumerate(G_vertices):
            if j <= i:
                continue
            if u1 == w1:
                continue
            if u2 == w2:
                continue
            if F.adjmat[u1][w1] != H.adjmat[u2][w2]:
                continue
            g_adj[i].add(j)
            g_adj[j].add(i)

    n = len(G_vertices)
    m = sum(len(lst) for lst in g_adj.values()) // 2
    print(n, m, 1)
    for v in range(n):
        print(" ".join(str(1 + w) for w in sorted(g_adj[v])))


def read_instance(filename):
    with open(filename, "r") as f:
        lines = [line.strip().split() for line in f.readlines()]
    n = int(lines[0][0])
    edge_count = int(lines[0][1])
    G = Graph(n, lines[1])
    for e in lines[2:]:
        G.add_edge(int(e[0]), int(e[1]))
    return G


if __name__ == "__main__":
    G = read_instance(sys.argv[1])
    H = read_instance(sys.argv[2])
    reduce(G, H)
