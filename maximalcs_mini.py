import sys
import itertools


class Graph(object):
    def __init__(self, n, labels):
        self.n = n
        self.adjmat = [[0] * n for _ in range(n)]
        self.labels = labels

    def add_edge(self, v, w):
        self.adjmat[v][w] = 1
        self.adjmat[w][v] = 1


class LabelClass(object):
    def __init__(self, is_adjacent):
        self.G_nodes = []
        self.H_nodes = []
        self.is_adjacent = is_adjacent
        self.X_count = 0

    def __repr__(self):
        return f"<{self.G_nodes}, {self.H_nodes}, {self.is_adjacent}, {self.X_count}>"


def refine_label_classes(G, H, label_classes, v, w, X):
    new_label_classes = []
    for lc in label_classes:
        label_to_new_lc = [LabelClass(lc.is_adjacent), LabelClass(True)]
        for u in lc.G_nodes:
            edge_label = G.adjmat[v][u]
            label_to_new_lc[edge_label].G_nodes.append(u)
            label_to_new_lc[edge_label].X_count += X[u]
        for u in lc.H_nodes:
            edge_label = H.adjmat[w][u]
            label_to_new_lc[edge_label].H_nodes.append(u)
        for new_lc in label_to_new_lc:
            if new_lc.G_nodes and new_lc.H_nodes:
                new_label_classes.append(new_lc)
    return new_label_classes


def select_label_class(connected, label_classes, assignment_count):
    if connected and assignment_count > 0:
        candidates = [
            lc
            for lc in label_classes
            if lc.is_adjacent and len(lc.G_nodes) > lc.X_count
        ]
    else:
        candidates = [lc for lc in label_classes if len(lc.G_nodes) > lc.X_count]
    if not candidates:
        return None
    return min(
        candidates, key=lambda lc: max(len(lc.G_nodes) - lc.X_count, len(lc.H_nodes))
    )


def search(G, H, connected, label_classes, assignments, X):
    # Returns number of maximal CISs found in this and its recursive calls
    label_class = select_label_class(connected, label_classes, len(assignments))
    if label_class is None:
        if connected and assignments:
            is_maximal = not any(lc.X_count and lc.is_adjacent for lc in label_classes)
        else:
            is_maximal = not any(lc.X_count for lc in label_classes)
        if is_maximal:
            print(assignments)
            return 1
        else:
            return 0
    for i, v in enumerate(label_class.G_nodes):
        if not X[v]:
            break
    del label_class.G_nodes[i]
    H_nodes = label_class.H_nodes[:]
    count = 0
    for w in H_nodes:
        label_class.H_nodes[:] = [u for u in H_nodes if u != w]
        assignments[v] = w
        new_label_classes = refine_label_classes(G, H, label_classes, v, w, X)
        count += search(G, H, connected, new_label_classes, assignments, X)
        del assignments[v]
    label_class.G_nodes.append(v)
    label_class.H_nodes[:] = H_nodes
    X[v] = 1
    label_class.X_count += 1
    count += search(G, H, connected, label_classes, assignments, X)
    X[v] = 0
    return count


def find_maximal_common_subgraphs(G, H, connected=False):
    G_labels = set(G.labels)
    H_labels = set(H.labels)
    label_classes = {label: LabelClass(False) for label in G_labels & H_labels}
    for v in range(G.n):
        label = G.labels[v]
        if label in label_classes:
            label_classes[label].G_nodes.append(v)
    for v in range(H.n):
        label = H.labels[v]
        if label in label_classes:
            label_classes[label].H_nodes.append(v)
    return search(G, H, connected, label_classes.values(), {}, [0] * G.n)


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
    print(find_maximal_common_subgraphs(G, H, connected=True))
