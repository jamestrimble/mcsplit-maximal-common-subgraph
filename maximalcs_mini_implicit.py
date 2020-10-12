import sys


class Graph(object):
    def __init__(self, n, labels):
        self.n = n
        self.adjmat = [[0] * n for _ in range(n)]
        self.adj_lists = [[] for _ in range(n)]
        self.labels = labels

    def add_edge(self, v, w):
        if not self.adjmat[v][w]:
            self.adjmat[v][w] = 1
            self.adjmat[w][v] = 1
            self.adj_lists[v].append(w)
            self.adj_lists[w].append(v)


class LabelClass(object):
    __slots__ = ['G_nodes', 'H_nodes', 'is_adjacent', 'X_count']

    def __init__(self, is_adjacent):
        self.G_nodes = []
        self.H_nodes = []
        self.is_adjacent = is_adjacent
        self.X_count = 0

    def __repr__(self):
        return f"<{self.G_nodes}, {self.H_nodes}, {self.is_adjacent}, {self.X_count}>"


def make_adjacent_label_classes(G, H, left, right, X, D_G, D_H):
    new_label_classes = []
    label_to_new_lc = {}
    for u in left:
        label = G.labels[u]
        if label not in label_to_new_lc:
            label_to_new_lc[label] = LabelClass(True)
        new_lc = label_to_new_lc[label]
        new_lc.G_nodes.append(u)
        new_lc.X_count += X[u]
    for u in right:
        label = H.labels[u]
        if label not in label_to_new_lc:
            continue
        new_lc = label_to_new_lc[label]
        new_lc.H_nodes.append(u)
    for new_lc in label_to_new_lc.values():
        if new_lc.G_nodes and new_lc.H_nodes:
            new_label_classes.append(new_lc)
    return new_label_classes


def refine_label_classes(G, H, label_classes, v, w, X):
    new_label_classes = []
    for lc in label_classes:
        if not lc.G_nodes or not lc.H_nodes:
            # an optimisation
            continue
        label_to_new_lc = [LabelClass(lc.is_adjacent), LabelClass(True)]
        G_adjrow = G.adjmat[v]
        for u in lc.G_nodes:
            new_lc = label_to_new_lc[G_adjrow[u]]
            new_lc.G_nodes.append(u)
            new_lc.X_count += X[u]
        H_adjrow = H.adjmat[w]
        for u in lc.H_nodes:
            new_lc = label_to_new_lc[H_adjrow[u]]
            new_lc.H_nodes.append(u)
        for new_lc in label_to_new_lc:
            if new_lc.G_nodes and new_lc.H_nodes:
                new_label_classes.append(new_lc)
    return new_label_classes


def select_label_class(label_classes):
    for lc in label_classes:
        if lc.is_adjacent and len(lc.G_nodes) > lc.X_count:
            return lc
    return None


def search(G, H, label_classes, assignments, X, D_G, D_H):
    # Returns number of maximal CISs found in this and its recursive calls
    label_class = select_label_class(label_classes)
    if label_class is None:
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
    left = [u for u in G.adj_lists[v] if u in D_G]
#    print("assignments", assignments)
#    print("left", left)
#    print("D_G", sorted(D_G))
    for u in left:
        D_G.remove(u)
#    print("v!", v)
    for w in H_nodes:
        label_class.H_nodes[:] = [u for u in H_nodes if u != w]
        assignments[v] = w
        right = [u for u in H.adj_lists[w] if u in D_H]
        for u in right:
            D_H.remove(u)
        new_label_classes = (
            refine_label_classes(G, H, label_classes, v, w, X) +
            make_adjacent_label_classes(G, H, left, right, X, D_G, D_H)
        )
        count += search(G, H, new_label_classes, assignments, X, D_G, D_H)
        del assignments[v]
        for u in right:
            D_H.add(u)
    for u in left:
        D_G.add(u)
    label_class.G_nodes.append(v)
    label_class.H_nodes[:] = H_nodes
    X[v] = 1
    label_class.X_count += 1
    count += search(G, H, label_classes, assignments, X, D_G, D_H)
    X[v] = 0
    return count


def start_search(G, H, label_classes):
    # Returns number of maximal CISs found
    if not label_classes:
        # An important edge case.  TODO write in thesis about how other algorithms handle it
        print({})
        return 1
    count = 0
    X = [0] * G.n
    D_G = set(range(G.n))
    D_H = set(range(H.n))
    for label_class in label_classes:
#        print("! ", label_class.H_nodes)
        for v in label_class.G_nodes:
            D_G.remove(v)
#            print("!!!!", sorted(D_G))
            left = G.adj_lists[v]
            for u in G.adj_lists[v]:
                D_G.remove(u)
            for w in label_class.H_nodes:
#                print("w", w)
                right = H.adj_lists[w]
                D_H.remove(w)
                for u in H.adj_lists[w]:
                    D_H.remove(u)
                lcs = make_adjacent_label_classes(G, H, left, right, X, D_G, D_H)
#                print("!!!", sorted(D_G))
#                print("   ", v, w)
#                print(lcs)
                count += search(G, H, lcs, {v: w}, X, D_G, D_H)
                D_H.add(w)
                for u in H.adj_lists[w]:
                    D_H.add(u)
            D_G.add(v)
            for u in G.adj_lists[v]:
                D_G.add(u)
            X[v] = 1
    return count


def find_maximal_common_subgraphs(G, H):
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
    return start_search(G, H, label_classes.values())


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
    print(find_maximal_common_subgraphs(G, H))
