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
    __slots__ = ['G_nodes', 'H_nodes', 'X_count']

    def __init__(self):
        self.G_nodes = []
        self.H_nodes = []
        self.X_count = 0

    def __repr__(self):
        return f"<{self.G_nodes}, {self.H_nodes}, {self.X_count}>"


def make_adjacent_label_classes(G, H, left, right, X):
    new_label_classes = []
    label_to_new_lc = {}
    for u in left:
        label = G.labels[u]
        if label not in label_to_new_lc:
            label_to_new_lc[label] = LabelClass()
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
        label_to_new_lc = [LabelClass(), LabelClass()]
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
        if len(lc.G_nodes) > lc.X_count:
            return lc
    return None


def set_all(vals, bools):
    for val in vals:
        bools[val] = True


def unset_all(vals, bools):
    for val in vals:
        bools[val] = False


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
    left = [u for u in G.adj_lists[v] if D_G[u]]
    unset_all(left, D_G)
    for w in H_nodes:
        label_class.H_nodes[:] = [u for u in H_nodes if u != w]
        assignments[v] = w
        right = [u for u in H.adj_lists[w] if D_H[u]]
        unset_all(right, D_H)
        new_label_classes = (
            refine_label_classes(G, H, label_classes, v, w, X) +
            make_adjacent_label_classes(G, H, left, right, X)
        )
        count += search(G, H, new_label_classes, assignments, X, D_G, D_H)
        del assignments[v]
        set_all(right, D_H)
    set_all(left, D_G)
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
    D_G = [True] * G.n
    D_H = [True] * H.n
    for label_class in label_classes:
        for v in label_class.G_nodes:
            left = G.adj_lists[v]
            D_G[v] = False
            unset_all(G.adj_lists[v], D_G)
            for w in label_class.H_nodes:
                right = H.adj_lists[w]
                D_H[w] = False
                unset_all(H.adj_lists[w], D_H)
                lcs = make_adjacent_label_classes(G, H, left, right, X)
                count += search(G, H, lcs, {v: w}, X, D_G, D_H)
                D_H[w] = True
                set_all(H.adj_lists[w], D_H)
            D_G[v] = True
            set_all(G.adj_lists[v], D_G)
            X[v] = 1
    return count


def find_maximal_common_subgraphs(G, H):
    G_labels = set(G.labels)
    H_labels = set(H.labels)
    label_classes = {label: LabelClass() for label in G_labels & H_labels}
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
