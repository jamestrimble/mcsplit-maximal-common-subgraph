import sys
import itertools
import networkx as nx


def _label_function(label):
    if label is None:
        return lambda data: 1
    if callable(label):
        return label
    return lambda data: data[label]


def get_edge_label(G, edge_label_fun, node_a, node_b):
    if G.has_edge(node_a, node_b):
        return edge_label_fun(G.edges[node_a, node_b])
    else:
        return None


class LabelClass(object):
    __slots__ = ["G_nodes", "H_nodes", "is_adjacent"]

    def __init__(self, is_adjacent):
        self.G_nodes = []
        self.H_nodes = []
        self.is_adjacent = is_adjacent

    def __repr__(self):
        return f"<{self.G_nodes}, {self.H_nodes}, {self.is_adjacent}>"


class PartitioningMCISFinder(object):
    __slots__ = ["G", "H", "connected", "node_label_fun", "edge_label_fun", "count"]

    def __init__(self, G, H, connected, node_label_fun, edge_label_fun):
        self.G = G
        self.H = H
        self.connected = connected
        self.node_label_fun = node_label_fun
        self.edge_label_fun = edge_label_fun

    def refine_label_classes(self, label_classes, v, w):
        new_label_classes = []
        for lc in label_classes:
            label_to_new_lc = {}
            for u in lc.G_nodes:
                edge_label = get_edge_label(self.G, self.edge_label_fun, v, u)
                if edge_label not in label_to_new_lc:
                    is_adjacent = lc.is_adjacent or self.G.has_edge(v, u)
                    label_to_new_lc[edge_label] = LabelClass(is_adjacent)
                label_to_new_lc[edge_label].G_nodes.append(u)
            for u in lc.H_nodes:
                if u == w:
                    continue
                edge_label = get_edge_label(self.H, self.edge_label_fun, w, u)
                if edge_label in label_to_new_lc:
                    label_to_new_lc[edge_label].H_nodes.append(u)
            for new_lc in label_to_new_lc.values():
                if new_lc.H_nodes:
                    new_label_classes.append(new_lc)
        return new_label_classes

    def select_label_class(self, label_classes, assignment_count):
        if self.connected and assignment_count > 0:
            candidates = [lc for lc in label_classes if lc.is_adjacent]
        else:
            candidates = label_classes
        if not candidates:
            return None
        return min(candidates, key=lambda lc: max(len(lc.G_nodes), len(lc.H_nodes)))

    def add_to_rejected(self, rejected, newly_rejected):
        set_H = set(newly_rejected.H_nodes)
        for i, r in enumerate(rejected):
            if set(r.H_nodes) == set_H:
                newly_rejected.G_nodes.extend(r.G_nodes)
                return rejected[:i] + rejected[i + 1 :] + [newly_rejected]
        return rejected + [newly_rejected]

    def search(self, label_classes, rejected, assignments):
        label_class = self.select_label_class(label_classes, len(assignments))
        if label_class is None:
            if assignments:
                if self.connected:
                    is_maximal = not any(reject.is_adjacent for reject in rejected)
                else:
                    is_maximal = not rejected
                if is_maximal:
                    self.count += 1
                    print(self.count)
            return
        v = label_class.G_nodes.pop()
        H_nodes = label_class.H_nodes[:]
        for w in H_nodes:
            label_class.H_nodes[:] = [u for u in H_nodes if u != w]
            assignments[v] = w
            new_label_classes = self.refine_label_classes(label_classes, v, w)
            new_rejected = self.refine_label_classes(rejected, v, w)
            self.search(new_label_classes, new_rejected, assignments)
            del assignments[v]
        label_class.H_nodes[:] = H_nodes
        keys = set(assignments.keys())
        vals = set(assignments.values())
        new_label_classes = [lc for lc in label_classes if lc.G_nodes]
        newly_rejected = LabelClass(label_class.is_adjacent)
        newly_rejected.G_nodes.append(v)
        newly_rejected.H_nodes.extend(H_nodes)
        new_rejected = self.add_to_rejected(rejected, newly_rejected)
        self.search(new_label_classes, new_rejected, assignments)

    def get_node_labels(self, G):
        return (self.node_label_fun(G.nodes[node]) for node in G.nodes())

    def find_maximal_common_subgraphs(self):
        G_labels = set(self.get_node_labels(self.G))
        H_labels = set(self.get_node_labels(self.H))
        label_classes = {label: LabelClass(False) for label in G_labels & H_labels}
        for v in self.G.nodes():
            label = self.node_label_fun(self.G.nodes[v])
            if label in label_classes:
                label_classes[label].G_nodes.append(v)
        for v in self.H.nodes():
            label = self.node_label_fun(self.H.nodes[v])
            if label in label_classes:
                label_classes[label].H_nodes.append(v)
        self.count = 0
        self.search(label_classes.values(), [], {})
        return self.count


def check_valid_graph_types(G, H):
    if G.is_directed() or H.is_directed():
        msg = "not implemented for directed graphs"
        raise nx.NetworkXNotImplemented(msg)
    if G.is_multigraph() or H.is_multigraph():
        msg = "not implemented for multigraphs"
        raise nx.NetworkXNotImplemented(msg)
    if nx.number_of_selfloops(G) or nx.number_of_selfloops(H):
        msg = "not implemented for graphs with self-loops"
        raise nx.NetworkXNotImplemented(msg)


def check_valid_labels(G, H, node_label, edge_label):
    """Raise an exception if labels are invalid

    An exception is raised if any edge label is None or if accessing any label
    raises an exception.
    """
    for graph in G, H:
        for edge in graph.edges():
            if edge_label(graph.edges[edge]) is None:
                raise ValueError(f"Edge {edge} has label None")
        for node in graph.nodes():
            # Make sure that accessing the node does not raise an exception
            node_label(graph.nodes[node])


def maximal_common_subgraphs(G, H, connected=False, node_label=None, edge_label=None):
    check_valid_graph_types(G, H)
    node_label = _label_function(node_label)
    edge_label = _label_function(edge_label)
    check_valid_labels(G, H, node_label, edge_label)
    return PartitioningMCISFinder(
        G, H, connected, node_label, edge_label
    ).find_maximal_common_subgraphs()


def read_instance(filename):
    with open(filename, "r") as f:
        lines = [line.strip().split() for line in f.readlines()]
    n = int(lines[0][0])
    edge_count = int(lines[0][1])
    labels = lines[1]
    e_lines = lines[2:]
    G = nx.Graph()
    for v in range(n):
        G.add_node(v, label=labels[v])
    for e in e_lines:
        v, w = int(e[0]), int(e[1])
        G.add_edge(v, w, xyz=(v, w))
    return G


def main(G0, G1):
    print(maximal_common_subgraphs(G0, G1, connected=True, node_label="label"))


if __name__ == "__main__":
    main(read_instance(sys.argv[1]), read_instance(sys.argv[2]))
