
class MAGE:

    def __init__(self):
        self.ani_score_gca = None
        self.ani_score_gcf = None
        self.ani_score_self = None
        self.handle_to_genome = None
        self.ani_graph = None

    def load_ani_data(self, path_gca, path_gcf, path_self):
        self.ani_score_gca = self.read_ani_out(path_gca)
        self.ani_score_gcf = self.read_ani_out(path_gcf)
        self.ani_score_self = self.read_ani_out(path_self)

    @staticmethod
    def read_ani_out(f):
        ani_score = {}
        with open(f, 'r') as fh:
            line = fh.readline()
            while line:
                _p = line.split('\t')
                g1, g2, score, k1, k2 = _p
                if g1 not in ani_score:
                    ani_score[g1] = {}
                ani_score[g1][g2] = (score, k1, k2)
                line = fh.readline()
        return ani_score

    @staticmethod
    def get_gtdb_genome(s):
        return s.split('/')[-1][:-4]

    @staticmethod
    def get_handle_name(handle_ref):
        return handle_to_genome[handle_ref.split('/')[-1]]

    @staticmethod
    def load_data2(ani_scores, nodes, edges):
        for e1 in ani_scores:
            node1 = self.get_handle_name(e1)
            nodes.add(node1)
            for e2 in ani_scores[e1]:
                node2 = self.get_gtdb_genome(e2)
                nodes.add(node2)
                if e1 != e2:
                    score, p1, p2 = ani_scores[e1][e2]
                    p_1_2 = (node1, node2)
                    p_2_1 = (node2, node1)
                    if not p_1_2 in edges and not p_2_1 in edges:
                        edges[p_1_2] = [score]
                    elif p_1_2 in edges:
                        edges[p_1_2].append(score)
                    elif p_2_1 in edges:
                        edges[p_2_1].append(score)
                    else:
                        print('this should never print')

    def build_graph(self):
        nodes = set()
        edges = {}
        load_data2(self.ani_score_gca, nodes, edges)
        print(len(nodes), len(edges))
        load_data2(self.ani_score_gcf, nodes, edges)
        print(len(nodes), len(edges))

        nodes_self = set()
        for e1 in self.ani_score_self:
            node1 = get_handle_name(e1)
            nodes_self.add(node1)
            nodes.add(node1)
            for e2 in self.ani_score_self[e1]:
                node2 = get_handle_name(e2)
                nodes_self.add(node2)
                nodes.add(node2)
                if e1 != e2:
                    score, p1, p2 = ani_score_self[e1][e2]
                    p_1_2 = (node1, node2)
                    p_2_1 = (node2, node1)
                    if not p_1_2 in edges and not p_2_1 in edges:
                        edges[p_1_2] = [score]
                    elif p_1_2 in edges:
                        edges[p_1_2].append(score)
                    elif p_2_1 in edges:
                        edges[p_2_1].append(score)
                    else:
                        print('this should never print')

        print(len(nodes), len(edges))

        import networkx as nx
        _edges = []
        for p in edges:
            e1, e2 = p
            scores = edges[p]
            total = 0
            for s in scores:
                total += float(s)
            score = total / len(scores)
            _edges.append([e1, e2, score])

        self.ani_graph = nx.Graph()
        self.ani_graph.add_nodes_from(nodes)
        self.ani_graph.add_weighted_edges_from(_edges)

        return self.ani_graph

    @staticmethod
    def walk(graph, start, distance):
        nodes = set()
        MAGE._walk_collect(graph, start, distance, nodes)
        return nodes

    @staticmethod
    def _walk_collect(graph, pos, distance, collect):
        collect.add(pos)
        for dst, attr in G[pos].items():
            v = attr['weight']
            # print(pos, dst, v, v >= distance)
            if dst not in collect:
                if v >= distance:
                    MAGE._walk_collect(graph, dst, distance, collect)

    def walk_graph(self, centroids_nodes, min_ani=75, max_ani=100):
        data = {}
        for walk_distance in range(min_ani, max_ani):
            centroids = set(centroids_nodes)
            centroid_nodes = {}
            for n in G.nodes:
                if n in centroids:
                    centroid_nodes[n] = {}
            print(len(centroid_nodes))
            visited = set()
            clusters = []
            for node_c in centroid_nodes:
                if node_c not in visited:
                    visits = MAGE.walk(self.ani_graph, node_c, walk_distance)
                    visited |= visits
                    clusters.append(visits)
            data[walk_distance] = {
                'clusters': clusters,
                'visited': visited
            }

        return data

    @staticmethod
    def compute_distance_singletons(data, centroids):

        distance_singletons = {}
        for walk_distance in range(70, 100):
            d = data[walk_distance]
            singletons = set()
            for c in d['clusters']:
                sz = len(c)
                centroids_in_c = centroids | c
                if sz == 1:
                    singletons |= c
                elif centroids_in_c == c:
                    # print(c)
                    pass
            distance_singletons[walk_distance] = singletons
            # if c in centroids
            print(walk_distance, len(d['clusters']), len(d['visited']), len(singletons))

    @staticmethod
    def compute_all_genomes(data):
        all_genomes = set()
        for i in data:
            for c in data[i]['clusters']:
                all_genomes |= c

    def aaaaa(self, all_genomes, metadata_gtdb_gca, metadata_gtdb_gcf):
        genome_to_file = {}
        for k in all_genomes:
            gca = None
            gcf = None
            if k in metadata_gtdb_gca['h_to_genome']:
                gca = metadata_gtdb_gca['h_to_genome'][k]
            if k in metadata_gtdb_gcf['h_to_genome']:
                gcf = metadata_gtdb_gcf['h_to_genome'][k]
            if gca is None and gcf is None:
                # print(k, gca, gcf)
                pass
            elif (not gca is None and gcf is None) or (gca is None and not gcf is None):
                d = None
                if gca is None:
                    d = set(gcf)
                else:
                    d = set(gca)
                if len(d) == 1:
                    pass
                else:
                    print(k, gca, gcf)
            else:
                print(k, gca, gcf)

        return genome_to_file


