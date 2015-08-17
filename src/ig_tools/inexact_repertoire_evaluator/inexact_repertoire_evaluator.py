from collections import defaultdict
from data_structs import ConnectionType
from data_structs import Connection
import sys


class RepertoireEvaluator:
    def __init__(self, repertoires):
        self.repertoires = repertoires
        self.neighbour_clusters = None
        self.connections = defaultdict(dict)
        self.components = []
        self.close_components = []
        self.original_repertoire = None

    # -------------------------------------------------------------------------
    # process neighbour clusters
    # -------------------------------------------------------------------------

    def build_neighbour_connections(self, rep_id1, rep_id2, neighbours):
        for cluster_id1, (cluster2_ids, score) in neighbours.items():
            for cluster_id2 in cluster2_ids:
                if score == 0:
                    continue
                if cluster_id1 in self.connections[rep_id1][rep_id2] and \
                        cluster_id2 in self.connections[rep_id1][rep_id2][cluster_id1]:
                    self.connections[rep_id1][rep_id2][cluster_id1][cluster_id2].type = ConnectionType.close_each
                else:
                    self.connections[rep_id1][rep_id2][cluster_id1][cluster_id2] = \
                        Connection((rep_id1, cluster_id1), (rep_id2, cluster_id2), set(), ConnectionType.close)
                    self.connections[rep_id2][rep_id1][cluster_id2][cluster_id1] = \
                        self.connections[rep_id1][rep_id2][cluster_id1][cluster_id2]

    def add_neighbour_clusters(self, neighbour_clusters):
        self.neighbour_clusters = neighbour_clusters
        for rep_id1 in xrange(len(self.repertoires)):
            for rep_id2 in xrange(len(self.repertoires)):
                self.connections[rep_id1][rep_id2] = defaultdict(dict)

        for rep_id1, item in enumerate(self.neighbour_clusters):
            for rep_id2, neighbours in enumerate(item):
                if not neighbours:
                    continue
                self.build_neighbour_connections(rep_id1, rep_id2, neighbours)

    # -------------------------------------------------------------------------
    # process reads
    # -------------------------------------------------------------------------

    def process_reads_for_reps(self, rep_id1, rep1, rep_id2, rep2):
        for cluster_id1, reads in rep1.cluster_reads.items():
            for read in reads:
                cluster_id2 = rep2.read_clusters[read]
                if cluster_id2 == -1:
                    continue
                if cluster_id1 not in self.connections[rep_id1][rep_id2] or \
                                cluster_id2 not in self.connections[rep_id1][rep_id2][cluster_id1]:
                    self.connections[rep_id1][rep_id2][cluster_id1][cluster_id2] = \
                        Connection((rep_id1, cluster_id1), (rep_id2, cluster_id2), set(), ConnectionType.distant)
                    self.connections[rep_id2][rep_id1][cluster_id2][cluster_id1] = \
                        self.connections[rep_id1][rep_id2][cluster_id1][cluster_id2]
                self.connections[rep_id1][rep_id2][cluster_id1][cluster_id2].reads.add(read)

    def add_reads(self):
        for rep_id1, rep1 in enumerate(self.repertoires):
            for rep_id2, rep2 in enumerate(self.repertoires):
                if rep_id2 <= rep_id1:
                    continue
                if not self.repertoires[rep_id1].cluster_reads or \
                        not self.repertoires[rep_id2].cluster_reads:
                    continue
                self.process_reads_for_reps(rep_id1, rep1, rep_id2, rep2)

    # -------------------------------------------------------------------------
    # build connectivity components
    # -------------------------------------------------------------------------

    def build_connectivity_components_impl(self, close_only):
        components = [None] * len(self.repertoires)
        for rep_id, rep in enumerate(self.repertoires):
            components[rep_id] = dict((cluster_id, -1) for cluster_id in rep.clusters)

        curr_cluster = 0
        for rep_id1, rep1 in enumerate(self.repertoires):
            for cluster_id1 in rep1.clusters:
                if components[rep_id1][cluster_id1] != -1:
                    continue
                curr_cluster += 1
                components[rep_id1][cluster_id1] = curr_cluster
                queue = []
                for rep_id2 in xrange(len(self.repertoires)):
                    if rep_id1 == rep_id2:
                        continue
                    queue = [(rep_id2, cluster_id2) for cluster_id2, edge
                             in self.connections[rep_id1][rep_id2][cluster_id1].items()
                             if (close_only and edge.conn_type != ConnectionType.distant) or not close_only]
                i = 0
                while i != len(queue):
                    curr_rep_id, curr_cluster_id = queue[i]
                    i += 1
                    if components[curr_rep_id][curr_cluster_id] == curr_cluster:
                        continue
                    components[curr_rep_id][curr_cluster_id] = curr_cluster
                    for rep_id2 in xrange(len(self.repertoires)):
                        if curr_rep_id == rep_id2:
                            continue
                        queue += [(rep_id2, cluster_id2) for cluster_id2, edge
                                  in self.connections[curr_rep_id][rep_id2][curr_cluster_id].items()
                                  if ((close_only and edge.conn_type != ConnectionType.distant) or not close_only)
                                  and (rep_id2, cluster_id2) not in queue]
        return components

    def build_connectivity_components(self):
        self.close_components = self.build_connectivity_components_impl(True)
        self.components = self.build_connectivity_components_impl(False)

    # -------------------------------------------------------------------------
    # additional functions
    # -------------------------------------------------------------------------

    def get_neighbour_clusters(self, rep_id1, rep_id2, cluster_no1):
        if cluster_no1 not in self.neighbour_clusters[rep_id1][rep_id2]:
            return None
        return self.neighbour_clusters[rep_id1][rep_id2][cluster_no1][0]

    def is_isolated(self, rep_id, cluster_id, rep_ids=None):
        if not rep_ids:
            rep_ids = range(len(self.repertoires))
        for other_rep_id in rep_ids:
            if other_rep_id == rep_id:
                continue
            if not self.get_neighbour_clusters(rep_id, other_rep_id, cluster_id):
                return True
        return False

    def get_maximum_component_id(self):
        return max(max(v.values() for v in self.components))

    def are_neighbours(self, rep_id1, cluster_id1, rep_id2, cluster_id2):
        return cluster_id2 in self.neighbour_clusters[rep_id1][rep_id2][cluster_id1][0] or \
               cluster_id1 in self.neighbour_clusters[rep_id2][rep_id1][cluster_id2][0]

    def get_non_trivial_components(self):
        non_trivial_components = []
        component_sizes = [0] * (self.get_maximum_component_id() + 1)
        for rep_id, item in enumerate(self.components):
            for cluster_id, component_id in item.items():
                component_sizes[component_id] += 1

        for component_id, size in enumerate(component_sizes):
            if size > 2:
                non_trivial_components.append(component_id)
        return non_trivial_components

    def get_components_max_cluster_size(self):
        components_max_cluster_size = [0] * (self.get_maximum_component_id() + 1)
        for rep_id, item in enumerate(self.components):
            for cluster_id, component_id in item.items():
                components_max_cluster_size[component_id] = max(components_max_cluster_size[component_id],
                                                                self.repertoires[rep_id].clusters[cluster_id].size)
        return components_max_cluster_size
