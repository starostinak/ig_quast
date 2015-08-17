from collections import defaultdict 
from itertools import combinations
from string_utils import align_sequences
from data_structs import ConnectionType

RATE_CUTOFF = 0.9
LENGTH_CUTOFF = 300

class Group:
    def __init__(self, members, share_rate, trusted):
        self.members = members
        self.share_rate = share_rate
        self.trusted = trusted

class Metrics:
    def __init__(self, evaluator, size_cutoff, max_errors, rep_ids = None):
        self.evaluator = evaluator

        self.size_cutoff = size_cutoff
        self.max_errors = max_errors
        self.rep_ids = rep_ids

        self.isolated_clusters = [0] * len(evaluator.repertoires)
        self.isolated_cluster_sizes = [[] for i in xrange(len(evaluator.repertoires))]
        self.isolated_cluster_lengths = [[] for i in xrange(len(evaluator.repertoires))]
        self.big_isolated_clusters = []

        self.trusted_groups = [0 for i in xrange(6)]
        self.untrusted_groups = [0 for i in xrange(6)]
        self.trusted_groups_nt = [0 for i in xrange(6)]
        self.untrusted_groups_nt = [0 for i in xrange(6)]
        self.big_untrusted_groups = [0 for i in xrange(6)]

        self.suspected_groups = []
        self.has_reads = True
        for rep_id in self.rep_ids:
            if not self.evaluator.repertoires[rep_id].read_clusters:
                self.has_reads = False

    # -------------------------------------------------------------------------
    # count isolated clusters
    # -------------------------------------------------------------------------

    def add_isolated_cluster(self, rep_id, cluster_id):
        self.isolated_clusters[rep_id] += 1
        self.isolated_cluster_lengths[rep_id].append(\
            self.evaluator.repertoires[rep_id].get_cluster_seq_length(cluster_id))
        self.isolated_cluster_sizes[rep_id].append(\
            self.evaluator.repertoires[rep_id].clusters[cluster_id].size)
        if self.evaluator.repertoires[rep_id].clusters[cluster_id].size > self.size_cutoff:
            self.big_isolated_clusters.append((rep_id, cluster_id))

    def count_isolated_clusters(self):
        for rep_id, rep in enumerate(self.evaluator.repertoires):
            if self.rep_ids and rep_id not in self.rep_ids:
                continue
            for cluster_no in rep.clusters:
                if self.evaluator.is_isolated(rep_id, cluster_no, self.rep_ids):
                    self.add_isolated_cluster(rep_id, cluster_no)

    # -------------------------------------------------------------------------
    # count trusted, untrusted and ideal groups
    # -------------------------------------------------------------------------

    def convert_raw_rate(self, rate):
        return int(rate * 10) - 5

    def add_group(self, group):
        if group.trusted:
            self.trusted_groups[self.convert_raw_rate(group.share_rate)] += 1
            for rep_id, cluster_id in group.members:
                if self.evaluator.repertoires[rep_id].clusters[cluster_id].size != 1:
                    self.trusted_groups_nt[self.convert_raw_rate(group.share_rate)] += 1
                    break
        else:
            self.untrusted_groups[self.convert_raw_rate(group.share_rate)] += 1
            nt = False
            for rep_id, cluster_id in group.members:
                if self.evaluator.repertoires[rep_id].clusters[cluster_id].size != 1 and not nt:
                    self.untrusted_groups_nt[self.convert_raw_rate(group.share_rate)] += 1
                    nt = True
                if self.evaluator.repertoires[rep_id].clusters[cluster_id].size \
                        > self.size_cutoff:
                    self.big_untrusted_groups[self.convert_raw_rate(group.share_rate)] += 1
                    if group.share_rate >= RATE_CUTOFF:
                        self.suspected_groups.append(group)
                    break

    def check_inside_trusted(self, edge):
        rep_id1, cluster_id1 = edge.begin
        rep_id2, cluster_id2 = edge.end
        if self.evaluator.close_components[rep_id1][cluster_id1] != \
                self.evaluator.close_components[rep_id2][cluster_id2]:
            return False
        if edge.conn_type == ConnectionType.distant and \
                align_sequences(self.evaluator.repertoires[rep_id1].clusters[cluster_id1].seq,
                                self.evaluator.repertoires[rep_id2].clusters[cluster_id2].seq,
                                self.max_errors).errors > self.max_errors:
            return False
        return True

    def count_trusted_untrusted_groups(self):
        rep_id1 = self.rep_ids[0]
        raw_groups = []
        for cluster_id1 in self.evaluator.repertoires[rep_id1].clusters:
            trusted = True
            group_reads = self.evaluator.repertoires[rep_id1].cluster_reads[cluster_id1] if self.has_reads else None
            cluster_size1 = self.evaluator.repertoires[rep_id1].clusters[cluster_id1].size
            min_group_cluster_size = cluster_size1
            max_groups_cluster_size = cluster_size1
            group_members = [(rep_id1, cluster_id1)]
            for rep_id2 in self.rep_ids[1:]:
                for cluster_id2, edge in self.evaluator.connections[rep_id1][rep_id2][cluster_id1].items():
                    cluster_size2 = self.evaluator.repertoires[rep_id2].clusters[cluster_id2].size
                    shared = float(len(edge.reads))
                    if not self.has_reads:
                        shared = float(min(min_group_cluster_size, cluster_size2))
                    share_rate = min(shared / cluster_size2, shared / max_groups_cluster_size)
                    if share_rate > 0.5:
                        if self.has_reads:
                            group_reads &= edge.reads
                        else:
                            min_group_cluster_size = shared
                        group_members.append((rep_id2, cluster_id2))
                        max_groups_cluster_size = max(max_groups_cluster_size, cluster_size2)
                        if not self.check_inside_trusted(edge):
                            trusted = False
                if self.has_reads and float(len(group_reads)) / max_groups_cluster_size <= 0.5:
                    group_members = []
                    break
            if len(group_members) == len(self.rep_ids):
                if self.has_reads:
                    raw_groups.append(Group(group_members, float(len(group_reads)) / max_groups_cluster_size, trusted))
                else:
                    raw_groups.append(Group(group_members,
                                            float(min_group_cluster_size) / max_groups_cluster_size, trusted))

        for group in raw_groups:
            self.add_group(group)

    # -------------------------------------------------------------------------
    # small additional functions
    # -------------------------------------------------------------------------

    def get_ideal_groups(self):
        return self.trusted_groups[-1]

    def get_trusted_groups(self):
        return sum(self.trusted_groups[int(RATE_CUTOFF * 10) - 5:])

    def get_untrusted_groups(self):
        return sum(self.untrusted_groups[int(RATE_CUTOFF * 10) - 5:])

    def get_ideal_groups_nt(self):
        return self.trusted_groups_nt[-1]

    def get_trusted_groups_nt(self):
        return sum(self.trusted_groups_nt[int(RATE_CUTOFF * 10) - 5:])

    def get_untrusted_groups_nt(self):
        return sum(self.untrusted_groups_nt[int(RATE_CUTOFF * 10) - 5:])

    def get_big_untrusted_groups(self):
        return sum(self.big_untrusted_groups[int(RATE_CUTOFF * 10) - 5:])

    def get_short_isolated_clusters(self, rep_id):
        return len([l for l in self.isolated_cluster_lengths[rep_id] if l < LENGTH_CUTOFF])

    def get_clusters_number(self, rep_id):
        rep = self.evaluator.repertoires[rep_id]
        return len(rep.clusters)

    def get_singletons_number(self, rep_id):
        rep = self.evaluator.repertoires[rep_id]
        return len([c for c in rep.clusters.values() if c.size == 1])

    def get_max_cluster_size(self, rep_id):
        rep = self.evaluator.repertoires[rep_id]
        return max(c.size for c in rep.clusters.values())

    def get_avg_cluster_size(self, rep_id):
        rep = self.evaluator.repertoires[rep_id]
        return float(sum(c.size for c in rep.clusters.values())) / len(rep.clusters.values())

    def get_min_isolated_cluster_size(self, rep_id):
        if not self.isolated_cluster_sizes[rep_id]:
            return 0
        return min(self.isolated_cluster_sizes[rep_id])

    def get_avg_isolated_cluster_size(self, rep_id):
        if not self.isolated_cluster_sizes[rep_id]:
            return 0
        return float(sum(self.isolated_cluster_sizes[rep_id])) / \
               len(self.isolated_cluster_sizes[rep_id])

    def get_max_isolated_cluster_size(self, rep_id):
        if not self.isolated_cluster_sizes[rep_id]:
            return 0
        return max(self.isolated_cluster_sizes[rep_id])

    def get_trivial_isolated_clusters(self, rep_id):
        return self.isolated_cluster_sizes[rep_id].count(1)

    '''
    def calculate_component_stats(self, log):
        components = dict() # {component_id -> [[sizes of clusters for each repertoire]]}
        for rep_id, component in enumerate(self.evaluator.components):
            rep = self.evaluator.repertoires[rep_id]
            for cluster_id, component_id in component.items():
                if component_id not in components:
                    components[component_id] = [[] for i in range(len(self.evaluator.repertoires))]
                components[component_id][rep_id].append(rep.clusters[cluster_id].size)
        log.info('Number of components ' + str(len(components)))
        num_of_clusters_in_comp = [sum(len(rep_clusters) for rep_clusters in comp_clusters) for comp_clusters in components.values()]
        points = []
        for r0, r1 in components.values():
            if r0 and r1:
                points.append([sum(r0), sum(r1)])

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import pylab
        plt.plot(*zip(*points), marker='*', label=[self.evaluator.repertoires[0].name, self.evaluator.repertoires[1].name], ls='')
        plt.xlabel(self.evaluator.repertoires[0].name)
        plt.ylabel(self.evaluator.repertoires[1].name)
        pylab.legend()
        plt.savefig('sum_of_sizes.png')

        log.info('Min number of clusters in component ' + str(min(num_of_clusters_in_comp)))
        log.info('Avg number of clusters in component ' + str(sum(num_of_clusters_in_comp) / float(len(num_of_clusters_in_comp))))
        log.info('Max number of clusters in component ' + str(max(num_of_clusters_in_comp)))
        for rep_id in range(len(self.evaluator.repertoires)):
            log.info('For ' + self.evaluator.repertoires[rep_id].name + ':')
            log.info('Min number of clusters in component ' + str(min(len(a[rep_id]) for a in components.values() if a[rep_id])))
            log.info('Avg number of clusters in component ' + str(sum([len(a[rep_id]) for a in components.values() if a[rep_id]]) / float(len([c for c, val in components.items() if val[rep_id]]))))
            log.info('Max number of clusters in component ' + str(max(len(a[rep_id]) for a in components.values() if a[rep_id])))
            log.info('')
    ''' 

    def evaluate(self):
        self.count_isolated_clusters()
        self.count_trusted_untrusted_groups()

class ClusterSizeId:
    def __init__(self, cluster_id, size):
        self.cluster_id = cluster_id
        self.size = size

class AlignmentStats:
    def __init__(self, identity, num_mismatch_pos, num_gap_pos, max_mismatches, max_gaps, min_mismatches, min_gaps):
        self.identity = identity
        self.num_mismatch_pos = num_mismatch_pos
        self.num_gap_pos = num_gap_pos
        self.max_mismatches = max_mismatches
        self.max_gaps = max_gaps
        self.min_mismatches = min_mismatches
        self.min_gaps = min_gaps

class ComponentMetrics:
    def __init__(self, evaluator, size_cutoff = 10):
        self.evaluator = evaluator
        self.size_cutoff = size_cutoff
        self.components = dict() # {component_id -> [[ClusterSizeIds], [..]]

    def get_number_of_clusters_in_components(self, rep_id = -1):
        if rep_id < 0:
            return dict((comp_id, sum(len(rep_clusters) for rep_clusters in comp_clusters))
                    for comp_id, comp_clusters in self.components.items())
        else:
            return dict((comp_id, len(comp_clusters[rep_id])) for comp_id, comp_clusters in self.components.items())

    def compute_aln_stats(self, alignment):
        aln_length = alignment.get_alignment_length()
        start_positions = [next(i for i, c in enumerate(s)if c != '-') for s in alignment]
        end_positions = [next(len(s) - i - 1 for i, c in enumerate(reversed(s)) if c != '-') for s in alignment]
        identical_count = 0
        distances = {} # {(seq1_id, seq2_id) -> (mismatches_num, gaps_num)}
        mismatch_positions = [False] * aln_length
        gap_positions = [False] * aln_length
        identical_positions = [True] * aln_length
        for i in range(len(alignment)):
            for j in range(i + 1, len(alignment)):
                distances[(i, j)] = [0, 0]
                for pos in range(aln_length):
                    if start_positions[i] > pos or end_positions[i] < pos or \
                            start_positions[i] > pos or end_positions[i] < pos:
                        continue
                    if alignment[i][pos] == alignment[j][pos]:
                        continue
                    else:
                        identical_positions[pos] = False
                        if alignment[i][pos] == '-' or alignment[j][pos] == '-':
                            gap_positions[pos] = True
                            distances[(i, j)][1] += 1
                        else:
                            mismatch_positions[pos] = True
                            distances[(i, j)][0] += 1
        num_mismatch_pos = sum(mismatch_positions)
        num_gap_pos = sum(gap_positions)
        num_identical_pos = sum(identical_positions)
        aln_stats = AlignmentStats(num_identical_pos / float(aln_length),
                                   num_mismatch_pos,
                                   num_gap_pos,
                                   max(a[0] for a in distances.values()),
                                   max(a[1] for a in distances.values()),
                                   min(a[0] for a in distances.values()),
                                   min(a[1] for a in distances.values()))
        return aln_stats

    def get_alignment_stats_for_component(self, component_id, rep_id_to_process = -1):
        from Bio.Align.Applications import MuscleCommandline
        from Bio.Align.Applications import ClustalwCommandline
        from Bio import AlignIO
        from StringIO import StringIO
        muscle_cline = MuscleCommandline(clwstrict=True, diags=True, maxiters=2)
        to_align_str = ''
        for rep_id, clusters_list in enumerate(self.components[component_id]):
            if rep_id_to_process >= 0 and rep_id_to_process != rep_id:
                continue
            for value in clusters_list:
                cluster_id = value.cluster_id
                to_align_str += '>cluster___' + str(cluster_id) + '\n'
                to_align_str += self.evaluator.repertoires[rep_id].clusters[cluster_id].seq + '\n'
        handle = StringIO()
        handle.write(to_align_str)
        stdout, stderr = muscle_cline(stdin=handle.getvalue())
        alignment = AlignIO.read(StringIO(stdout), 'clustal')
        aln_stats = self.compute_aln_stats(alignment)
        return aln_stats

    def get_alignment_stats(self, rep_id = -1):
        alignments_stats = []
        num_of_clusters = self.get_number_of_clusters_in_components(rep_id)
        for component_id, rep_clusters in self.components.items():
            if num_of_clusters[component_id] < self.size_cutoff:
                continue
            alignments_stats.append(self.get_alignment_stats_for_component(component_id, rep_id))
        return alignments_stats

    def evaluate(self):
        for rep_id, component in enumerate(self.evaluator.components):
            rep = self.evaluator.repertoires[rep_id]
            for cluster_id, component_id in component.items():
                if component_id not in self.components:
                    self.components[component_id] = [[] for i in range(len(self.evaluator.repertoires))]
                self.components[component_id][rep_id].append(ClusterSizeId(cluster_id, rep.clusters[cluster_id].size))

def evaluate_metrics(evaluator, size_cutoff, max_errors):
    metrics = {}
    for i in range(2, len(evaluator.repertoires) + 1):
        for comb in combinations(range(len(evaluator.repertoires)), i):
            metrics[comb] = Metrics(evaluator, size_cutoff, max_errors, comb)
    
    for metric in metrics.values():
        metric.evaluate()

    return metrics
