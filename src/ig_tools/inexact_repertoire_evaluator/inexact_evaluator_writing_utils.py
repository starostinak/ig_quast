import os.path
import dot_utils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import string_utils
from collections import defaultdict
from data_structs import ConnectionType
import sys

GLUING_THRESHOLD = 5

class GluedSingleton:
    def __init__(self, rep_id, cluster_id, conn_type, clusters_count):
        self.rep_id = rep_id
        self.cluster_id = cluster_id
        self.conn_type = conn_type
        self.clusters_count = clusters_count

class OutputFormat:
    def __init__(self, format_single_int, format_single_float, format_table_reps,
            format_mult_int, format_mult_float):
        self.format_single_int = format_single_int
        self.format_single_float = format_single_float
        self.format_table_reps = format_table_reps
        self.format_mult_int = format_mult_int
        self.format_mult_float = format_mult_float

def print_repertoires_header(handler, repertoires):
    for i, rep in enumerate(repertoires):
        handler.write(str(i + 1) + '.' + rep.name + '\n')
    handler.write('\n')

def write_single_metrics(handler, metric, fmt):
    handler.write('General metrics:\n')
    handler.write(fmt.format_single_int.format('#ideal groups', metric.get_ideal_groups()))
    handler.write(fmt.format_single_int.format('#trusted groups', metric.get_trusted_groups()))
    handler.write(fmt.format_single_int.format('#untrusted groups(>' +
                  str(metric.max_errors) + ' errors)', metric.get_untrusted_groups()))
    handler.write(fmt.format_single_int.format('#non-trivial ideal groups',
        metric.get_ideal_groups_nt()))
    handler.write(fmt.format_single_int.format('#non-trivial trusted groups',
        metric.get_trusted_groups_nt()))
    if metric.has_reads:
        handler.write(fmt.format_single_int.format('#non-trivial untrusted groups(>' +
                      str(metric.max_errors) + ' errors)', metric.get_untrusted_groups_nt()))
        if metric.size_cutoff:
            handler.write(fmt.format_single_int.format('#big untrusted groups(>' +
                          str(metric.max_errors) + ' errors, size > ' +
                          str(metric.size_cutoff) + ')', metric.get_big_untrusted_groups()))
    handler.write('\n')

    handler.write('Metrics for repertoires:\n')

    curr_reps = dict((rep_id, rep) for rep_id, rep in ((rep_id, metric.evaluator.repertoires[rep_id]) for rep_id in metric.rep_ids))
    rep_names = [os.path.basename(rep.name) for rep_id, rep in curr_reps.items()]

    handler.write(fmt.format_table_reps.format(*rep_names))
    handler.write(fmt.format_mult_int.format('#clusters',
        *(metric.get_clusters_number(rep_id) for rep_id in curr_reps)))
    handler.write(fmt.format_mult_int.format('#singletons',
        *(metric.get_singletons_number(rep_id) for rep_id in curr_reps)))
    handler.write(fmt.format_mult_int.format('max constructed cluster size',
        *(metric.get_max_cluster_size(rep_id) for rep_id in curr_reps)))
    handler.write(fmt.format_mult_float.format('avg constructed cluster size',
        *(metric.get_avg_cluster_size(rep_id) for rep_id in curr_reps)))
    handler.write(fmt.format_mult_int.format('#isolated clusters',
        *(metric.isolated_clusters[rep_id] for rep_id in curr_reps)))
    handler.write(fmt.format_mult_int.format('#short isolated clusters(<300bp)',
        *(metric.get_short_isolated_clusters(rep_id) for rep_id in curr_reps)))
    handler.write(fmt.format_mult_int.format('min isolated cluster size',
        *(metric.get_min_isolated_cluster_size(rep_id) for rep_id in curr_reps)))
    handler.write(fmt.format_mult_float.format('avg isolated cluster size',
        *(metric.get_avg_isolated_cluster_size(rep_id) for rep_id in curr_reps)))
    handler.write(fmt.format_mult_int.format('max isolated cluster size',
        *(metric.get_max_isolated_cluster_size(rep_id) for rep_id in curr_reps)))
    handler.write(fmt.format_mult_int.format('#trivial isolated clusters(size = 1)',
        *(metric.get_trivial_isolated_clusters(rep_id) for rep_id in curr_reps)))

def write_metrics(filename, all_metrics, ext):
    if filename == 'stdout':
        handler = sys.stdout
    else:
        handler = open(filename, 'w')
    for rep_ids, metrics in sorted(all_metrics.items(), key=lambda (key, value): len(key), reverse=True):
        '''
        for rep_id in rep_ids:
            handler.write(str(rep_id + 1) + '.' + metrics.evaluator.repertoires[rep_id].name + '\n')
        '''
        if ext == 'txt':
            handler.write('\n')
            out_fmt = OutputFormat('{:<45s}{:d}\n', '{:45s}{:.3f\n}',
                ' ' * 45 + '{:s}\t' * len(metrics.rep_ids) + '\n',
                '{:<45s}' + '{:<8d}' * len(metrics.rep_ids) + '\n',
                '{:<45s}' + '{:<8.3f}' * len(metrics.rep_ids) + '\n')
            write_single_metrics(handler, metrics, out_fmt)
            handler.write('\n-------------------------------------------------------------------------------------------\n')
        elif ext == 'csv':
            out_fmt = OutputFormat('{:s}\t{:d}\n', '{:s}\t{:.3f\n}',
                '\t{:s}\t' * len(metrics.rep_ids) + '\n',
                '{:s}\t' + '{:d}\t' * len(metrics.rep_ids) + '\n',
                '{:s}\t' + '{:.3f}\t' * len(metrics.rep_ids) + '\n')
            write_single_metrics(handler, metrics, out_fmt)
        handler.write('\n')
    if filename != 'stdout':
        handler.close()

def get_glued_singletons(evaluator):
    glued_clusters = []
    glued_singletons = []
    for rep_id1, item in evaluator.connections.items():
        for rep_id2, connections in item.items():
            for cluster_id1, adj_list in connections.items():
                adj_singletons = [[] for i in xrange(4)]
                size1 = evaluator.repertoires[rep_id1].clusters[cluster_id1].size
                if size1 == 1:
                    continue
                for cluster_id2, edge in adj_list.items():
                    if edge.conn_type == ConnectionType.close and edge.begin == (rep_id1, cluster_id1):
                        continue
                    if not edge.reads:
                        continue
                    if len(evaluator.connections[rep_id2][rep_id1][cluster_id2]) > 1:
                        continue
                    size2 = evaluator.repertoires[rep_id2].clusters[cluster_id2].size
                    if size2 == 1:
                        adj_singletons[edge.conn_type].append((rep_id2, cluster_id2))
                for conn_type, clusters_list in enumerate(adj_singletons):
                    if len(clusters_list) >= GLUING_THRESHOLD:
                        glued_clusters += clusters_list
                        glued_singletons.append(GluedSingleton(rep_id1, cluster_id1,
                            conn_type, len(clusters_list)))
    return glued_clusters, glued_singletons

def write_cluster_graphs_from(output_prefix, evaluator, non_trivial_components,
                              glued_clusters, glued_singletons, size_cutoff, first, max_count):
    handlers = {}
    component_max_reads_count = evaluator.get_components_max_cluster_size()
    last_ind = first
    opened = 0
    for component_id in non_trivial_components:
        if component_id <= first:
            continue
        if opened == max_count:
            break
        last_ind = component_id
        if component_max_reads_count[component_id] <= size_cutoff:
            continue
        filename = os.path.join(output_prefix, 'component' + str(component_id) + '.dot')
        handlers[component_id] = open(filename, 'w')
        dot_utils.write_dot_header(handlers[component_id])
        opened += 1

    for rep_id, rep in enumerate(evaluator.repertoires):
        for cluster_id, cluster in rep.clusters.items():
            component_id = evaluator.components[rep_id][cluster_id]
            if component_id not in handlers:
                continue
            if (rep_id, cluster_id) in glued_clusters:
                continue
            dot_utils.write_cluster_in_dot(handlers[component_id], rep_id, cluster_id, cluster.size)

    for glued_singleton in glued_singletons:
        component_id = evaluator.components[glued_singleton.rep_id][glued_singleton.cluster_id]
        if component_id not in handlers:
            continue
        dot_utils.write_glued_singleton_cluster_in_dot(handlers[component_id], glued_singleton)
        dot_utils.write_edge_to_glued_singleton_in_dot(handlers[component_id], glued_singleton)

    for rep_id1, item in evaluator.connections.items():
        for rep_id2, connection in item.items():
            if rep_id1 >= rep_id2:
                continue
            for cluster_id1, adj_list in connection.items():
                for cluster_id2, edge in adj_list.items():
                    component_id = evaluator.components[rep_id1][cluster_id1]
                    if component_id not in handlers:
                        continue
                    if (rep_id1, cluster_id1) in glued_clusters or (rep_id2, cluster_id2) in glued_clusters:
                        continue
                    dot_utils.write_edge_in_dot(handlers[component_id], edge)

    for handler in handlers.values():
        dot_utils.write_dot_ending(handler)
        handler.close()

    return last_ind

def write_cluster_graphs(output_prefix, evaluator, size_cutoff):
    non_trivial_components = evaluator.get_non_trivial_components()
    if not non_trivial_components:
        return
    last_ind = -1

    glued_clusters, glued_singletons = [], []
    if len(evaluator.repertoires) == 2:
        glued_clusters, glued_singletons = get_glued_singletons(evaluator)

    while last_ind < max(non_trivial_components):
        last_ind = write_cluster_graphs_from(output_prefix, evaluator,
            non_trivial_components, glued_clusters, glued_singletons, size_cutoff, last_ind, 1000)

def write_single_alignment(handler, evaluator, rep_id1, cluster_id1,
                           rep_id2, cluster_id2, for_group):
    seq1 = evaluator.repertoires[rep_id1].clusters[cluster_id1].seq
    seq2 = evaluator.repertoires[rep_id2].clusters[cluster_id2].seq
    size2 = evaluator.repertoires[rep_id2].clusters[cluster_id2].size
    alignment = string_utils.align_sequences(seq1, seq2)
    handler.write('    ' +
        ((str(rep_id1 + 1) + '.' + str(cluster_id1) + ' - ') if for_group else '') +
        str(rep_id2 + 1) + '.' + str(cluster_id2) +
        ('' if for_group else '(size = ' + str(size2) + ')') + ' shares ' +
        str(len(evaluator.connections[rep_id1][rep_id2][cluster_id1][cluster_id2].reads)) +
        ': score = ' + str(alignment.score) + ', shift = ' + str(alignment.shift) +
        ', errors = ' + str(alignment.errors) + ', strand = ' + str(alignment.strand) + '\n')

def write_groups(filename, groups, evaluator):
    handler = open(filename, 'w')
    print_repertoires_header(handler, evaluator.repertoires)
    for i, group in enumerate(groups):
        handler.write('Group #' + str(i + 1) + '\n')
        handler.write('    Members: ')
        for rep_id, cluster_id in group.members:
            handler.write(str(rep_id + 1) + '.' + str(cluster_id) + '(size = ' +
                str(evaluator.repertoires[rep_id].clusters[cluster_id].size) + ')')
            if (rep_id, cluster_id) != group.members[-1]:
                handler.write(', ')
        handler.write('\n')
        handler.write('    Shared rate: ' + str(group.share_rate) + '\n')
        handler.write('\n')
        for i, (rep_id1, cluster_id1) in enumerate(group.members):
            for rep_id2, cluster_id2 in group.members[i+1:]:
                write_single_alignment(handler, evaluator, rep_id1, cluster_id1, rep_id2, cluster_id2, True)
        handler.write('\n')
    handler.close()

def write_big_isolated_clusters(filename, clusters_list, evaluator):
    handler = open(filename, 'w')
    print_repertoires_header(handler, evaluator.repertoires)
    for rep_id1, cluster_id1 in clusters_list:
        handler.write('Cluster ' + str(rep_id1 + 1) + '.' + str(cluster_id1) +
            '(size = ' + str(evaluator.repertoires[rep_id1].clusters[cluster_id1].size) + '):')
        for rep_id2 in xrange(len(evaluator.repertoires)):
            if rep_id1 == rep_id2:
                continue
            handler.write('\n')
            for cluster_id2, edge in evaluator.connections[rep_id1][rep_id2][cluster_id1].items():
                write_single_alignment(handler, evaluator, rep_id1, cluster_id1, rep_id2, cluster_id2, False)
        handler.write('\n-------------------------------------------------------------------------------------------\n\n')

def print_aln_stats_metrics(handler, aln_stats):
    if not aln_stats:
        handler.write('No components to follow conditions for identity metrics\n')
    else:
        handler.write('Number of components, taken into account for identity metrics is ' +
                      str(len(aln_stats)) + '\n')

        identities = [s.identity for s in aln_stats]
        handler.write('Min percent of identity in component ' +
                      str(min(identities)) + '\n')
        handler.write('Avg percent of identity in component ' +
                      str(sum(identities) / float(len(identities))) + '\n')
        handler.write('Max percent of identity in component ' +
                      str(max(identities)) + '\n')

        nums_mismatch_pos = [s.num_mismatch_pos for s in aln_stats]
        handler.write('Min number of mismatch positions in component ' +
                      str(min(nums_mismatch_pos)) + '\n')
        handler.write('Avg number of mismatch positions in component ' +
                      str(sum(nums_mismatch_pos) / float(len(nums_mismatch_pos))) + '\n')
        handler.write('Max number of mismatch positions in component ' +
                      str(max(nums_mismatch_pos)) + '\n')

        nums_gap_pos = [s.num_gap_pos for s in aln_stats]
        handler.write('Min number of gap positions in component ' +
                      str(min(nums_gap_pos)) + '\n')
        handler.write('Avg number of gap positions in component ' +
                      str(sum(nums_gap_pos) / float(len(nums_gap_pos))) + '\n')
        handler.write('Max number of gap positions in component ' +
                      str(max(nums_gap_pos)) + '\n')

        maxs_mismatches = [s.max_mismatches for s in aln_stats]
        handler.write('Min max number of mismatches between two sequences in component ' +
                      str(min(maxs_mismatches)) + '\n')
        handler.write('Avg max number of mismatches between two sequences in component ' +
                      str(sum(maxs_mismatches) / float(len(maxs_mismatches))) + '\n')
        handler.write('Max max number of mismatches between two sequences in component ' +
                      str(max(maxs_mismatches)) + '\n')

        maxs_gaps = [s.max_gaps for s in aln_stats]
        handler.write('Min max number of gaps between two sequences in component ' +
                      str(min(maxs_gaps)) + '\n')
        handler.write('Avg max number of gaps between two sequences in component ' +
                      str(sum(maxs_gaps) / float(len(maxs_gaps))) + '\n')
        handler.write('Max max number of gaps between two sequences in component ' +
                      str(max(maxs_gaps)) + '\n')

        mins_mismatches = [s.min_mismatches for s in aln_stats]
        handler.write('Min min number of mismatches between two sequences in component ' +
                      str(min(mins_mismatches)) + '\n')
        handler.write('Avg min number of mismatches between two sequences in component ' +
                      str(sum(mins_mismatches) / float(len(mins_mismatches))) + '\n')
        handler.write('Max min number of mismatches between two sequences in component ' +
                      str(max(mins_mismatches)) + '\n')

        mins_gaps = [s.min_gaps for s in aln_stats]
        handler.write('Min min number of gaps between two sequences in component ' +
                      str(min(mins_gaps)) + '\n')
        handler.write('Avg min number of gaps between two sequences in component ' +
                      str(sum(mins_gaps) / float(len(mins_gaps))) + '\n')
        handler.write('Max min number of gaps between two sequences in component ' +
                      str(max(mins_gaps)) + '\n')

def write_component_stats(filename, comp_metrics):
    handler = open(filename, 'w')
    handler.write('Number of components ' + str(len(comp_metrics.components)) + '\n')

    num_of_clusters_in_comp = [v for v in comp_metrics.get_number_of_clusters_in_components().values() if v]
    handler.write('Min number of clusters in component ' +
                  str(min(num_of_clusters_in_comp)) + '\n')
    handler.write('Avg number of clusters in component ' +
                  str(sum(num_of_clusters_in_comp) / float(len(num_of_clusters_in_comp))) + '\n')
    handler.write('Max number of clusters in component ' +
                  str(max(num_of_clusters_in_comp)) + '\n')
    aln_stats = comp_metrics.get_alignment_stats()
    print_aln_stats_metrics(handler, aln_stats)
    handler.write('\n')
    for rep_id in range(len(comp_metrics.evaluator.repertoires)):
        num_clusters_in_comp_for_rep = [v for v in comp_metrics.get_number_of_clusters_in_components(rep_id).values() if v]
        handler.write('For ' + comp_metrics.evaluator.repertoires[rep_id].name + ':\n')
        handler.write('Min number of clusters in component ' +
                 str(min(num_clusters_in_comp_for_rep)) + '\n')
        handler.write('Avg number of clusters in component ' +
                 str(sum(num_clusters_in_comp_for_rep) / float(len(num_clusters_in_comp_for_rep))) + '\n')
        handler.write('Max number of clusters in component ' +
                 str(max(num_clusters_in_comp_for_rep)) + '\n')
        aln_stats = comp_metrics.get_alignment_stats(rep_id)
        print_aln_stats_metrics(handler, aln_stats)
        handler.write('\n')
    handler.close()

def draw_component_sizes_distr(filename, comp_metrics):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pylab
    points = []
    # will fail if more than two repertoires
    for r0, r1 in comp_metrics.components.values():
        if r0 and r1:
            points.append([sum(v.size for v in r0),
                           sum(v.size for v in r1)])
    plt.plot(*zip(*points), marker=',', ls='')
    plt.xlabel(comp_metrics.evaluator.repertoires[0].name)
    plt.ylabel(comp_metrics.evaluator.repertoires[1].name)
    plt.savefig(filename)
    plt.gcf().clear()

