#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import getopt
import os
import logging
import shutil
import datetime
from time import gmtime, strftime

import init
from init import AbnormalFinishMsg

import drawing_utils
import clusters_utils
import files_utils
import inexact_metrics_utils
from inexact_repertoire_evaluator import RepertoireEvaluator
import reading_utils
import inexact_evaluator_writing_utils
import somatic_mutation_search

class Params:
    max_repertoires_number = 2

    constructed_clusters = None
    original_clusters = ["", ""]
    assembled_barcodes = ["", ""]

    draw_hist = True
    output_basename = ""
    output_dir = ""
    stats_dir = "stats"

    deep_analysis = False
    blast_output = ""
    igblast_output = ""
    min_cluster_size = 10

    file_metrics = {}
    file_metrics_csv = {}

    test_dataset = False
    log = ""

    size_cutoff = 20
    tau = 3
    # max_mismatches = 4
    # max_gaps = 0

    threads_num = 4

    original_sizes_list = []
    constructed_sizes_list = []
    original_sizes_file = ""
    constructed_sizes_file = []
    draw_cluster_graphs = False

    evaluator = None
    inexact_metrics = dict()

    # temporary 
    # neighbour_clusters_file = ""
    matches_files = []
    rerun = False

class BaseOptions:
    long_options = ("skip-drawing test-original test-general test-mult-cmp help adv-analysis adv-output tau= isol-min-size= adv-min-size= Rc= Rr= Bc= Br= out= threads-num= rerun".split() +
                    ["c" + str(i) + "=" for i in
                        xrange(1, Params.max_repertoires_number + 1)] +
                    ["r" + str(i) + "=" for i in
                        xrange(1, Params.max_repertoires_number + 1)])
    # temporary option -n for debugging
    short_options = "R:n:"

def usage(log):
    log.info("./ig_quast.py [options]")
    log.info("\nRepertoire options:")
    log.info("  --Rc\t\t\t<filename>\tCLUSTERS file with original (ideal) repertoire")
    log.info("  --Rr\t\t\t<filename>\tRCM file with original (ideal) repertoire")
    log.info("  --Bc\t\t\t<filename>\tCLUSTERS file with repertoire corresponding to assembled barcodes")
    log.info("  --Br\t\t\t<filename>\tRCM file with for repertoire corresponding to assembled barcodes")
    log.info("  --c<#>\t\t<filename>\tCLUSTERS file with constructed repertoire number <#> (<#> = 1,2)")
    log.info("  --r<#>\t\t<filename>\tRCM file with constructed repertoire number <#> (<#> = 1,2)")

    log.info("\nBasic options:")
    log.info("  --test-original\t\t\t\truns test dataset")
    log.info("  --test-general\t\t\t\truns test dataset")
    log.info("  --test-mult-cmp\t\t\t\truns test dataset")
    log.info("  --out\t\t\t<filename>\toutput directory")
    log.info("  --help\t\t\t\tprints help")

    log.info("\n\"Multiple repertoires comparison\" options:")
    log.info("  --tau\t\t<int>\t\tmaximal allowed distance between clusters [default: %i]" % Params.tau)
    log.info("  --threads-num\t\t<int>\t\tnumber of threads to use [default: %i]" % Params.threads_num)
    log.info("  --isol-min-size\t\t<int>\t\tsize cutoff for isolated clusters comparison and drawing graphs [default: %i]" % Params.size_cutoff)
    log.info("  --adv-output\t\t<int>\t\toutput cluster graphs")

    log.info("\n\"Single repertoire analysis\" options:")
    log.info("  --adv-analysis\t\t\tenables advanced analysis of mutated groups")
    log.info("  --adv-min-size\t<int>\t\tminimal size of clusters that will be analyzed by advanced analysis of mutated groups [default: 10]")

def ParseCommandLine(params, log):
    constructed_clusters_raw = [["", ""] for i 
        in range(params.max_repertoires_number + 1)]
    # preparing command line arguments
    options, tmp = getopt.gnu_getopt(sys.argv, BaseOptions.short_options, BaseOptions.long_options)
    # reading options
    for opt, arg in options:
        if opt == '--Rc':
            params.original_clusters[0] = arg
        elif opt == '--Rr':
            params.original_clusters[1] = arg
        elif opt == '--Bc':
            params.assembled_barcodes[0] = arg
        elif opt == '--Br':
            params.assembled_barcodes[1] = arg
        elif opt == '--skip-drawing':
            params.draw_hist = False
        elif opt == '--test-original':
            constructed_clusters_raw[1] = ["ig_test_dataset/constructed.clusters.fa", "ig_test_dataset/constructed.rcm"]
            params.original_clusters = ["ig_test_dataset/ideal.clusters.fa", "ig_test_dataset/ideal.rcm"]
            params.test_dataset = True
        elif opt == '--test-general':
            constructed_clusters_raw[1] = ["ig_test_dataset/constructed.clusters.fa", "ig_test_dataset/constructed.rcm"]
            params.min_cluster_size = 10
            params.deep_analysis = True
            params.test_dataset = True
        elif opt == '--test-mult-cmp':
            constructed_clusters_raw[1] = ["ig_test_dataset/constructed.clusters.fa", "ig_test_dataset/constructed.rcm"]
            constructed_clusters_raw[2] = ["ig_test_dataset/ideal.clusters.fa", "ig_test_dataset/ideal.rcm"]
            params.test_dataset = True
        elif opt == '--out':
            params.output_dir = arg
        elif opt == '--help':
            usage(log)
            sys.exit(0)
        elif opt == '--tau':
            params.tau = int(arg)
        elif opt == '--threads-num':
            params.threads_num = int(arg)
        elif opt == '--isol-min-size':
            params.size_cutoff = int(arg)
        elif opt == '-n':
            params.matches_files = arg
        elif opt[:-1] == '--c':
            num = int(opt[-1])
            constructed_clusters_raw[num][0] = os.path.abspath(arg)
        elif opt[:-1] == '--r':
            num = int(opt[-1])
            constructed_clusters_raw[num][1] = os.path.abspath(arg)
        elif opt == '--adv-analysis':
            params.deep_analysis = True
        elif opt == '--adv-min-size':
            params.min_cluster_size = int(arg)
        elif opt == '--adv-output':
            params.draw_cluster_graphs = True
        elif opt == '--rerun':
            params.rerun = True
            
    if len(options) == 0 and len(tmp) == 1 :
        usage(log)
        sys.exit(1)

    params.constructed_clusters = [(clusters, rcm) for clusters, rcm 
        in constructed_clusters_raw if clusters or rcm]

    if len(params.constructed_clusters) == 0:
        log.info("ERROR: Invalid command line, you must provide at least one CLUSTERS or RCM file with constructed clusters")
        sys.exit(1)

def PrepareOutputDir(base_dir, params):
    if params.output_dir == "":
        base_output_dir = os.path.join(base_dir, "ig_quast_results")
        if not os.path.exists(base_output_dir):
            os.makedirs(base_output_dir)
        date_time = strftime("%Y-%m-%d_%H-%M-%S", gmtime())
        cur_output_dir = os.path.join(base_output_dir, date_time)
        params.output_dir = cur_output_dir
        os.makedirs(cur_output_dir)
    params.output_basename = os.path.join(params.output_dir, "results")
    params.stats_dir = os.path.join(params.output_dir, "stats")
    if not os.path.exists(params.stats_dir):
        os.makedirs(params.stats_dir)

def CreateRepertoireSubdirs(params, log):
    if len(params.evaluator.repertoires) == 1 or params.original_clusters == ["", ""]:
        return
    for rep in params.evaluator.repertoires:
        curr_dir = os.path.join(params.output_dir, rep.name)
        if not os.path.exists(curr_dir):
            os.makedirs(curr_dir)

# --------------------------------------------------------------------------------

def RunExactEvaluator(params, log):
    for rep in params.evaluator.repertoires:
        clusters, rcm = rep.clusters_filename, rep.rcm_filename
        if len(params.constructed_clusters) == 1:
            out_name = os.path.join(params.output_dir, "metrics")
        else:
            out_name = os.path.join(params.output_dir, rep.name, "metrics")
        command_line = init.PathToBins.run_repertoire_evaluator_tool + " " + \
            params.original_clusters[1] + " " + params.original_clusters[0] + " " + \
            rcm + " " + clusters + " " + out_name
        error_code = os.system(command_line + " 2>&1 | tee -a " + params.log)

        if error_code != 0:
            AbnormalFinishMsg(log, "repertoire_evaluator")
            sys.exit(1)

        params.file_metrics[rep.name] = os.path.join(out_name + ".txt")
        params.file_metrics_csv[rep.name] = os.path.join(out_name + ".csv")

        if os.path.exists(params.file_metrics[rep.name]) and os.path.exists(params.file_metrics_csv[rep.name]):
            log.info("\n* Computed metrics for " + clusters + " in TXT were written to " + params.file_metrics[rep.name])
            log.info("* Computed metrics for " + clusters + " in CSV were written to " + params.file_metrics_csv[rep.name])
        else:
            log.info("ERROR: Computed metrics were not found")
            sys.exit(1)

    if len(params.file_metrics) > 1:
        MergeMetricsIntoSingleFile(params, log)
        params.file_metrics['general'] = os.path.join(params.output_dir, "metrics.txt")
        params.file_metrics_csv['general'] = os.path.join(params.output_dir, "metrics.csv")
        if os.path.exists(params.file_metrics['general']) and os.path.exists(params.file_metrics_csv['general']):
            log.info("\n* Merged computed metrics in TXT were written to " + params.file_metrics['general'])
            log.info("* Merged computed metrics in CSV were written to " + params.file_metrics_csv['general'])
        else:
            log.info("ERROR: Computed metrics were not found")
            sys.exit(1)
    else:
        params.file_metrics['general'] = params.file_metrics.values()[0]
        params.file_metrics_csv['general'] = params.file_metrics_csv.values()[0]

def RunGeneralEvaluator(params, log):
    # temporary solution
    rcm_file = params.constructed_clusters[0][1]
    ids, rcm = clusters_utils.ParseRCM(rcm_file)

    result = list()
    sizes = ConvertRcmToClusterSizes(rcm)

    result.append("#constructed clusters")
    result.append(len(rcm))

    result.append("#constructed singletons")
    result.append(len([s for s in rcm if len(rcm[s]) == 1]))

    result.append("max constructed cluster")
    result.append(drawing_utils.MaxValue(sizes))

    result.append("avg cluster size")
    result.append(drawing_utils.AverageValue(sizes))

    result.append("# clusters (>=10)")
    result.append(drawing_utils.NumValuesGreaterThresholds(sizes, 10))

    result.append("# clusters (>=50)")
    result.append(drawing_utils.NumValuesGreaterThresholds(sizes, 50))

    result.append("# clusters (>=100)")
    result.append(drawing_utils.NumValuesGreaterThresholds(sizes, 100))

    result.append("# clusters (>=500)")
    result.append(drawing_utils.NumValuesGreaterThresholds(sizes, 500))

    result.append("# clusters (>=1000)")
    result.append(drawing_utils.NumValuesGreaterThresholds(sizes, 1000))

    clusters_utils.PrintMetricsDict(result, log)
    params.file_metrics['general'] = os.path.join(params.output_dir, "metrics.txt")
    params.file_metrics_csv['general'] = os.path.join(params.output_dir, "metrics.csv")
    params.file_metrics[params.evaluator.repertoires[0].name] = params.file_metrics['general']
    params.file_metrics_csv[params.evaluator.repertoires[0].name] = params.file_metrics_csv['general']
    clusters_utils.PrintMetricsDictToFile(result, params.file_metrics['general'])
    clusters_utils.PrintMetricsDictToCSV(result, params.file_metrics_csv['general'])
    if os.path.exists(params.file_metrics['general']) and os.path.exists(params.file_metrics_csv['general']):
        log.info("\n* Computed metrics in TXT were written to " + params.file_metrics['general'])
        log.info("* Computed metrics in CSV were written to " + params.file_metrics_csv['general'])
    else:
        log.info("ERROR: Computed metrics were not found")
        sys.exit(1)

'''
def RunRepertoiresComparer(log, params):
    log.info(".. Building neighbour clusters")
    neighbour_clusters_output = os.path.join(params.stats_dir, "neighbour_clusters.txt")
    if params.neighbour_clusters_file != "":
        neighbour_clusters_output = params.neighbour_clusters_file
    else:
        command_line = init.PathToBins.run_repertoire_comparer_tool + \
                ' --max_mismatches ' + str(params.max_mismatches) + ' --max_gaps ' + str(params.max_gaps) + \
                ' --threads ' + str(params.threads_num) + ' --out ' + neighbour_clusters_output + ' ' + \
                ' '.join(rep.clusters_filename for rep in params.evaluator.repertoires)
        exit_code = os.system(command_line + ' 2>&1 | tee -a ' + params.log)
        if exit_code != 0:
            AbnormalFinishMsg(log, "repertoire_comparer")
            sys.exit(-1)
    if not os.path.exists(neighbour_clusters_output):
        log.info('ERROR: cannot find neighbour clusters file') 
        sys.exit(-1)

    log.info("* Neighbour clusters were written to " + neighbour_clusters_output)
    return reading_utils.read_neighbour_clusters(neighbour_clusters_output, params.evaluator.repertoires)
'''

def RunIgMatcher(log, params):
    log.info(".. Building neighbour clusters")
    matches_files = [os.path.join(params.stats_dir, params.evaluator.repertoires[0].name), 
                     os.path.join(params.stats_dir, params.evaluator.repertoires[1].name)]
    aln_filenames = []
    cluster_num_to_ids = []
    for rep in params.evaluator.repertoires:
        aln_filenames.append(os.path.join(params.stats_dir, rep.name))
        command_line = init.PathToBins.run_ig_kplus_vj_finder_tool + \
                ' -i ' + rep.clusters_filename + \
                ' -o ' + aln_filenames[-1] + '.vj.fa' + \
                ' -b ' + os.path.join(params.stats_dir, rep.name + '.bad.fa') + ' -S'
        if not params.rerun:
            exit_code = os.system(command_line + ' 2>&1 | tee -a ' + params.log)
            if exit_code != 0:
                AbnormalFinishMsg(log, 'ig_kplus_vj_finder')
                sys.exit(-1)
        if not os.path.exists(aln_filenames[-1] + '.vj.fa'):
            log.info('ERROR: cannot find alignment file ' + aln_filenames[-1] + '.vj.fa')
        command_line = init.PathToBins.run_ig_trie_compressor_tool + \
                ' -i ' + aln_filenames[-1] + '.vj.fa' + \
                ' -o ' + aln_filenames[-1] + '.fa'
        if not params.rerun:
            exit_code = os.system(command_line + ' 2>&1 | tee -a ' + params.log)
            if exit_code != 0:
                AbnormalFinishMsg(log, 'ig_trie_compressor')
                sys.exit(-1)
        if not os.path.exists(aln_filenames[-1] + '.fa'):
            log.info('ERROR: cannot find alignment file ' + aln_filenames[-1] + '.fa')
            sys.exit(-1)
        aln_filenames[-1] += '.fa'
        cluster_num_to_ids.append(reading_utils.read_cluster_numbers_to_ids(aln_filenames[-1]))

    command_line = init.PathToBins.run_ig_matcher_tool + \
            ' -i ' + aln_filenames[0] + \
            ' -I ' + aln_filenames[1] + \
            ' --tau ' + str(params.tau) + ' --threads ' + str(params.threads_num) + \
            ' -o ' + matches_files[0] + \
            ' -O ' + matches_files[1]

    if not params.rerun:
        exit_code = os.system(command_line + ' 2>&1 | tee -a ' + params.log)
        if exit_code != 0:
            AbnormalFinishMsg(log, "ig_matcher")
            sys.exit(-1)
    if any(not os.path.exists(filename) for filename in matches_files):
        log.info('ERROR: cannot find neighbour clusters file') 
        sys.exit(-1)

    log.info("* Neighbour clusters were written to " + str(matches_files))
    return reading_utils.read_neighbour_clusters(matches_files, params.evaluator.repertoires, cluster_num_to_ids)

def RunInexactEvaluator(params, log):
    neighbour_clusters = RunIgMatcher(log, params)
    log.info(".. Analizing neighbour clusters")
    params.evaluator.add_neighbour_clusters(neighbour_clusters)
    log.info(".. Analyzing reads")
    params.evaluator.add_reads()
    log.info(".. Building connectivity components")
    params.evaluator.build_connectivity_components()
    log.info(".. Calculating metrics")
    params.inexact_metrics = inexact_metrics_utils.evaluate_metrics(params.evaluator, 
        params.size_cutoff, params.tau)
    params.file_metrics['general'] = os.path.join(params.output_dir, "metrics.txt")
    params.file_metrics_csv['general'] = os.path.join(params.output_dir, "metrics.csv")
    inexact_evaluator_writing_utils.write_metrics('stdout', params.inexact_metrics, 'txt')
    inexact_evaluator_writing_utils.write_metrics(params.file_metrics['general'], params.inexact_metrics, 'txt')
    log.info("* Computed metrics in TXT were written to " + params.file_metrics['general'])
    inexact_evaluator_writing_utils.write_metrics(params.file_metrics_csv['general'], params.inexact_metrics, 'csv')
    log.info("* Computed metrics in CSV were written to " + params.file_metrics_csv['general'])

    if params.draw_cluster_graphs:
        log.info("\n==== Writing advanced output")
        graphs_dir = os.path.join(params.output_dir, "cluster_graphs")
        if os.path.exists(graphs_dir):
            shutil.rmtree(graphs_dir)
        os.makedirs(graphs_dir)
        inexact_evaluator_writing_utils.write_cluster_graphs(graphs_dir, params.evaluator, params.size_cutoff)
        log.info("* Cluster graphs were written to " + graphs_dir + " directory")

    main_metrics = sorted(params.inexact_metrics.items(), 
                          key=lambda (key, value): len(key), reverse=True)[0][1]
    suspected_groups_file = os.path.join(params.stats_dir, "big_untrusted_groups.txt")
    inexact_evaluator_writing_utils.write_groups(suspected_groups_file, 
        main_metrics.suspected_groups, params.evaluator)
    log.info("* Big untrusted groups were written to " + suspected_groups_file)

    big_isolated_clusters_file = os.path.join(params.stats_dir, "big_isolated_clusters.txt")
    inexact_evaluator_writing_utils.write_big_isolated_clusters(big_isolated_clusters_file, 
        main_metrics.big_isolated_clusters, params.evaluator)
    log.info("* Big isolated clusters were written to " + suspected_groups_file)

    main_metrics.get_barcode_metrics(params.stats_dir)

    '''
    component_metrics = inexact_metrics_utils.ComponentMetrics(params.evaluator, params.size_cutoff)
    component_metrics.evaluate()
    inexact_evaluator_writing_utils.write_component_stats(
        os.path.join(params.output_dir, 'component_metrics.txt'), component_metrics)
    inexact_evaluator_writing_utils.draw_component_sizes_distr(
        os.path.join(params.output_dir, 'comp_cluster_sizes.png'), component_metrics)
    '''


# -----------------------------------------------------------------------------

def ConvertRcmToClusterSizes(rcm):
    sizes = list()
    for key in rcm:
        sizes.append(len(rcm[key]))
    return sizes

def ReadRepertoires(params, log):
    repertoires, read_names = reading_utils.read_repertoires(params.constructed_clusters, params.assembled_barcodes)
    params.evaluator = RepertoireEvaluator(repertoires)
    log.info('\n==== Repertoires:')
    for rep in params.evaluator.repertoires:
        log.info(rep.clusters_filename + ' ==> ' + rep.name)
    log.info('')


def MergeMetricsIntoSingleFile(params, log):
    fields = []
    metrics = {}
    for rep_name, file_metrics in params.file_metrics_csv.items():
        metrics[rep_name] = {}
        handler = open(file_metrics)
        for l in handler:
            name, value = l.strip().split('\t')
            metrics[rep_name][name] = value
            if name not in fields:
                fields.append(name)
    clusters_utils.PrintMultipleMetricsDictToFile(metrics, fields, os.path.join(params.output_dir, "metrics.txt")) 
    clusters_utils.PrintMultipleMetricsDictToCSV(metrics, fields, os.path.join(params.output_dir, "metrics.csv")) 

def ReadIdentityMetrics(filename):
    handler = open(filename)
    percentages, cluster_lengths = [], []
    for l in handler:
        perc, cluster_len = l.strip().split('\t')
        percentages.append(float(perc))
        cluster_lengths.append(int(cluster_len))
    return percentages, cluster_lengths

# -------------------------------------------------------------------------

def DrawIdentityMetrics(params, log):
    for rep_name, file_metrics in params.file_metrics.items():
        if rep_name == 'general':
            continue
        file_size_perc = ".".join(file_metrics.split(".")[:-1]) + ".size.perc"
        percentages, cluster_lengths = ReadIdentityMetrics(file_size_perc)
        out_basename = os.path.basename(file_size_perc).split(".")[0].split("_", 1)
        hist_name, plot_name = '', ''
        if len(params.evaluator.repertoires) == 1:
            hist_name = os.path.join(params.output_dir, "identity_perc_distribution.png")
            plot_name = os.path.join(params.output_dir, "identity_perc_length.png")
        else:
            hist_name = os.path.join(params.output_dir, rep_name, "identity_perc_distribution.png")
            plot_name = os.path.join(params.output_dir, rep_name, "identity_perc_length.png")
        drawing_utils.DrawIdentityPercentageDistribution(percentages, hist_name)
        if os.path.exists(hist_name):
            log.info("* Histogram of identity percentage distribution were written to " + hist_name)
        else:
            log.info("Histogram of identity percentage distribution")

        drawing_utils.DrawIdentityToLengthDistribution(percentages, cluster_lengths, plot_name)
        if os.path.exists(plot_name):
            log.info("* Plot of identity percentage to sequence length distribution were written to " + plot_name)
        else:
            log.info("Plot of identity percentage to sequence length distribution")

def DrawClusterSizesHist(params, log):
    log.info("\n==== Drawing comparative histograms for all constructed cluster sizes distribution")
    all_cluster_sizes = []
    rep_names = []
    hist_all = os.path.join(params.output_dir, "all_cluster_sizes.png")
    for rep in params.evaluator.repertoires:
        all_cluster_sizes.append([cluster.size for cluster in rep.clusters.values()])
        rep_names.append(rep.name)

    if params.evaluator.original_repertoire:
        all_cluster_sizes.append([cluster.size for cluster in params.evaluator.original_repertoire.clusters.values()])
        rep_names.append("Original repertoire")

    drawing_utils.DrawAnyClusterSizesHist(all_cluster_sizes, rep_names, hist_all)
    if os.path.exists(hist_all):
        log.info("* Histogram of all clusters sizes distribution was written to " + hist_all)
    else:
        log.info("ERROR: files with histograms of all clusters sizes distribution was not found")
        sys.exit(1)

    hist_all_nt = os.path.join(params.output_dir, "nt_cluster_sizes.png")
    drawing_utils.DrawAnyClusterSizesHist([[sz for sz in sizes if sz > 1] for sizes in 
                                            all_cluster_sizes], 
                                          rep_names, hist_all_nt)
    if os.path.exists(hist_all_nt):
        log.info("* Histogram of all non-trivial clusters sizes distribution was written to " + hist_all_nt)
    else:
        log.info("Histogram of all non-trivial clusters sizes distribution was not found")

    hist_all_big = os.path.join(params.output_dir, "big_cluster_sizes.png")
    drawing_utils.DrawAnyClusterSizesHist([[sz for sz in sizes if sz > params.size_cutoff] for sizes in 
                                           all_cluster_sizes], 
                                          rep_names, hist_all_big)
    if os.path.exists(hist_all_big):
        log.info("* Histogram of all big clusters sizes distribution was written to " + hist_all_big)
    else:
        log.info("Histogram of all big clusters sizes distribution was not found")

    for rep_ids, metrics in params.inexact_metrics.items():
        isolated_cluster_sizes = [sz for i, sz in enumerate(metrics.isolated_cluster_sizes)
                                  if i in metrics.rep_ids]
        curr_rep_names = [name for i, name in enumerate(rep_names) if i in metrics.rep_ids]
        hist_isol_all = os.path.join(params.output_dir, 
            "isolated_cluster_sizes_" + "".join(str(i + 1) for i in rep_ids) + ".png")
        drawing_utils.DrawAnyClusterSizesHist(isolated_cluster_sizes, curr_rep_names, hist_isol_all)
        if os.path.exists(hist_isol_all):
            log.info("* Histogram of isolated cluster sizes distribution for " + 
                ", ".join(curr_rep_names) + " repertoires was written to " + hist_isol_all)
        else:
            log.info("Histogram of isolated cluster sizes distribution for " + 
                ", ".join(curr_rep_names) + " repertoires was not found")

        hist_isol_nt = os.path.join(params.output_dir, 
            "isolated_nt_cluster_sizes_" + "".join(str(i + 1) for i in rep_ids) + ".png")
        drawing_utils.DrawAnyClusterSizesHist([[sz for sz in sizes if sz > 1] for sizes in 
                                               isolated_cluster_sizes],
                                              curr_rep_names, hist_isol_nt)
        if os.path.exists(hist_isol_nt):
            log.info("* Histogram of non-trivial isolated cluster sizes distribution for " + 
                ", ".join(curr_rep_names) + " repertoires was written to " + hist_isol_nt)
        else:
            log.info("Histogram of non-trivial isolated cluster sizes distribution for " + 
                ", ".join(curr_rep_names) + " repertoires was not found")

        hist_isol_big = os.path.join(params.output_dir, 
            "isolated_big_cluster_sizes_" + "".join(str(i + 1) for i in rep_ids) + ".png")
        drawing_utils.DrawAnyClusterSizesHist([[sz for sz in sizes if sz > params.size_cutoff] for sizes in 
                                               isolated_cluster_sizes],
                                              curr_rep_names, hist_isol_big)
        if os.path.exists(hist_isol_big):
            log.info("* Histogram of big isolated cluster sizes distribution for " + 
                ", ".join(curr_rep_names) + " repertoires was written to " + hist_isol_big)
        else:
            log.info("Histogram of big isolated cluster sizes distribution for " + 
                ", ".join(curr_rep_names) + " repertoires was not found")

def DrawClusterLengthHist(params, log):
    log.info("\n==== Drawing comparative histograms for all constructed cluster lengths distribution")
    all_cluster_lengths = []
    rep_names = []
    hist_all = os.path.join(params.output_dir, "constructed_cluster_lengths.png")
    for rep in params.evaluator.repertoires:
        all_cluster_lengths.append([len(cluster.seq) for cluster in rep.clusters.values()])
        rep_names.append(rep.name)

    if params.evaluator.original_repertoire:
        all_cluster_lengths.append([len(cluster.seq) for cluster in params.evaluator.original_repertoire.clusters.values()])
        rep_names.append("Original repertoire")

    drawing_utils.DrawClusterLengthsHist(all_cluster_lengths, rep_names, hist_all)
    if os.path.exists(hist_all):
        log.info("* Histogram of all clusters lengths distribution was written to " + hist_all)
    else:
        log.info("ERROR: files with histograms of all clusters lengths distribution was not found")
        sys.exit(1)

    for rep_ids, metrics in params.inexact_metrics.items():
        isolated_cluster_lengths = [sz for i, sz in enumerate(metrics.isolated_cluster_lengths)
                                  if i in metrics.rep_ids]
        curr_rep_names = [name for i, name in enumerate(rep_names) if i in metrics.rep_ids]
        hist_isol = os.path.join(params.output_dir, 
            "isolated_cluster_lengths_" + "".join(str(i+1) for i in rep_ids) + ".png")
        drawing_utils.DrawClusterLengthsHist(isolated_cluster_lengths, curr_rep_names, hist_isol)
        if os.path.exists(hist_isol):
            log.info("* Histogram of isolated cluster lengths distribution for " + 
                ", ".join(curr_rep_names) + " repertoires was written to " + hist_isol)
        else:
            log.info("* Histogram of isolated cluster lengths distribution for " + 
                ", ".join(curr_rep_names) + " repertoires was not found")

def DrawClusterGroupsDistribution(params, log):
    log.info("\n==== Drawing histograms for cluster groups")
    for rep_ids, metrics in params.inexact_metrics.items():
        if not metrics.has_reads:
            continue
        curr_rep_names = [rep.name for i, rep in 
                          enumerate(params.evaluator.repertoires) if i in rep_ids]
        chart_fname = os.path.join(params.output_dir, 
            "cluster_groups_" + "_".join(params.evaluator.repertoires[i].name for i in rep_ids) + ".png")
        data = [metrics.trusted_groups, metrics.untrusted_groups, metrics.big_untrusted_groups]
        labels = ["trusted_groups", "untrusted_groups", "big_untrusted_groups"]
        drawing_utils.DrawClusterGroupsBarChart(data, labels, chart_fname)
        if os.path.exists(chart_fname):
            log.info("* Histogram of all groups distribution for " + 
                ", ".join(curr_rep_names) + " repertoires was written to " + chart_fname)
        else:
            log.info("ERROR: files with histograms of all groups distribution for " +
                ", ".join(curr_rep_names) + " was not found")
            sys.exit(1)
        nt_chart_fname = os.path.join(params.output_dir, 
            "cluster_groups_nt_" + "_".join(params.evaluator.repertoires[i].name for i in rep_ids) + ".png")
        data = [metrics.trusted_groups_nt, metrics.untrusted_groups_nt, metrics.big_untrusted_groups]
        labels = ["trusted_groups_nt", "untrusted_groups_nt", "big_untrusted_groups"]
        drawing_utils.DrawClusterGroupsBarChart(data, labels, nt_chart_fname)
        if os.path.exists(nt_chart_fname):
            log.info("* Histogram of all non-trivial groups distribution for " + 
                ", ".join(curr_rep_names) + " repertoires was written to " + nt_chart_fname)
        else:
            log.info("ERROR: files with histograms of all non-trivial groups distribution for " +
                ", ".join(curr_rep_names) + " was not found")
            sys.exit(1)

def DrawHistograms(params, log):
    if not params.draw_hist:
        log.info("\nDrawing histograms was skippped")
        return 

    log.info("\n==== Drawing histograms")

    if params.original_clusters != ["", ""]:
        DrawIdentityMetrics(params, log)

    DrawClusterSizesHist(params, log)
    DrawClusterLengthHist(params, log)

    if params.original_clusters == ["", ""] and \
            (len(params.constructed_clusters) >= 2 or params.assembled_barcodes[0]):
        DrawClusterGroupsDistribution(params, log)

# -------------------------------------------------------------------

def CheckParams(params, log):
    for constructed_clusters_fa, constructed_clusters_rcm in \
            params.constructed_clusters:
        if params.original_clusters == ["", ""] and len(params.constructed_clusters) == 1 and \
                not constructed_clusters_rcm:
            log.info("ERROR: RCM file for CLUSTERS " + constructed_clusters_fa + " must be given")
            sys.exit(1)
        
        elif params.original_clusters == ["", ""] and len(params.constructed_clusters) != 1 and \
                not constructed_clusters_fa:
            log.info("ERROR: CLUSTERS file for RCM " + constructed_clusters_rcm + " must be given")
            sys.exit(1)
        elif params.original_clusters != ["", ""] and (not constructed_clusters_fa or not constructed_clusters_rcm):
            log.info("ERROR: CLUSTERS and RCM files for constructed clusters must be given")
            sys.exit(1)
        if constructed_clusters_fa and not os.path.exists(constructed_clusters_fa):
            log.info("ERROR: CLUSTERS with constructed repertoire " + constructed_clusters_fa + " was not found")
            sys.exit(1)
        if constructed_clusters_rcm and not os.path.exists(constructed_clusters_rcm):
            log.info("ERROR: RCM with constructed repertoire " + constructed_clusters_rcm + " was not found")
            sys.exit(1)

    if params.original_clusters != ["", ""]:
        if params.original_clusters[0] == "":
            log.info("ERROR: only RCM file for ideal repertoire have given, please provide CLUSTERS file for ideal repertoire")
            sys.exit(1)
        if params.original_clusters[1] == "":
            log.info("ERROR: only CLUSTERS file for ideal repertoire have given, please provide RCM file for ideal repertoire")
            sys.exit(1)
        if not os.path.exists(params.original_clusters[0]) or \
                not os.path.exists(params.original_clusters[1]):
            if not os.path.exists(params.original_clusters[0]):
                log.info("ERROR: CLUSTERS with ideal repertoire " + params.original_clusters[0] + " was not found")
            if not os.path.exists(params.original_clusters[1]):
                log.info("ERROR: RCM with ideal repertoire " + params.original_clusters[1] + " was not found")
            sys.exit(1)
    if params.tau < 0:
        log.info("ERROR: value for --tau must be >= 0")
        sys.exit(1)
    if params.size_cutoff < 0:
        log.info("ERROE: value for --size-cutoff must be >= 0")
        sys.exit(1)
    if params.min_cluster_size < 0:
        log.info("ERROE: value for --adv-min-size must be >= 0")
        sys.exit(1)


def PrintMainOutput(params, log):
    log.info("\nMain output files:")
    log.info("Computed metrics in TXT were written to " + params.file_metrics['general'])
    log.info("Computed metrics in CSV were written to " + params.file_metrics_csv['general'])

# ----------------------------------------------------------------------------------------

def RunBlast(clusters_fa, output_dir, log):   
    blast_output = os.path.join(output_dir, "blast.output") 
    log.info("\n== BLAST starts")
    command_line = init.RunBlast() + " -db " + init.BlastDB() + " -query " + clusters_fa + " -num_alignments 3 -outfmt 7 > " + blast_output
    log.info("BLAST command line: " + command_line)
    err_code = os.system(command_line)
    if err_code != 0:
        log.info("ERROR: BLAST finished abnormally")
        sys.exit(1)
    log.info("* BLAST output was written to " + blast_output)
    return blast_output

def RunIgBlast(clusters_fa, output_dir, log):
    igblast_output = os.path.join(output_dir, "igblast.output")
    log.info("\n== IgBlast starts")
    igblast_command_line = init.RunIgblast() + " -germline_db_V "+ init.ig_blast_directory +"database/human_gl_V -germline_db_J "+ init.ig_blast_directory +"database/human_gl_J -germline_db_D "+ init.ig_blast_directory +\
"database/human_gl_D -query "+ clusters_fa + " -show_translation -auxiliary_data auxilary_file -num_alignments 10  -outfmt \"7 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore slen\" -domain_system imgt > " + igblast_output
    os.environ['IGDATA'] = init.ig_blast_directory
    print("IgBlast command line: " + igblast_command_line)
    err_code = os.system(igblast_command_line)
    if err_code != 0:
        log.info("ERROR: IgBlast finished abnormally")
        sys.exit(1)
    log.info("* IgBlast output was written to " + igblast_output)
    return igblast_output

def PrepareDAOutputDir(clusters_fa, num_repertoire, params):
    output_dir = os.path.join(params.output_dir, "advanced_analysis_" + str(num_repertoire))
    os.makedirs(output_dir)
    return output_dir

def RunSomaticMutationSearch(params, clusters_fa, blast_output, igblast_output, output_dir, log):
    search_dir = os.path.join(output_dir, "somatic_search_results")
    log.info("\n== Somatic mutation search starts")
    log.info("Results will be written to " + search_dir)
    somatic_mutation_search.SomaticMutationSearch(clusters_fa, blast_output, igblast_output, search_dir, params.min_cluster_size, log)
 
def RunDeepAnalysisForSingleRepertoire(clusters_fa, num_repertoire, params, log):
    log.info("\n== Analysis of somatic mutations for repertoire " + clusters_fa)
    output_dir = PrepareDAOutputDir(clusters_fa, num_repertoire, params)  
    log.info("Results will be written to " + output_dir)  
    blast_output = RunBlast(clusters_fa, output_dir, log)
    igblast_output = RunIgBlast(clusters_fa, output_dir, log)
    RunSomaticMutationSearch(params, clusters_fa, blast_output, igblast_output, output_dir, log)

def RunDeepAnalysis(params, log):
    log.info("\n==== Analysis of somatic mutations starts")
    num_repertoire = 1
    for repertoire in params.constructed_clusters:
        clusters_fa = repertoire[0]
        if clusters_fa != "":
            print("Processing repertoire # " + str(num_repertoire) + ": " + clusters_fa)
            RunDeepAnalysisForSingleRepertoire(clusters_fa, num_repertoire, params, log)
            num_repertoire += 1

# ----------------------------------------------------------------------------------------

def RunIgQUAST(params, log):
    log.info("\n======== IgQUAST starts")

    # run repertoire evaluator
    log.info("\n==== Computation of metrics")
    ReadRepertoires(params, log)
    CreateRepertoireSubdirs(params, log)
    if params.original_clusters != ["", ""]:
        params.evaluator.original_repertoire = reading_utils.read_repertoires(
            [params.original_clusters])[0][0]
        RunExactEvaluator(params, log)
    elif len(params.constructed_clusters) == 1 and not params.assembled_barcodes[0]:
        RunGeneralEvaluator(params, log)
    else:
        RunInexactEvaluator(params, log)

    # print histograms
    DrawHistograms(params, log)

    if params.deep_analysis:
        RunDeepAnalysis(params, log)

    '''
     output clusters sizes
    CreateSizesOfClusters(params, log)
    '''

    log.info("\n======== IgQUAST ends")
    PrintMainOutput(params, log)
    log.info("\nThank you for using IgQUAST!")

def main():

    sys.stdout=init.Unbuffered(sys.stdout)

    # prepare log
    log = logging.getLogger('ig_quast')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    params = Params()
    ParseCommandLine(params, log)
    
    # preparation of output directory
    PrepareOutputDir(".", params)

    log_filename = os.path.join(params.output_dir, "ig_quast.log")
    if os.path.exists(log_filename):
        os.remove(log_filename)
    log_handler = logging.FileHandler(log_filename, mode='a')
    log.addHandler(log_handler)
    params.log = log_filename
    log.info("\nLog will be written to " + log_filename + "\n")

    init.PrintCommandLine(sys.argv, log)

    CheckParams(params, log)

    try:
        RunIgQUAST(params, log)
    except (KeyboardInterrupt):
        log.info("\nIgQUAST was interrupted!")
    except BaseException:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught. Please contact us and send .log file")

    log.info("\nLog was written to " + log_filename)

if __name__ == '__main__':
    main()
