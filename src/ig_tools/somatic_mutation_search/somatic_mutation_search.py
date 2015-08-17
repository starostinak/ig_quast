#!/usr/bin/env python

import os
from Bio import pairwise2
import sys
import shutil
import logging

import blast_parser

sys.path.append("../python_utils")

import drawing_utils
import igblast_utils

# ----------------------------------------------------------

class Cluster:
    size = 0
    seq = ""

    def __init__(self, size, seq):
        self.size = size
        self.seq = seq

class Repertoire:
    id_cluster = dict()

    def ExtractClustersFromFile(self, fname):
        if not os.path.exists(fname):
            print("File with clusters " + fname + " was not found")
            sys.exit(1)
        lines = open(fname, "r").readlines()
        num_clusters = len(lines) / 2
        for i in range(0, num_clusters):
            header = lines[i * 2].strip()
            splits = header.split('___')
            cluster = Cluster(int(splits[3]), lines[i * 2 + 1].strip())
            self.id_cluster[int(splits[1])] = cluster

    def GetSeqById(self, cluster_id):
        if cluster_id in self.id_cluster:
            return self.id_cluster[cluster_id].seq
        return ""

    def GetFullNameById(self, cluster_id):
        if not cluster_id in self.id_cluster:
            print "ERROR: repertoire does not contain cluster with id " + str(cluster_id)
            sys.exit(1)
        size = self.id_cluster[cluster_id].size
        return "cluster___" + str(cluster_id) + "___size___" + str(size)

    def NonTrivialClusters(self, size = 1):
        nt_clusters = list()
        for cluster_id in self.id_cluster:
            if self.id_cluster[cluster_id].size >= size:
                nt_clusters.append(cluster_id)
        return nt_clusters

# ----------------------------------------------------------

class MutatedCandidates:
    candidates = dict()

    def GetCandidatesFor(self, cluster_id):
        if cluster_id in self.candidates:
            return self.candidates[cluster_id]
        return set()

    def MainClusters(self):
        return self.candidates.keys()

    def InitByBlastOutput(self, blast_output):
        match_map = dict()
        for cluster_id in blast_output.match_map:
            matches = blast_output.match_map[cluster_id]
            for m in matches:
                if m not in match_map:
                    match_map[m] = list()
                match_map[m].append(cluster_id)

#        print(match_map)

        for m in match_map:
            cluster_ids = match_map[m]
            for id1 in cluster_ids:
                if id1 not in self.candidates:
                    self.candidates[id1] = set()
                for id2 in cluster_ids:
                    if id1 != id2:
                        self.candidates[id1].add(id2)

# ----------------------------------------------------------

class AlignmentSearcher:
    def Align(self, s1, s2):
        return pairwise2.align.globalms(s1, s2, 1, -.1, -10, -.5)

def GetCDRsForBlocks(igblast_blocks):
    cdr1_start = sys.maxint
    cdr1_end = -1
    cdr2_start = sys.maxint
    cdr2_end = -1
    cdr3_start = sys.maxint
    cdr3_end = -1

    for block in igblast_blocks:
        cdr1_start = min(block.alignment_summary.cdr1.from_ind, cdr1_start)
        cdr1_end = max(block.alignment_summary.cdr1.to_ind, cdr1_end)

        cdr2_start = min(block.alignment_summary.cdr2.from_ind, cdr2_start)
        cdr2_end = max(block.alignment_summary.cdr2.to_ind, cdr2_end)

        cdr3_start = min(block.alignment_summary.cdr3.from_ind, cdr3_start)
        cdr3_end = max(block.alignment_summary.cdr3.to_ind, cdr3_end)

    return [cdr1_start, cdr1_end, cdr2_start, cdr2_end, cdr3_start, cdr3_end]

#--------------------------------------------------
#def AlignStringsOverlap(align_str1, align_str2):
#   tail1 = 0
#   for i in range(0, align_str1):
#        if align_str1[i] != '-':
#            tail1 = i
#            break
#
#    tail2 = len(align_str2)
#    for i in range(0, align_str2):
#        if align_str2[len(align_str2) - i - 1] != '-':
#            tail2 = len(align_str2) - i - 1
#            break
#
#    if len(align_str2) - tail1 != tail2 + 1
#        print("DEBUG ME!")
#        print("Str1 - " + align_str1)
#        print("Str2 - " + align_str2)
#        print("Tail1 - " + str(tail1))
#        print("Tail2 - " + str(tail2))
#        sys.exit(1)
#    return tail1, tail2
#
#def OverlapIsGood(align_str1, align_str2, tail1, tail2):
#    overlap_size = align_str2 - tail1
#    start1 = 0
#    start2 = tail1
#    for i in range(0, overlap_size):
#        if align_str1[i + start1] != align_str2[i + start2]:
#            return False
#    return True
#
#def AlignmentContainsShifts(align):
#    tail11, tail12 = AlignStringsOverlap(align[0], align[1])
#    tail21, tail22 = AlignStringsOverlap(align[1], align[0])
#    return OverlapIsGood(align[0], align[1], tail11, tail12) or OverlapIsGood(align[1], align[0], tail21, tail22)
#
#def AlignmentIsGood(align):
#    return AlignmentContainsShifts(align)
#def GetBestAlignment(alignments):
#    for align in alignments:

def AlignWithoutGaps(align_str1, align_str2):
    return len([s for s in align_str1 if s == '-']) == 0 and len([s for s in align_str2 if s == '-']) == 0

def GetNumberOfMismathces(align_str):
    num_mismatches = 0
    for s in align_str:
        if s == '_':
            num_mismatches += 1
    return num_mismatches
    
# ----------------------------------------------------------    

class SomaticMulationsSearcher:
    reperoire = Repertoire()
    candidated = MutatedCandidates()
    igblast_output = igblast_utils.RawIgblastOutput()

    output_dir = ""

    def CleanOutputDir(self):
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)

    def __init__(self, repertoire, candidates, igblast_output, output_dir):
        self.repertoire = repertoire
        self.candidates = candidates
        self.igblast_output = igblast_output
        self.output_dir = os.path.abspath(output_dir)
        self.CleanOutputDir()

    def GetAlignmentSymbolsString(self, align1, align2):
        align = ""
        for i in range(len(align1)):
            if align1[i] == align2[i]:
                align += "."
            else:
                if align1[i] == '-' or align2[i] == '-':
                    align += "!"
                else:
                    align += "_"
        return align

    def GetAlignmentNumbersStr(self, align_str, start_num):
        align_num_str = list()
        for s in align_str:
            if s == '.':
                align_num_str.append(start_num)
            else:
                align_num_str.append(start_num + 4)
        return align_num_str, start_num + 5

    def AlignmentIsGood(self, align_str):
        return len([s for s in align_str if s == "!"]) == 0

    def ProcessAlignment(self, cluster1, cluster2, align_str):
        output_dirname = os.path.join(self.output_dir, "cluster_" + str(cluster1))
        if not os.path.exists(output_dirname):
            os.mkdir(output_dirname)
        output_dirname = os.path.join(output_dirname, "pairwise")
        if not os.path.exists(output_dirname):
            os.mkdir(output_dirname)
        output_fname = os.path.join(output_dirname,  "cluster_" + str(cluster1) + "_cluster_" + str(cluster2))
        num_list = self.GetAlignmentNumbersStr(align_str, 0)[0]

        # draw figure
        settings = drawing_utils.GetGraphicalSettings(xlabel = "Alignment position", ylabel = "", title = "", output_filename = output_fname + ".png", show_yaxis = False, label = [""])

        # cdrs reconstruction
        cluster1_name = self.repertoire.GetFullNameById(cluster1)
        cluster2_name = self.repertoire.GetFullNameById(cluster2)

        cdrs = GetCDRsForBlocks([self.igblast_output.GetBlockByName(cluster1_name), self.igblast_output.GetBlockByName(cluster2_name)])
        drawing_utils.DrawSomaticMutations([range(0, len(num_list))], [num_list], cdrs, settings)

    def OutputAlignmentStats(self, align_str):
        #mismatches_num_fhandler = open(os.path.join(self.output_dir, "mismatches_num.txt"), "a")
        #mismatches_num_fhandler.write(str(len([s for s in align_str if s == '_'])) + "\n")

        mismatches_pos_fhandler = open(os.path.join(self.output_dir, "polymorphism_positions.txt"), "a")
        for i in range(0, len(align_str)):
            if align_str[i] == '_':
                mismatches_pos_fhandler.write(str(float(i) / len(align_str)) + "\n")

    def DrawMismatchesHistogram(self, log):
        mismatches_pos_fname = os.path.join(self.output_dir, "polymorphism_positions.txt")
        mismatches_pos = drawing_utils.ReadFloatGraphicalData(mismatches_pos_fname)
        output_fname = os.path.join(self.output_dir, "polymorphism_positions.png")
        settings = drawing_utils.GetGraphicalSettings(xlabel = "Relative position of mutation", ylabel = "Mutation frequency", title = "", output_filename = output_fname)
        drawing_utils.DrawMutationHistogram(mismatches_pos.all_keys, settings)

        log.info("\n* Relative positions of polymorphisms were written to " + mismatches_pos_fname)
        log.info("* Histogram of relative positions of polymorphisms was written to " + output_fname)

    def ProcessPartClusterAlign(self, cluster_id, X, Y, align_list, output_fname):
        labels = ["Cluster " + str(align_list[len(align_list) - i - 1][1]) for i in range(0, len(align_list))]
        settings = drawing_utils.GetGraphicalSettings(xlabel = "Alignment positions", ylabel = "", title = "", output_filename = output_fname, label = labels, show_yaxis = False, draw_legend = True, xmin_shift = 10, xmax_shift = 200, ymin_shift = 2, ymax_shift = 2, legend_loc = "center right")

        igblast_blocks = [self.igblast_output.GetBlockByName(self.repertoire.GetFullNameById(align[1])) for align in align_list]
        igblast_blocks.append(self.igblast_output.GetBlockByName(self.repertoire.GetFullNameById(cluster_id)))
        cdrs = GetCDRsForBlocks(igblast_blocks)
        drawing_utils.DrawSomaticMutations(X, Y, cdrs, settings)

        for align in align_list:
            self.OutputAlignmentStats(align[2])

    def ProcessClusterAlign(self, cluster_id, align_list):
        X = list()
        Y = list()
        seed = 0
        for align in align_list:
            y, seed = self.GetAlignmentNumbersStr(align[2], seed)
            x = range(1, len(y) + 1)
            X.insert(0, x)
            Y.insert(0, y)

        max_align = 10
        if len(align_list) < max_align:    
            output_fname = os.path.join(self.output_dir, "cluster_" + str(cluster_id) + "/cluster_" + str(cluster_id) + ".png")
            self.ProcessPartClusterAlign(cluster_id, X, Y, align_list, output_fname)
        else:
            num_plots = len(align_list) / max_align 
            if len(align_list) % max_align != 0:
                num_plots += 1
            for i in range(0, num_plots):
                start_ind = i * max_align
                end_ind = min((i + 1) * max_align, len(align_list))
                X_part = [X[j] for j in range(start_ind, end_ind)] 
                Y_part = [Y[j] for j in range(start_ind, end_ind)] 
                align_part = [align_list[j] for j in range(start_ind, end_ind)]
                output_fname = os.path.join(self.output_dir, "cluster_" + str(cluster_id) + "/cluster_" + str(cluster_id) + "_part_" + str(i + 1) + ".png")
                self.ProcessPartClusterAlign(cluster_id, X_part, Y_part, align_part, output_fname)
        return os.path.join(self.output_dir, "cluster_" + str(cluster_id))

    def GetBestAlignment(self, alignments):
        for i in range(0, len(alignments)):
            align_str1 = alignments[i][0]
            align_str2 = alignments[i][1]
            score = alignments[i][2]

            if score > 330:
                return i
        return -1

    def Search(self, min_cluster_size, log):
        aligner = AlignmentSearcher()
        main_clusters = self.repertoire.NonTrivialClusters(min_cluster_size) #self.candidates.MainClusters()
        log.info("Repertoire contains " + str(len(main_clusters)) + " clusters of size >=" + str(min_cluster_size))        

        for cluster1 in main_clusters:
            seq1 = self.repertoire.GetSeqById(cluster1)
            candidates = self.candidates.GetCandidatesFor(cluster1)

            aligns = list()
            for cluster2 in candidates:
                seq2 = self.repertoire.GetSeqById(cluster2)

                alignment = aligner.Align(seq1, seq2)
                if len(alignment) == 0:
                    continue
            
                best_align_ind = self.GetBestAlignment(alignment)
                if best_align_ind == -1:
                    continue

                align_seq1 = alignment[best_align_ind][0]
                align_seq2 = alignment[best_align_ind][1]

                if len(align_seq1) != len(align_seq2):
                    log.info("ERROR: alingment for clusters" + str(cluster1) + " & " + str(cluster2) + " is incorrect!")
                align = self.GetAlignmentSymbolsString(align_seq1, align_seq2)

                aligns.append([cluster1, cluster2, align])
                self.ProcessAlignment(cluster1, cluster2, align)

            if len(aligns) <= 1:
                continue
            output_cluster_dir = self.ProcessClusterAlign(cluster1, aligns)
            log.info("Results of " + str(len(aligns)) + " alignments for cluster " + str(cluster1) + " were written to " + output_cluster_dir)
        self.DrawMismatchesHistogram(log)

# ------------------------------------------------------------------------

def SomaticMutationSearch(clusters_fa, blast_fname, igblast_fname, output_dir, min_cluster_size, log):
    blast_output = blast_parser.SimpleBlastOutput(blast_fname)

    repertoire = Repertoire()
    repertoire.ExtractClustersFromFile(clusters_fa)

    mutated_candidates = MutatedCandidates()
    mutated_candidates.InitByBlastOutput(blast_output)

    igblast_output = igblast_utils.ParseIgBlastOutput(igblast_fname, log)

    searcher = SomaticMulationsSearcher(repertoire, mutated_candidates, igblast_output, output_dir)
    searcher.Search(min_cluster_size, log)

def main(): 
    if len(sys.argv) != 5:
        print "Invalid input parameters."
        print "./somatic_mutation_search.py repertoire.clusters.fa blast.output igblast.output output_dir"
        sys.exit(1)

    # prepare log
    log = logging.getLogger('somatic_mutation_searcher')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    clusters_fa = sys.argv[1] 
    blast_fname = sys.argv[2] 
    igblast_output = sys.argv[3] 

    blast_output = blast_parser.SimpleBlastOutput(blast_fname)

    repertoire = Repertoire()
    repertoire.ExtractClustersFromFile(clusters_fa)

    mutated_candidates = MutatedCandidates()
    mutated_candidates.InitByBlastOutput(blast_output)

    igblast_output = igblast_utils.ParseIgBlastOutput(igblast_output, log)

    output_dir = sys.argv[4] 
    searcher = SomaticMulationsSearcher(repertoire, mutated_candidates, igblast_output, output_dir)
    searcher.Search(log)

if __name__ == '__main__':
    main()

