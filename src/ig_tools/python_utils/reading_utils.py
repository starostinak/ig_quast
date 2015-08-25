from Bio import SeqIO
from repertoire_structs import Repertoire, Cluster
from collections import defaultdict
import sys
import os.path

def read_repertoire(fasta, rcm, read_ids):
    clusters = read_fasta(fasta)
    cluster_reads, read_clusters = None, None
    if rcm:
        cluster_reads, read_clusters = read_rcm(rcm, read_ids)
    return clusters, cluster_reads, read_clusters

def read_repertoires(repertoires_list, barcode_files):
    repertoires = []
    read_ids = {}
    if barcode_files[0]:
        clusters, cluster_reads, read_clusters = read_repertoire(barcode_files[0], barcode_files[1], read_ids)
        repertoires.append(Repertoire(barcode_files[0], barcode_files[1], clusters, 
                                       cluster_reads, read_clusters, 'assembled_barcodes'))
    for i, (fasta, rcm) in enumerate(repertoires_list):
        clusters, cluster_reads, read_clusters = read_repertoire(fasta, rcm, read_ids)
        name = os.path.basename(fasta).split('.', 1)[0] + \
            ('_' + str(i + 1) if len(repertoires_list) > 1 else '')
        repertoires.append(Repertoire(fasta, rcm, clusters, cluster_reads, read_clusters, name))
    read_names = [''] * len(read_ids.keys())
    for name, id in read_ids.items():
        read_names[id] = name
    return repertoires, read_names

def read_fasta(fasta):
    clusters = {}
    fasta_records = SeqIO.parse(open(fasta), 'fasta')
    for i, rec in enumerate(fasta_records):
        fields = rec.id.split('___')
        if len(fields) != 4:
            print "ERROR: wrong format of cluster id " + rec.id + ", must be in format cluster___id___size___sz"
            sys.exit(1)
        clusters[int(fields[1])] = Cluster(int(fields[3]), rec.seq)
    return clusters

def read_rcm(rcm, read_ids):
    fill_read_ids = True
    read_clusters = []
    cluster_reads = {}
    if read_ids:
        fill_read_ids = False
        read_clusters = [-1] * len(read_ids.keys())
    handle = open(rcm)
    for line in handle:
        fields = line.strip().split('\t')
        read_id_str = fields[0].split('___')[0]
        read_id_str = read_id_str.split(',')[0]
        if len(fields) != 2:
            print "ERROR: wrong line in rcm: " + line.strip() + ", must be in format read_id\\tcluster_no"
            sys.exit(1)
        cluster_no = int(fields[1])
        if fill_read_ids or read_id_str not in read_ids:
            read_ids[read_id_str] = len(read_clusters)
            read_id = len(read_clusters)
            read_clusters.append(-1)
        else:
            read_id = read_ids[read_id_str]
        if cluster_no not in cluster_reads:
            cluster_reads[cluster_no] = set()
        cluster_reads[cluster_no].add(read_id)
        read_clusters[read_id] = cluster_no
    return cluster_reads, read_clusters

def read_cluster_numbers_to_ids(filtered_clusters_fa):
    cluster_num_to_id = {}
    fasta_records = SeqIO.parse(open(filtered_clusters_fa), 'fasta')
    for i, rec in enumerate(fasta_records):
        fields = rec.id.split('___')
        if len(fields) != 4:
            print "ERROR: wrong format of cluster id " + rec.id + ", must be in format cluster___id___size___sz"
            sys.exit(1)
        cluster_num_to_id[i] = int(fields[1])
    return cluster_num_to_id

def read_neighbour_clusters(matches_filenames, repertoires, cluster_num_to_ids):
    cluster_pairs = [[None] * len(repertoires) for i in xrange(len(repertoires))]
    for rep_id1, filename in enumerate(matches_filenames):
        handle = open(filename)
        rep_id2 = 1 - rep_id1
        for cluster_no1, l in enumerate(handle):
            fields = l.strip().split(' ')
            if not fields:
                continue
            score = fields[0]
            cluster_nums = fields[1:]
            rep1, cluster1 = rep_id1, cluster_num_to_ids[rep_id1][cluster_no1]
            for cluster_no2 in cluster_nums:
                cluster_no2 = int(cluster_no2)
                rep2, cluster2 = rep_id2, cluster_num_to_ids[rep_id2][cluster_no2]
                if not cluster_pairs[rep1][rep2]:
                    cluster_pairs[rep1][rep2] = {}
                if cluster1 not in cluster_pairs[rep1][rep2]:
                    cluster_pairs[rep1][rep2][cluster1] = [set(), 0]
                cluster_pairs[rep1][rep2][cluster1][0].add(cluster2)
                cluster_pairs[rep1][rep2][cluster1][1] = -int(score)
        
    return cluster_pairs
