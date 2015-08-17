from Bio import SeqIO
from repertoire_structs import Repertoire, Cluster
from collections import defaultdict
import sys
import os.path

def read_repertoires(repertoires_list):
    repertoires = []
    read_ids = {}
    for i, (fasta, rcm) in enumerate(repertoires_list):
        clusters = read_fasta(fasta)
        cluster_reads, read_clusters = None, None
        if rcm:
            cluster_reads, read_clusters = read_rcm(rcm, read_ids)
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
    for rec in fasta_records:
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

def read_neighbour_clusters(filename, repertoires):
    cluster_pairs = [[None] * len(repertoires) for i in xrange(len(repertoires))]
    handle = open(filename)
    handle.readline()

    for l in handle:
        fields = l.strip().split(' ')
        rep1, cluster1 = fields[0].split('.')
        rep1, cluster1 = int(rep1) - 1, int(cluster1)
        rep2, cluster2 = fields[2].split('.')
        rep2, cluster2 = int(rep2) - 1, int(cluster2)
        score = int(fields[5][:-1])
        if not cluster_pairs[rep1][rep2]:
            cluster_pairs[rep1][rep2] = {}
        if cluster1 not in cluster_pairs[rep1][rep2]:
            cluster_pairs[rep1][rep2][cluster1] = [set(), 0]
        cluster_pairs[rep1][rep2][cluster1][0].add(cluster2)
        cluster_pairs[rep1][rep2][cluster1][1] = score
    
    return cluster_pairs
