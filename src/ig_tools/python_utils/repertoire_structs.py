class Repertoire:
    def __init__(self, clusters_filename, rcm_filename, clusters, cluster_reads = None, read_clusters = None, name = None):
        self.clusters_filename = clusters_filename
        self.rcm_filename = rcm_filename
        self.clusters = clusters
        self.cluster_reads = cluster_reads
        self.read_clusters = read_clusters
        self.name = name if name else clusters_filename

    def __str__(self):
        return ('Clusters: ' + str(self.clusters) + ', \n reads: ' + 
                ('None' if not self.cluster_reads else self.cluster_reads))

    def __repr__(self):
        return str(self)

    def get_cluster_seq_length(self, cluster_id):
        return len(self.clusters[cluster_id].seq)

class Cluster:
    def __init__(self, size, seq):
        self.size = size
        self.seq = seq

    def __str__(self):
        return 'size ' + str(self.size) + ', seq ' + str(self.seq)

    def __repr__(self):
        return str(self)
