class RowCluster:

    def __init__(self, cluster_id, size, is_core, members, function, ec):
        self.cluster_id = cluster_id
        self.size = size
        self.is_core = is_core
        self.members = members
        self.function = function
        self.ec = ec


class RowFeature:

    def __init__(self, genome_id, contig_id, feature_id, length, start, end, strand, annotation):
        self.genome_id = genome_id
        self.feature_id = feature_id
        self.contig_id = contig_id
        self.length = length
        self.start = start
        self.end = end
        self.strand = strand
        self.annotation = annotation