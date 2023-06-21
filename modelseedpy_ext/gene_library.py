class GeneLibrary:

    def __init__(self, re):
        self.re = re
        self.genes = {}
        self.gene_id_to_genome_id = {}

    def index_genome(self, genome_id, genome):
        for f in genome.features:
            h = self.re.protein_store.get_sequence_hash(f.seq)
            if h not in self.genes:
                self.genes[h] = {'seq': f.seq, 'proteins': {}, 'annotation': {}}
            if genome_id not in self.genes[h]['proteins']:
                self.genes[h]['proteins'][genome_id] = set()
            self.genes[h]['proteins'][genome_id].add(f.id)

            if f.id not in self.gene_id_to_genome_id:
                self.gene_id_to_genome_id[f.id] = set()
            self.gene_id_to_genome_id[f.id].add(genome_id)

    def get_duplicate_gene_names(self):
        multi = set()
        for g_id in self.gene_id_to_genome_id:
            if len(self.gene_id_to_genome_id[g_id]) > 1:
                multi.add(g_id)

        return multi

    def get_freq(self):
        freq = {}
        for h in self.genes:
            proteins = self.genes[h]['proteins']
            count = 0
            for genome_id in proteins:
                count += len(proteins[genome_id])
            if count not in freq:
                freq[count] = 0
            freq[count] += 1

        return freq
