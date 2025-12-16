from modelseedpy_ext.re.core.genome import REAssembly


class IdentityVault:

    def __init__(self, mongo_database):
        self.mongo_database = mongo_database

    def get_assembly_hash(self, filename):
        doc = self.mongo_database['file_to_h'].find_one({'_id': filename})
        if doc is None:
            assembly = REAssembly.from_fasta(filename)
            self.mongo_database['file_to_h'].insert_one({'_id': filename, 'h': assembly.hash_value})
            doc = self.mongo_database['file_to_h'].find_one({'_id': filename})

        if doc:
            return doc['h']
        else:
            raise ValueError('unable to hash')


class NonRedundantGenomeVault:

    def __init__(self, mongo_database):
        self.mongo_database = mongo_database
        """
        from modelseedpy_ext.ani.nr_expansion import NrExp
        from modelseedpy_ext.ani.skani import ProgramSkani
        from modelseedpy_ext.ani import GenomeReference, ANIMethodSkani

        ani_method = ANIMethodSkani(ProgramSkani('/opt/skani/skani'), 80)
        nr_expansion = NrExp()
        nr_expansion.ani_radius = 95
        nr_expansion.ani_method = ani_method
        nr_expansion.MOCK_PRE_COMP_DATA = datasets
        nr_expansion.rep_clusters = {i: {} for i in gtdb_r220_rep}
        nr_expansion.rep_exp_libraries = {'gtdb_r220_rep': None}
        print(len(nr_expansion.rep_clusters))
        """
        pass

    def add_representative_genome(self, genome: REAssembly, type_designation, type_designation_source, method):
        if type(genome) is not REAssembly:
            raise TypeError('genome type must be REAssembly')

        contig_set_h = genome.hash_value
        self.mongo_database['representative_genomes'].insert_one({
            '_id': contig_set_h,
            'type_designation': type_designation,
            'type_designation_source': type_designation_source,
            'method': method,
        })
        self.mongo_database[f'genome_clusters_95_{contig_set_h}'].insert_one({
            '_id': contig_set_h,
            'ani': 100,
            'method': 'identity',
        })

    def get_genome_cluster(self, genome: REAssembly, ani: int):
        if type(genome) is not REAssembly:
            raise TypeError('genome type must be REAssembly')

        contig_set_h = genome.hash_value

        for doc in self.mongo_database[f'genome_clusters_{ani}_{contig_set_h}'].find():
            print(doc)
        pass


class DataVault:

    def __init__(self):
        pass

    def get_ncbi_record(self):
        pass

    def get_kbase_record(self):
        pass

    def get_assembly(self):
        pass


class SeedVault:

    def __init__(self, mongo_database):
        self.identity = IdentityVault(mongo_database)
        self.nr_genomes = NonRedundantGenomeVault(mongo_database)
