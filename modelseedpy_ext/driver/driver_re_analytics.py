import logging
from tqdm import tqdm
from modelseedpy_ext.re.hash_seq import HashSeqList, HashSeq
from modelseedpy import MSGenome
from modelseedpy_ext.utils import sha_hex

logger = logging.getLogger(__name__)


def _process_contigs(contigs):
    hash_list = HashSeqList()
    contig_h_d = []
    for contig in contigs.features:
        seq = HashSeq(contig.seq)
        hash_list.append(seq)
        seq_h = seq.hash_value
        contig_h_d.append([seq_h, contig.id, contig.description])
    return {
        'genome_h': hash_list.hash_value,
        'contig_h': contig_h_d
    }


def process_files(files):
    h_to_contig = {}
    h_to_genome = {}
    contig_descs = {}
    genome_to_contig_h = {}
    genome_to_contig_id = {}
    genome_error = set()
    file_to_hash = {}
    for g in tqdm(files):
        if g not in genome_to_contig_h and g not in genome_error:
            if g.endswith('.fna.gz') or True:
                try:
                    path_fna = g
                    genome_contigs = MSGenome.from_fasta2(path_fna, split=' ')
                    hvals = _process_contigs(genome_contigs)
                    file_to_hash[path_fna] = hvals['genome_h']
                    if hvals['genome_h'] not in h_to_genome:
                        h_to_genome[hvals['genome_h']] = set()
                        genome_to_contig_h[hvals['genome_h']] = []
                        genome_to_contig_id[hvals['genome_h']] = []
                    h_to_genome[hvals['genome_h']].add(path_fna)
                    for contig_h, contig_id, des in hvals['contig_h']:
                        genome_to_contig_h[hvals['genome_h']].append(contig_h)
                        genome_to_contig_id[hvals['genome_h']].append(contig_id)
                        if contig_h not in h_to_contig:
                            h_to_contig[contig_h] = set()
                        h_to_contig[contig_h].add(contig_id)
                        if contig_id not in contig_descs:
                            contig_descs[contig_id] = set()
                        contig_descs[contig_id].add(des)
                except Exception as exxx:
                    print(exxx)
                    genome_error.add(g)
                #if len(h_to_genome) > 10:
                #    return file_to_hash, h_to_genome, contig_descs, h_to_contig, genome_to_contig_h, genome_to_contig_id, genome_error
            else:
                genome_error.add(g)
        else:
            genome_error.add(g)
    return file_to_hash, h_to_genome, contig_descs, h_to_contig, genome_to_contig_h, genome_to_contig_id, genome_error


def _process_old_style(gtdb_path, w_path='/scratch/fliu/data/gtdb/GCA/'):
    h_to_contig = {}
    h_to_genome = {}
    genome_to_contig_h = {}
    counter = 0
    genome_error = set()

    import os

    for x in os.listdir(gtdb_path):
        x_path = gtdb_path + '/' + x
        if os.path.isdir(x_path):
            for y in os.listdir(x_path):
                y_path = x_path + '/' + y
                if os.path.isdir(y_path):
                    for z in os.listdir(y_path):
                        z_path = y_path + '/' + z
                        if os.path.isdir(z_path):
                            for g in os.listdir(z_path):
                                if g.startswith('GCA') and g not in genome_to_contig_h and g not in genome_error:
                                    for filename in os.listdir(z_path + '/' + g):
                                        if filename.endswith('.fna.gz'):
                                            try:
                                                genome_to_contig_h[g] = set()
                                                path_fna = f'{z_path}/{g}/{filename}'
                                                genome = MSGenome.from_fasta2(path_fna, split=' ')
                                                w_group = counter // 10000

                                                h_seq = HashSeqList()
                                                for f in genome.features:
                                                    hash_seq = HashSeq(f.seq)
                                                    _h = hash_seq.hash_value
                                                    genome_to_contig_h[g].add(_h)
                                                    if _h not in h_to_contig:
                                                        h_to_contig[_h] = set()
                                                    h_to_contig[_h].add(f.id)
                                                    h_seq.append(hash_seq)
                                                hash_contig_set = h_seq.hash_value

                                                if hash_contig_set not in h_to_genome:
                                                    h_to_genome[hash_contig_set] = set()
                                                    d_group = f'{w_path}/genome_set_{w_group}'
                                                    if not os.path.exists(d_group):
                                                        os.mkdir(d_group)
                                                    genome.to_fasta(f'{d_group}/{hash_contig_set}.fna')
                                                    counter += 1
                                                h_to_genome[hash_contig_set].add(g)
                                            except Exception as exp:
                                                genome_error.add(g)

    return h_to_genome, h_to_contig, genome_to_contig_h, counter, genome_error


class KeggGenomeTable:

    def __init__(self, re):
        self.re = re
        col_kegg_genomes = re.db['kegg_genome']
        self.kegg_genomes = {}
        for o in col_kegg_genomes.fetchAll(rawResults=True):
            self.kegg_genomes[o['_key']] = o

        logger.info(f'kegg_genome [{len(self.kegg_genomes)}]')

        self.kegg_genome_ids = {x['_id'] for x in self.kegg_genomes.values()}

        self.kegg_genomes_to_ncbi_taxonomy = self.re.get_link('kegg_genome_has_ncbi_taxonomy', self.kegg_genome_ids)
        self.ncbi_taxonomy = self.get_docs('ncbi_taxonomy', set(self.kegg_genomes_to_ncbi_taxonomy.values()))

        logger.info(f'KEGG ncbi_taxonomy [{len(self.ncbi_taxonomy)}]')

        self.kegg_genome_has_ncbi_assembly_gca_acc = self.re.get_link('kegg_genome_has_ncbi_assembly_gca_acc',
                                                                      self.kegg_genome_ids)
        self.kegg_genome_has_ncbi_assembly_gcf_acc = self.re.get_link('kegg_genome_has_ncbi_assembly_gcf_acc',
                                                                      self.kegg_genome_ids)

        logger.info(f'KEGG ncbi_assembly_gca_acc [{len(self.kegg_genome_has_ncbi_assembly_gca_acc)}]')
        logger.info(f'KEGG ncbi_assembly_gcf_acc [{len(self.kegg_genome_has_ncbi_assembly_gcf_acc)}]')

        self.ncbi_assembly_gca_to_gcf = self.re.get_link('ncbi_assembly_gca_to_gcf',
                                                         set(self.kegg_genome_has_ncbi_assembly_gca_acc.values()),
                                                         rev=False)
        self.ncbi_assembly_has_gca = self.re.get_link('ncbi_assembly_has_ncbi_assembly_gca_acc',
                                                      set(self.kegg_genome_has_ncbi_assembly_gca_acc.values()),
                                                      rev=True)

        logger.info(f'GCA to GCF [{len(self.ncbi_assembly_gca_to_gcf)}]')

        self.all_gcf = set(self.kegg_genome_has_ncbi_assembly_gcf_acc.values())
        self.all_gcf.update(self.ncbi_assembly_gca_to_gcf.keys())
        self.ncbi_assembly_has_gcf = self.re.get_link('ncbi_assembly_has_ncbi_assembly_gcf_acc',
                                                      self.all_gcf,
                                                      rev=True)

        logger.info(f'total GCF [{len(self.all_gcf)}]')
        logger.info(f'ncbi_assembly from GCA [{len(self.ncbi_assembly_has_gca)}]')
        logger.info(f'ncbi_assembly from GCF [{len(self.ncbi_assembly_has_gcf)}]')

        self.all_assemblies = set(self.ncbi_assembly_has_gca.values())
        self.all_assemblies.update(self.ncbi_assembly_has_gcf)
        self.ncbi_assembly = self.get_docs('ncbi_assembly', self.all_assemblies)

        logger.info(f'ncbi_assembly [{len(self.ncbi_assembly)}]')

        self.gca_to_gcf = {}
        for a, b in self.ncbi_assembly_gca_to_gcf.items():
            if b not in self.gca_to_gcf:
                self.gca_to_gcf[b] = a
            else:
                print(a)

        self.gca_to_assembly = {}
        for gca, a in self.ncbi_assembly_has_gca.items():
            if gca not in self.gca_to_assembly:
                self.gca_to_assembly[gca] = a
            else:
                print(gca, a)

    def get_docs(self, collection_id, filter_list):
        aql = f"""
        FOR d IN {collection_id}
            FILTER d._id IN @list
            RETURN d
        """
        cursor = self.re.db.AQLQuery(aql, rawResults=True, bindVars={'list': list(filter_list)})
        res = {}
        for doc in cursor:
            if doc['_id'] not in res:
                res[doc['_id']] = doc
            else:
                raise ValueError('!!!')
        return res

    @staticmethod
    def diagram():
        s = """
        kegg_genome -[kegg_genome_has_ncbi_taxonomy]-> ncbi_taxonomy
         |
         |
         |-----[kegg_genome_has_ncbi_assembly_gca_acc]--> ncbi_assembly_gca_acc
         |                                                         |
         |                                              [ncbi_assembly_gca_to_gcf]
         |                                                         |
         |                                                         v
         |-----[kegg_genome_has_ncbi_assembly_gcf_acc]--> ncbi_assembly_gcf_acc
                                                                   |  
                                                                   | 
                                                [ncbi_assembly_has_ncbi_assembly_gcf_acc]
                                                                   |
                                                                   v
                                                            ncbi_assembly
        """
        return s
