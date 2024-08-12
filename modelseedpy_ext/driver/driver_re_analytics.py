import logging

logger = logging.getLogger(__name__)


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
