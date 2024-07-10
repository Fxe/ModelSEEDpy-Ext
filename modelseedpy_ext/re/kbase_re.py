from cobrakbase.kbaseapi import KBaseAPI


class KBaseRE(KBaseAPI):

    def __init__(self, re, token=None, dev=False, config=None):
        super().__init__(token, dev, config)

    def get_from_ws(self, id_or_ref, workspace=None):
        super().get_from_ws(id_or_ref, workspace)


class NodeTaxonomy:

    def __init__(self, re, doc):
        self.doc = doc
        self.re = re

    @property
    def id(self):
        return self.doc['_key']

    @property
    def label(self):
        return self.doc['_id'].split('/', 1)[0]

    @property
    def rank(self):
        return self.doc['rank']

    @property
    def genomes(self):
        return

    @property
    def parent(self):
        doc = self.re.db['ncbi_taxonomy'][self.doc['parent_tax_id']]
        return NodeTaxonomy(self.re, doc)

    @property
    def child_taxa(self):
        res = []
        taxa_edges = self.re.db.collections['ncbi_taxonomy_has_parent_ncbi_taxonomy']
        for o in taxa_edges.fetchByExample({'_to': self.doc['_id']}, batchSize=10, rawResults=True):
            n = self.re.db['ncbi_taxonomy'][o['_from'].split('/')[1]]
            res.append(NodeTaxonomy(self.re, n))

        return res

    @property
    def lineage(self):
        cursor = self

        res = [cursor.scientific_name]

        max_it = 30
        it = 0
        while str(cursor.id) != str(cursor.doc['parent_tax_id']) and it < max_it:
            cursor = cursor.parent
            res.append(cursor.scientific_name)
            it += 1

        res.reverse()
        return '; '.join(res)

    @property
    def scientific_name(self):
        return self.doc['scientific name'][0][0]