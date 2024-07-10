from abc import ABC


def read_kegg_file(filename, encoding=None, debug=False):
    cursor = None
    with open(filename, 'r', encoding=encoding) as fh:
        data = {}
        line = fh.readline()
        while line:
            if debug:
                print(line[:-1])
            current = line[:12].strip()
            rest = line[12:]
            if current and current != cursor:
                cursor = current
                if cursor not in data:
                    data[cursor] = {'value': [], 'subfields': {}}
                data[cursor]['value'].append(rest)
            else:
                data[cursor]['value'].append(rest)

            if current == '///':
                return data
            line = fh.readline()


def transform_genomes(genomes):
    from modelseedpy_ext.re.etl.transform_graph import TransformGraph, Node
    graph = TransformGraph()
    copy = {'category', 'chromosome', 'code', 'comments', 'db_links', 'disease', 'lineage', 'name', 'statistics'}
    for g in genomes.values():
        data = {}
        for k in copy:
            if k in g:
                data[k] = g[k]
        n = graph.add_transform_node2(Node(g['entry'], 'kegg_genome', data=data))
        if 'taxonomy' in g and 'TAX:' in g['taxonomy']:
            node_ncbi = graph.add_transform_node2(
                Node(g['taxonomy'].split('TAX:')[1].strip(), 'ncbi_taxonomy', data={'_proxy': True}))
            graph.add_transform_edge2(n, node_ncbi, 'kegg_genome_has_ncbi_taxonomy')
        else:
            print(g['taxonomy'])
        if 'data_source' in g:
            if len(g['data_source']) == 2:
                ncbi_assembly, bio_proj = [s.strip() for s in g['data_source']]
                if bio_proj.startswith('BioProject: '):
                    bioproject_id = bio_proj.split(':')[1].strip()
                    node_ncbi = graph.add_transform_node2(Node(bioproject_id, 'ncbi_bioproject', data={'_proxy': True}))
                    graph.add_transform_edge2(n, node_ncbi, 'kegg_genome_has_ncbi_bioproject')
                else:
                    raise ValueError('invalid bioproject')
                if 'Assembly:' in ncbi_assembly:
                    # print(ncbi_assembly)
                    a, b = ncbi_assembly.split('Assembly:')
                    genome_id = b.strip().split(' ')[0].strip()
                    if ncbi_assembly.startswith('GenBank'):
                        node_ncbi = graph.add_transform_node2(Node(genome_id, 'ncbi_assembly_gca_acc'))
                        graph.add_transform_edge2(n, node_ncbi, 'kegg_genome_has_ncbi_assembly_gca_acc')
                    elif ncbi_assembly.startswith('RefSeq'):
                        node_ncbi = graph.add_transform_node2(Node(genome_id, 'ncbi_assembly_gcf_acc'))
                        graph.add_transform_edge2(n, node_ncbi, 'kegg_genome_has_ncbi_assembly_gcf_acc')
                    else:
                        raise ValueError('invalid Assembly')
                else:
                    # print(genome, ncbi_assembly, d['data_source'])
                    pass
            else:
                data['data_source'] = g['data_source']
    return graph


class AbstractKeggTransformer(ABC):

    def __init__(self, singles: dict, arrays: dict):
        self.singles = singles
        self.arrays = arrays

    def transform(self, data: dict):
        a, b = data['ENTRY']['value'][0].strip().split(' ', 1)
        a, b.strip()
        doc = {
            'entry': a
        }

        for k, v in self.singles.items():
            if v in data:
                doc[k] = data[v]['value'][0].strip()
        for k, v in self.arrays.items():
            if v in data:
                doc[k] = [x.strip() for x in data[v]['value']]

        return doc


class KeggKOTransformer(AbstractKeggTransformer):

    def __init__(self):
        super().__init__({
            'name': 'NAME',
            'symbol': 'SYMBOL',
        }, {
            'genes': 'GENES',
            'pathway': 'PATHWAY',
            'reaction': 'REACTION',
            'db_links': 'DBLINKS',
        })


class KeggMDTransformer(AbstractKeggTransformer):

    def __init__(self):
        super().__init__({
            'name': 'NAME',
            'class': 'CLASS',
        }, {
            'definition': 'DEFINITION',
            'pathway': 'PATHWAY',
            'reaction': 'REACTION',
            'compound': 'COMPOUND',
            'db_links': 'DBLINKS',
        })


class KeggGenomeTransformer(AbstractKeggTransformer):

    def __init__(self):
        super().__init__({
            'name': 'NAME',
            'taxonomy': 'TAXONOMY',
            'lineage': 'LINEAGE',
            'code': 'ORG_CODE',
            'chromosome': 'CHROMOSOME',
            'category': 'CATEGORY',
            'disease': 'DISEASE',
        }, {
            'data_source': 'DATA_SOURCE',
            'statistics': 'STATISTICS',
            'comments': 'COMMENT',
            'db_links': 'DBLINKS',
        })
