from modelseedpy_ext.re.etl.transform_graph import TransformGraph, Node
from modelseedpy_ext.utils import progress


class ETLNcbi:

    def __init__(self, base_dir):
        self.base_dir = base_dir

    def read_names(self, base_dir=None):
        ncbi_dump_folder = base_dir
        if ncbi_dump_folder is None:
            ncbi_dump_folder = self.base_dir

        node_names = {}
        with open(f'{ncbi_dump_folder}/names.dmp', 'r') as fh:
            line = fh.readline()
            while line:
                tax_id, name, unique_name, name_class = line.strip()[:-2].split('\t|\t')
                if tax_id not in node_names:
                    node_names[tax_id] = {}
                if name_class not in node_names[tax_id]:
                    node_names[tax_id][name_class] = []
                node_names[tax_id][name_class].append([name, unique_name])
                line = fh.readline()
        return node_names

    def read_div_nodes(self, graph, base_dir=None):
        ncbi_dump_folder = base_dir
        if ncbi_dump_folder is None:
            ncbi_dump_folder = self.base_dir

        div_nodes = {}
        with open(f'{ncbi_dump_folder}/division.dmp', 'r') as fh:
            line = fh.readline()
            while line:
                div_id, code, name, comments = line.strip()[:-2].split('\t|\t')
                node_div = graph.add_transform_node2(
                    Node(div_id, 'ncbi_taxonomy_division', data={'name': name, 'code': code}))
                if node_div._key not in div_nodes:
                    div_nodes[node_div._key] = node_div
                line = fh.readline()

        return div_nodes

    def read_code_nodes(self, graph, base_dir=None):
        ncbi_dump_folder = base_dir
        if ncbi_dump_folder is None:
            ncbi_dump_folder = self.base_dir

        genetic_code_nodes = {}
        with open(f'{ncbi_dump_folder}/gencode.dmp', 'r') as fh:
            line = fh.readline()
            while line:
                gen_id, abbreviation, name, cde, starts = line.strip()[:-2].split('\t|\t')
                gen_data = {
                    'name': name,
                    'cde': cde,
                    'starts': starts
                }
                if abbreviation:
                    gen_data['abbreviation'] = abbreviation
                node_gen = graph.add_transform_node2(Node(gen_id, 'ncbi_taxonomy_genetic_code', data=gen_data))
                if node_gen._key not in genetic_code_nodes:
                    genetic_code_nodes[node_gen._key] = node_gen
                line = fh.readline()

        return genetic_code_nodes

    def read_taxa_nodes(self, graph, node_names, div_nodes, genetic_code_nodes, base_dir=None):
        ncbi_dump_folder = base_dir
        if ncbi_dump_folder is None:
            ncbi_dump_folder = self.base_dir

        taxa_nodes = {}
        with open(f'{ncbi_dump_folder}/nodes.dmp', 'r') as fh:
            line = fh.readline()
            while line:
                tax_id, parent_tax_id, rank, embl_code, division_id, inherited_div, genetic_code_id, \
                  inherited_genetic_code, mitochondrial_genetic_code_id, inherited_mito_genetic_code, \
                  gb_hidden, hidden_subtree_root, comments = line.strip()[:-2].split('\t|\t')

                data = {
                    'parent_tax_id': parent_tax_id,
                    'rank': rank,
                    'inherited_div': inherited_div == '1',
                    'inherited_genetic_code': inherited_genetic_code == '1',
                    'inherited_mito_genetic_code': inherited_mito_genetic_code == '1',
                    'genbank_hidden': gb_hidden == '1',
                    'hidden_subtree_root': hidden_subtree_root == '1',
                }
                if embl_code:
                    data['embl_code'] = embl_code
                if comments:
                    data['comments'] = comments
                data.update(node_names[tax_id])
                node_taxa = graph.add_transform_node2(Node(tax_id, 'ncbi_taxonomy', data=data))
                if node_taxa._key not in taxa_nodes:
                    taxa_nodes[node_taxa._key] = node_taxa

                node_taxa_div = div_nodes[division_id]
                graph.add_transform_edge2(node_taxa, node_taxa_div, 'ncbi_taxonomy_has_ncbi_taxonomy_division')

                node_taxa_gen = genetic_code_nodes[genetic_code_id]
                graph.add_transform_edge2(node_taxa, node_taxa_gen, 'ncbi_taxonomy_has_ncbi_taxonomy_genetic_code')

                node_taxa_gen_mito = genetic_code_nodes[mitochondrial_genetic_code_id]
                graph.add_transform_edge2(node_taxa, node_taxa_gen_mito,
                                          'ncbi_taxonomy_has_mitochondrial_ncbi_taxonomy_genetic_code')

                line = fh.readline()

        return taxa_nodes

    @staticmethod
    def build_taxa_link(graph, taxa_nodes):
        for node_id, node in progress(taxa_nodes.items()):
            if node_id == node.data['parent_tax_id']:
                pass
            else:
                other_node = taxa_nodes[node.data['parent_tax_id']]
                graph.add_transform_edge2(node, other_node, 'ncbi_taxonomy_has_parent_ncbi_taxonomy')

    def etl(self):
        graph = TransformGraph()
        node_names = self.read_names()
        div_nodes = self.read_div_nodes(graph)
        genetic_code_nodes = self.read_code_nodes(graph)
        taxa_nodes = self.read_taxa_nodes(graph, node_names, div_nodes, genetic_code_nodes)
        self.build_taxa_link(graph, taxa_nodes)

        return graph
