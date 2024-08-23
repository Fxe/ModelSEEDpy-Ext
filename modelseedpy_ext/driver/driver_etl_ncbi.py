from modelseedpy_ext.re.etl.transform_graph import TransformGraph, Node
import os
import json


class DriverEtlNcbi:

    def __init__(self, re, extract_ncbi, assembly_cache_path='M:/biodb/ncbi/cache/assembly'):
        self.assembly_cache_path = assembly_cache_path
        self.re = re
        self.extract_ncbi = extract_ncbi

    def etl_ncbi_assembly(self, ncbi_ids):
        # scan gca
        gca_to_assembly = self.re.get_link("ncbi_assembly_has_ncbi_assembly_gca_acc",
                                           ['ncbi_assembly_gca_acc/' + x for x in ncbi_ids], rev=True)
        # scan gcf
        gcf_to_assembly = self.re.get_link("ncbi_assembly_has_ncbi_assembly_gcf_acc",
                                           ['ncbi_assembly_gcf_acc/' + x for x in ncbi_ids], rev=True)

        # get missing
        missing = ncbi_ids - set(gca_to_assembly.keys()) - set(gcf_to_assembly.keys())

        # run api for ids
        missing_ncbi_assembly_ids = set()
        for i in missing:
            res = self.extract_ncbi.esearch('assembly', i.split('/')[1])
            if res and 'ids' in res:
                missing_ncbi_assembly_ids |= res.get('ids', set())

        # run api for data
        collected_files = set()
        for assembly_id in missing_ncbi_assembly_ids:
            filename = f'{self.assembly_cache_path}/{assembly_id}.json'
            if not os.path.exists(filename):
                res = self.extract_ncbi.esummary('assembly', assembly_id)[0][0]
                if res and len(res) > 0:
                    res['_id'] = assembly_id
                    with open(filename, 'w') as fh:
                        fh.write(json.dumps(res))
                        collected_files.add(filename)
            else:
                collected_files.add(filename)
        # elt data
        data = {}
        for f in collected_files:
            with open(f, 'r') as fh:
                _d = json.load(fh)
                data[_d['_id']] = _d

        self._etl_assemblies(data)

    @staticmethod
    def _etl_assemblies(ncbi_assemblies):
        graph = TransformGraph()
        for i, d in ncbi_assemblies.items():
            tax_id = d['Taxid']  # strain or other
            tax_species_id = d['SpeciesTaxid']  # species level
            bio_sample_id = d['BioSampleId']
            synonym = d['Synonym']
            ncbi_gb_acc = synonym.get('Genbank')
            ncbi_rs_acc = synonym.get('RefSeq')
            ncbi_gb_rs_type = synonym.get('Similarity')
            # print(ncbi_gb_rs_type, ncbi_gb_acc, ncbi_rs_acc)

            node_assembly = graph.add_transform_node2(Node(i, 'ncbi_assembly', data=d))
            node_gb_acc = None
            node_rs_acc = None

            if tax_id:
                node_taxa = graph.add_transform_node2(Node(str(tax_id), 'ncbi_taxonomy', data={'_proxy': True}))
                graph.add_transform_edge2(node_assembly, node_taxa, 'ncbi_assembly_has_ncbi_taxonomy')
            if tax_species_id:
                node_taxa_species = graph.add_transform_node2(Node(str(tax_species_id), 'ncbi_taxonomy', data={'_proxy': True}))
                graph.add_transform_edge2(node_assembly, node_taxa_species, 'ncbi_assembly_has_ncbi_taxonomy_species')

            if ncbi_gb_acc is not None:
                node_gb_acc = graph.add_transform_node2(Node(ncbi_gb_acc, 'ncbi_assembly_gca_acc'))
                graph.add_transform_edge2(node_assembly, node_gb_acc, 'ncbi_assembly_has_ncbi_assembly_gca_acc')
            if ncbi_rs_acc is not None:
                node_rs_acc = graph.add_transform_node2(Node(ncbi_rs_acc, 'ncbi_assembly_gcf_acc'))
                graph.add_transform_edge2(node_assembly, node_rs_acc, 'ncbi_assembly_has_ncbi_assembly_gcf_acc')
            if node_gb_acc is not None and node_rs_acc is not None:
                graph.add_transform_edge2(node_gb_acc, node_rs_acc, 'ncbi_assembly_gca_to_gcf')
                # tax_doc = re.db['ncbi_taxonomy'][tax_id]
                # tax_species_doc = re.db['ncbi_taxonomy'][tax_species_id]
            # print(i, tax_id, tax_species_id, tax_doc['rank'], tax_doc['scientific name'][0][0], tax_species_doc['rank'], tax_species_doc['scientific name'][0][0])
        return graph