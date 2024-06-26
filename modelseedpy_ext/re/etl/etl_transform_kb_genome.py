import hashlib
import logging
from Bio.Seq import Seq
from modelseedpy import MSGenome
from modelseedpy_ext.re.etl.etl_transform_graph import ETLTransformGraph

logger = logging.getLogger(__name__)


def locate_feature_dna_sequence_in_contig(f, contigs):
    ll = []
    for l in f.data["location"]:
        contig_id, start, direction, seq_size = l
        if contig_id in contigs:
            off_set = seq_size
            sub_str_start = start - 1
            sub_str_end = sub_str_start + off_set
            if direction == "+":
                contig_sub_string = contigs[contig_id][0][
                    sub_str_start : (sub_str_start + seq_size)
                ]
                if contig_sub_string == f.dna_sequence:
                    ll.append(
                        [
                            contig_id,
                            f.id,
                            sub_str_start,
                            sub_str_start + seq_size,
                            False,
                            direction,
                        ]
                    )
                # print(f.id, contig_sub_string[:10], contig_sub_string[-10:])
                # print(f.id, f.dna_sequence[:10], f.dna_sequence[-10:])
            elif direction == "-":
                contig_sub_string = contigs[contig_id][1][
                    (sub_str_start - seq_size + 1) : (sub_str_start + 1)
                ]
                if contig_sub_string == f.dna_sequence[::-1]:
                    ll.append(
                        [
                            contig_id,
                            f.id,
                            sub_str_start - seq_size + 1,
                            sub_str_start + 1,
                            True,
                            direction,
                        ]
                    )
                # print(f.id, contig_sub_string[:10], contig_sub_string[-10:])
                # print(f.id, f.dna_sequence[::-1][:10], f.dna_sequence[::-1][-10:])
            else:
                print("not sure what this means:", direction)

            # print(f.id, l, f.functions, contig_sub_string == f.dna_sequence)
        else:
            print("contig not found:", contig_id)
    return ll


class ETLTransformKBaseGenome:

    def __init__(self, transform_info, transform_contigs, kb_re, kbase):
        self.transform_info = transform_info
        self.transform_contigs = transform_contigs
        self.kb_re = kb_re
        self.kbase = kbase

    def _process_contigs(self, assembly, handle_path='/scratch/fliu/data/kbase/handle/'):
        i = str(assembly.info).replace("/", "_")
        self.kbase.download_file_from_kbase2(
            assembly.fasta_handle_ref, f"{handle_path}/{i}"
        )
        assembly_contigs = MSGenome.from_fasta(f"{handle_path}/{i}")
        contig_set = [x.seq for x in assembly_contigs.features]

        return contig_set

    def _process_features(self, genome, n_genome_id, G):
        seq_dna_hashes = []
        seq_protein_hashes = []
        list_nodes_seq_dna = []
        list_nodes_seq_protein = []
        for f in genome.features:
            seq_dna = f.dna_sequence
            seq_protein = f.protein_translation
            h_dna = self.kb_re.dna_store.store_sequence(seq_dna)
            h_protein = self.kb_re.protein_store.store_sequence(seq_protein)
            seq_dna_hashes.append(h_dna)
            seq_protein_hashes.append(h_protein)

            node_seq_dna = G.add_transform_node(h_dna, "re_seq_dna")
            node_seq_protein = G.add_transform_node(h_protein, "re_seq_protein")

            list_nodes_seq_dna.append(node_seq_dna)
            list_nodes_seq_protein.append(node_seq_protein)

            G.add_transform_edge(n_genome_id, node_seq_dna.id, "ws_genome_has_seq_dna",
                                 {'feature_id': f.id, 'location': f.location})
            G.add_transform_edge(n_genome_id, node_seq_protein.id, "ws_genome_has_seq_protein",
                                 {'feature_id': f.id})

            seq_translation = str(Seq(seq_dna).translate())[:-1]
            if seq_translation == seq_protein:
                G.add_transform_edge(node_seq_dna.id, node_seq_protein.id,
                                     "re_seq_dna_has_translation")
            else:
                G.add_transform_edge(node_seq_dna.id, node_seq_protein.id,
                                     "re_seq_dna_has_partial_translation")

        hash_dna_set = hashlib.sha256(
            "_".join(sorted(seq_dna_hashes)).encode("utf-8")
        ).hexdigest()
        hash_protein_set = hashlib.sha256(
            "_".join(sorted(seq_protein_hashes)).encode("utf-8")
        ).hexdigest()

        node_seq_dna_set = G.add_transform_node(hash_dna_set, "re_seq_dna_set")
        node_seq_protein_set = G.add_transform_node(hash_protein_set, "re_seq_protein_set")

        for n in list_nodes_seq_dna:
            G.add_transform_edge(node_seq_dna_set.id, n.id, "re_seq_dna_set_has_seq_dna")
        for n in list_nodes_seq_protein:
            G.add_transform_edge(node_seq_protein_set.id, n.id, "re_seq_protein_set_has_seq_protein")

        G.add_transform_edge(n_genome_id, node_seq_dna_set.id,
                             "ws_genome_has_seq_dna_set")
        G.add_transform_edge(n_genome_id, node_seq_protein_set.id,
                             "ws_genome_has_seq_protein_set")

    def transform(self, genome):
        assembly = self.kbase.get_from_ws(genome.assembly_ref)

        contig_set = self._process_contigs(assembly)

        g_contigs = self.transform_contigs.transform(contig_set)

        g_genome_info = self.transform_info.transform(genome.info)
        g_assembly_info = self.transform_info.transform(assembly.info)

        n_contig_set_id = list(g_contigs.t_nodes['re_contig_set'])[0]

        n_assembly_id = list(g_assembly_info.t_nodes['kbase_object'])[0]
        n_genome_id = list(g_genome_info.t_nodes['kbase_object'])[0]

        g_compose = g_genome_info.concat(g_assembly_info)
        g_compose = g_compose.concat(g_contigs)
        g_compose.add_transform_edge(n_genome_id, n_assembly_id,
                                     "kbase_object_has_object_reference",
                                     {"value": "assembly_ref"})
        g_compose.add_transform_edge(n_assembly_id, n_contig_set_id,
                                     "kbase_object_has_re_contig_set")

        self._process_features(genome, n_genome_id, g_compose)

        return g_compose


class ETLTransformKBaseGenomeOld(ETLTransformGraph):
    def __init__(self, dna_store, protein_store):
        super().__init__()
        self.dna_store = dna_store
        self.protein_store = protein_store

    def transform(self, kb_genome, contigs):
        nodes = {}
        edges = {}

        def add_node(node_id, label, data=None):
            if label not in nodes:
                nodes[label] = {}
            if label in nodes:
                _node = self.build_node(node_id, label, data)
                if _node.id not in nodes[label]:
                    nodes[label][_node.id] = _node
                else:
                    # print('dup', _node.id)
                    pass

                return nodes[label][_node.id]

            logger.error("add_node error")
            return None

        def add_edge(node_from, node_to, label, data=None):
            if label not in edges:
                edges[label] = []
            if label in edges:
                _edge = self.transform_edge(node_from, node_to, data)
                edges[label].append(_edge)
                return _edge

            logger.error("add_edge error")
            return None

        _contigs = {}
        _contigs_hash = {}
        for contig in contigs.features:
            _contigs[contig.id] = [contig.seq, str(Seq(contig.seq).complement())]
            _contigs_hash[contig.id] = self.dna_store.store_sequence(contig.seq)

        print(len(_contigs))
        print(_contigs_hash)

        node_ws_genome_object = None
        if kb_genome.info:

            key = str(kb_genome.info).replace("/", ":")
            node_ws_genome_object = add_node(
                key,
                "ws_object_version",
                data={
                    "workspace_id": kb_genome.info.workspace_uid,
                    "object_id": kb_genome.info.uid,
                    "version": kb_genome.info.version,
                    "name": kb_genome.info.id,
                },
            )

        seq_dna_hashes = []
        seq_protein_hashes = []
        for f in kb_genome.features:
            seq_dna = f.dna_sequence
            seq_protein = f.protein_translation
            h_dna = self.dna_store.store_sequence(seq_dna)
            h_protein = self.protein_store.store_sequence(seq_protein)
            seq_dna_hashes.append(h_dna)
            seq_protein_hashes.append(h_protein)

            node_seq_dna = add_node(h_dna, "re_seq_dna", data=None)
            node_seq_protein = add_node(h_protein, ""
                                                   "", data=None)

            if node_ws_genome_object:
                add_edge(node_ws_genome_object, node_seq_dna, "ws_genome_has_seq_dna")
                add_edge(
                    node_ws_genome_object, node_seq_protein, "ws_genome_has_seq_protein"
                )

            seq_translation = str(Seq(seq_dna).translate())[:-1]
            if seq_translation == seq_protein:
                add_edge(node_seq_dna, node_seq_protein, "re_seq_dna_has_translation")
            else:
                add_edge(
                    node_seq_dna, node_seq_protein, "re_seq_dna_has_partial_translation"
                )

        hash_dna_set = hashlib.sha256(
            "_".join(sorted(seq_dna_hashes)).encode("utf-8")
        ).hexdigest()
        hash_protein_set = hashlib.sha256(
            "_".join(sorted(seq_protein_hashes)).encode("utf-8")
        ).hexdigest()

        print(hash_dna_set)
        print(hash_protein_set)

        node_seq_dna_set = add_node(hash_dna_set, "re_seq_dna_set", data=None)
        node_seq_protein_set = add_node(
            hash_protein_set, "re_seq_protein_set", data=None
        )

        add_edge(node_ws_genome_object, node_seq_dna_set, "ws_genome_has_seq_dna_set")
        add_edge(
            node_ws_genome_object, node_seq_protein_set, "ws_genome_has_seq_protein_set"
        )

        return nodes, edges


class ETLGenome(ETLTransformGraph):
    def __init__(self, transform_kbase_object, transform_contigs, kbase, handle_path='/home/fliu/kbase/cache/handle'):
        super().__init__()
        self.kbase = kbase
        self.transform_kbase_object = transform_kbase_object
        self.transform_contigs = transform_contigs
        self.handle_path = handle_path

    def etl(self, genome):
        nodes = {}
        edges = {}

        def add_node(node_id, label, data=None):
            if label not in nodes:
                nodes[label] = {}
            if label in nodes:
                _node = self.build_node(node_id, label, data)
                if _node.id not in nodes[label]:
                    nodes[label][_node.id] = _node
                else:
                    # print('dup', _node.id)
                    pass

                return nodes[label][_node.id]

            logger.error("add_node error")
            return None

        def add_edge(node_from, node_to, label, data=None):
            if label not in edges:
                edges[label] = []
            if label in edges:
                _edge = self.transform_edge(node_from, node_to, data)
                edges[label].append(_edge)
                return _edge

            logger.error("add_edge error")
            return None

        assembly = self.kbase.get_from_ws(genome.assembly_ref)

        g_nodes, g_edges = self.transform_kbase_object.transform(genome.info)
        if assembly:
            a_nodes, a_edges = self.transform_kbase_object.transform(assembly.info)
            node_genome = list(a_nodes["kbase_object"].values())[0]
            node_assembly = list(g_nodes["kbase_object"].values())[0]
            add_edge(
                node_genome,
                node_assembly,
                "kbase_object_has_object_reference",
                {"value": "assembly_ref"},
            )

            i = str(assembly.info).replace("/", "_")
            self.kbase.download_file_from_kbase2(
                assembly.fasta_handle_ref, f"{self.handle_path}/{i}"
            )
            assembly_contigs = MSGenome.from_fasta(f"{self.handle_path}/{i}")
            contig_set = [x.seq for x in assembly_contigs.features]

            c_nodes, c_edges = self.transform_contigs.transform(contig_set)

            node_contig_set = list(c_nodes["re_contig_set"].values())[0]

            add_edge(node_contig_set, node_assembly, "re_contig_set_in_object")
