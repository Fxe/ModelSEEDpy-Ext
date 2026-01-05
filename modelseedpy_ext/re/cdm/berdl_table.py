from modelseedpy_ext.re.core.genome import ProteinSequence


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


class Wut:

    def __init__(self):
        self.convert_location = None
        self.collect_ontology = None
        pass

    def wut(self, genome, m_to_c, col_bakta, data_psortb, genome_scale_model,
            representative_to_pan_cluster_id, cluster_core, g_all_proteins):
        data_features = []
        for f in genome.features:
            locations = f.location
            if len(locations) > 1:
                raise ValueError('!')
            contig, start, end, strand = self.convert_location(f.location[0])

            annotation = {}

            if f.seq:
                protein = ProteinSequence(f.seq)
                protein_hash = protein.hash_value
                mmseqs_cluster = m_to_c[protein_hash]
                annotation_bakta = col_bakta.find_one({'_id': protein_hash})
                ontology_bakta = []
                if annotation_bakta:
                    ontology_bakta = self.collect_ontology(annotation_bakta)
                for kind, value in ontology_bakta:
                    if kind not in annotation:
                        annotation[kind] = []
                    annotation[kind].append(value)

                batka_gene_names = annotation_bakta.get('genes')
                if batka_gene_names:
                    annotation['gene_names'] = batka_gene_names
                pan_genome_cluser_from_mmseqs_match = representative_to_pan_cluster_id.get(mmseqs_cluster)
                pan_genome_cluser_from_mmseqs = None
                is_core = None
                if pan_genome_cluser_from_mmseqs_match is not None:
                    flag_is_core = [x in cluster_core for x in pan_genome_cluser_from_mmseqs_match]
                    if len(pan_genome_cluser_from_mmseqs_match) == 1:
                        pan_genome_cluser_from_mmseqs = list(pan_genome_cluser_from_mmseqs_match)[0]
                        is_core = flag_is_core[0]
                    else:
                        pan_genome_cluser_from_mmseqs = '; '.join(pan_genome_cluser_from_mmseqs_match)
                        if set(flag_is_core) == 1:
                            is_core = flag_is_core[0]

                rast_annotation = g_all_proteins.features.get_by_id(protein_hash).ontology_terms.get('RAST', [])
                if len(rast_annotation) > 0:
                    annotation['rast'] = '\t'.join(rast_annotation)

                annotation['ke_pan_genome_is_core'] = is_core

                annotation['seq_aa'] = f.seq

                if pan_genome_cluser_from_mmseqs:
                    annotation['ke_pan_genome_cluster_id'] = pan_genome_cluser_from_mmseqs

                if protein_hash in data_psortb:
                    d = data_psortb[protein_hash]
                    annotation['psortb'] = d['Final_Localization']

                model_reactions = {}
                if f.id in genome_scale_model.genes:
                    reactions = genome_scale_model.genes.get_by_id(f.id).reactions
                    for r in reactions:
                        model_reactions[r.id] = r

                annotation['seed.reaction'] = {r.id: r.build_reaction_string() for r in model_reactions.values()}
            data_features.append(
                RowFeature(genome.info.id, contig, f.id, len(f.seq), start, end, strand, annotation))
