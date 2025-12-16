import polars as pl
from sqlalchemy import text
from sqlalchemy.orm import Session
import json
from Bio import SeqIO
import pandas as pd


class AlexeyIdentifierFetch:

    def __init__(self, engine, genome_name, gbff_records):
        self.engine = engine
        self.genome_name = genome_name
        self.alexey_genome_id, self.alexey_genome, self.df_alexey_contigs = self.sql_search_alexey_data(genome_name)
        self.d_alexey_contig = self.df_alexey_contigs.set_index('contig_id').transpose().to_dict()
        self.pairs = AlexeyIdentifierFetch.pair_gbff_to_alexey(gbff_records, self.d_alexey_contig)

        self.df_alexey_gene = pd.read_sql(f"SELECT * FROM browser_gene WHERE genome_id = {self.alexey_genome_id}",
                                          self.engine).set_index('id')
        self.contig_ids = {str(k): None for k in self.df_alexey_gene['contig_id']}

        self.contig_id_to_d = {}
        for feature_id, d in self.df_alexey_gene.iterrows():
            if d['contig_id'] not in self.contig_id_to_d:
                self.contig_id_to_d[d['contig_id']] = {}
            if d['type'] not in self.contig_id_to_d[d['contig_id']]:
                self.contig_id_to_d[d['contig_id']][d['type']] = {}
            if d['start'] not in self.contig_id_to_d[d['contig_id']][d['type']]:
                self.contig_id_to_d[d['contig_id']][d['type']][d['start']] = {}
            self.contig_id_to_d[d['contig_id']][d['type']][d['start']][feature_id] = d

    def sql_search_alexey_data(self, genome_name):
        df_alexey_genome = pd.read_sql(f"SELECT * FROM `browser_genome` WHERE name = '{genome_name}'", self.engine)
        if len(df_alexey_genome) != 1:
            raise ValueError(f'expected one and only one record. Returned records: {len(df_alexey_genome)}')
        alexey_genome = next(df_alexey_genome.iterrows())[1]
        alexey_genome_id = alexey_genome['id']

        df_alexey_contigs = pd.read_sql(f"SELECT * FROM browser_contig WHERE genome_id = {alexey_genome_id}",
                                        self.engine)

        return alexey_genome_id, alexey_genome, df_alexey_contigs

    def get_contigset_identifier(self, cdm_contigset):
        return cdm_contigset.hash_contigset, self.alexey_genome_id, 'contigset', 'GenomeDepot'

    def get_contigset_name(self, cdm_contigset):
        return cdm_contigset.hash_contigset, self.genome_name, 'contigset', 'GenomeDepot'

    @staticmethod
    def pair_gbff_to_alexey(gbff_records, d_alexey_contig):
        gff_record_alexey_pair = set()
        for gbff_contig in gbff_records:
            accession = gbff_contig.annotations.get('accessions', [None])[0]

            pair = None
            if gbff_contig.name in d_alexey_contig:
                pair = (gbff_contig.id, gbff_contig.name)
            elif gbff_contig.id in d_alexey_contig:
                pair = (gbff_contig.id, gbff_contig.id)
            elif accession and accession in d_alexey_contig:
                pair = (gbff_contig.id, accession)
            else:
                raise ValueError('unable to pair gbff contig')
            if pair and pair not in gff_record_alexey_pair:
                gff_record_alexey_pair.add(pair)
        return gff_record_alexey_pair


class EE:

    def __init__(self):
        self._data_contigset = {}
        self._data_contig = {}
        self._data_contigset_x_contig = {}

    def read_gbff_records_from_file(self, filename: str):
        if filename.endswith(".gbff"):
            with open(filename, 'r') as fh:
                return self.read_gbff_records(fh)
        elif filename.endswith(".gz"):
            import gzip
            from io import StringIO
            with gzip.open(filename, 'rb') as fh:
                return self.read_gbff_records(StringIO(fh.read().decode('utf-8')))

    def read_gbff_records(self, handler):
        gbff_records = []
        for record in SeqIO.parse(handler, "gb"):
            gbff_records.append(record)
        return gbff_records

    def etl(self, filename):

        gbff_records = self.read_gbff_records_from_file(filename)

        from modelseedpy_ext.re.hash_seq import HashSeqList, HashSeq

        num_features = 0
        hlist_contigs = HashSeqList()
        contig_len = 0
        for contig in gbff_records:
            features = [x for x in contig.features if x.type == 'gene']
            num_features += len(features)
            hlist_contigs.append(HashSeq(str(contig.seq)))
            contig_len += len(contig.seq)
            # print(contig.id, contig.name, len(contig.features))
        print(num_features)

        cdm_contigset = CDMContigSet(hlist_contigs.hash_value)
        cdm_contigset.size = contig_len

        if cdm_contigset.hash_contigset not in self._data_contigset:
            self._data_contigset[cdm_contigset.hash_contigset] = cdm_contigset
            self._data_contigset[cdm_contigset.hash_contigset].names = [_contigset_name]
            self._data_contigset[cdm_contigset.hash_contigset].identifiers = [_contigset_identifier]
        else:
            self._data_contigset[cdm_contigset.hash_contigset].names.append(_contigset_name)
            self._data_contigset[cdm_contigset.hash_contigset].identifiers.append(_contigset_identifier)

        for gbff_id, alexey_id in pairs:
            gbff_record = name_to_gbff_contig[gbff_id]
            d_alexey_conti_data = d_alexey_contig[alexey_id]

            cdm_contig = CDMContig(cdm_contigset.hash_contigset, str(gbff_record.seq))
            if cdm_contig.hash_contig not in _data_contig:
                _data_contig[cdm_contig.hash_contig] = cdm_contig
            pk = (cdm_contigset.hash_contigset, cdm_contig.hash_contig)
            _id = f'{hash_contigset}_{hash_contig}_{i}'
            _id_hash = _hash_string(_id)

            unique_names = set()
            unique_names.add((d_alexey_conti_data['name'], 'GenomeDepot'))
            unique_names.add((record_name, 'GenomeDepot'))
            unique_names.add((gbff_record.id, 'gbff'))
            unique_names.add((gbff_record.name, 'gbff'))
            names = [
                (cdm_contigset.hash_contigset, cdm_contig.hash_contig, x[0], 'contig', x[1]) for x in unique_names
            ]
            identifiers = [(cdm_contigset.hash_contigset, cdm_contig.hash_contig, d_alexey_conti_data['id'], 'contig',
                            'GenomeDepot')]

            if pk not in _data_contigset_x_contig:
                _data_contigset_x_contig[pk] = {'count': 1, 'names': [names], 'identifiers': [identifiers]}
            else:
                _data_contigset_x_contig[pk]['count'] += 1
                _data_contigset_x_contig[pk]['names'].append(names)
                _data_contigset_x_contig[pk]['identifiers'].append(identifiers)

    def process_features(self, genome_name, gbff_records, protocol_id):

        param = ', '.join([f"'{x}'" for x in contig_ids])
        df_response = pd.read_sql(
            f"SELECT id_entity, identifier FROM identifier WHERE name_table = 'contig' AND identifier IN ({param})",
            engine_cdm)

        for row_id, d in df_response.iterrows():
            if d['identifier'] in contig_ids:
                contig_ids[d['identifier']] = d['id_entity']
        assert len(list(filter(lambda o: o[1] is None, contig_ids.items()))) == 0, 'missing contig map'

        #assembly = REAssembly.from_fasta(_scan_data['path_contigs.fasta'])
        #hash_contigset = assembly.hash_value

        #gbff_records = []
        #with open(_scan_data['path_gbff'], 'r') as fh:
        #    for record in SeqIO.parse(fh, "gb"):
        #        gbff_records.append(record)

        name_to_gbff_contig = {}
        for gbff_contig in gbff_records:
            if gbff_contig.id not in name_to_gbff_contig:
                name_to_gbff_contig[gbff_contig.id] = gbff_contig
            else:
                raise ValueError('gbff record id not unique')
        pairs = pair_gbff_to_alexey(gbff_records, d_alexey_contig)

        feature_registry = {}
        feature_to_protein = {}
        for gbff_id, alexey_id in pairs:
            gbff_record = name_to_gbff_contig[gbff_id]
            d_alexey_conti_data = d_alexey_contig[alexey_id]

            d_alexey_lib = contig_id_to_d.get(d_alexey_conti_data['id'], {})

            cdm_contig = CDMContig(hash_contigset, str(gbff_record.seq))
            hash_contig = cdm_contig.hash_contig

            # print(hash_contigset, hash_contig)

            x = pd.read_sql(
                f"SELECT hash_contigset_x_contig FROM contigset_x_contig WHERE hash_contigset = '{hash_contigset}' AND hash_contig = '{hash_contig}'",
                engine_cdm)
            if len(x) != 1:
                raise ValueError('contigset_x_contig error')
            contigset_x_contig_id = list(x['hash_contigset_x_contig'])[0]
            # print(contigset_x_contig_id)

            register_features(contigset_x_contig_id, gbff_record, d_alexey_lib, feature_registry, feature_to_protein,
                              protocol_id)

        # check if new features do not exist in catalog

        found = set(feature_registry) & set(_data_feature)
        if len(found) != 0:
            raise ValueError(f'{genome_name}: unable to register features {found}')

        for feature_id in feature_registry:
            if feature_id not in _data_feature:
                _data_feature[feature_id] = feature_registry[feature_id]

        for feature_id in feature_to_protein:
            cdm_protein = feature_to_protein[feature_id]
            _data_feature_x_protein[feature_id] = cdm_protein.hash
            # if cdm_protein.hash not in all_hashes:
            #    print(f'{genome_name}: {feature_id} {cdm_protein.hash} not found')
            _data_protein[cdm_protein.hash] = cdm_protein

        return None






    def export(self):
        _table_name_names = ["id_entity", "name_table", "name", "source", "description"]
        _table_name_data = {k: [] for k in _table_name_names}

        _table_id_names = ["id_entity", "name_table", "identifier", "source", "description"]
        _table_id_data = {k: [] for k in _table_id_names}

        for id_entity, name_table, name, source, description in tqdm(self._data_feature_name):
            _table_name_data['id_entity'].append(id_entity)
            _table_name_data['name_table'].append(name_table)
            _table_name_data['name'].append(str(name))
            _table_name_data['source'].append(source)
            _table_name_data['description'].append(description)

        for id_entity, name_table, identifier, source, description in tqdm(self._data_feature_identifier):
            _table_id_data['id_entity'].append(id_entity)
            _table_id_data['name_table'].append(name_table)
            _table_id_data['identifier'].append(str(identifier))
            _table_id_data['source'].append(source)
            _table_id_data['description'].append(description)
        table_name = pa.Table.from_arrays(
            [pa.array(_table_name_data[k]) for k in _table_name_names],
            names=_table_name_names)
        table_identifier = pa.Table.from_arrays(
            [pa.array(_table_id_data[k]) for k in _table_id_names],
            names=_table_id_names)

        return table_name, table_identifier

    def export_feature_table(self):
        _table_feature_names = ["feature_id", "hash_contigset_x_contig", "source", "protocol_id", "type",
                                "start", "end", "strand", "cds_phase"]
        _table_feature_data = {k: [] for k in _table_feature_names}

        for feature_id, cdm_feature in tqdm(_data_feature.items()):
            contigset_x_contig_id = cdm_feature.hash_contig
            _table_feature_data['feature_id'].append(cdm_feature.id)
            _table_feature_data['hash_contigset_x_contig'].append(contigset_x_contig_id)
            _table_feature_data['source'].append('genomedepot')
            _table_feature_data['protocol_id'].append(cdm_feature.protocol)
            _table_feature_data['type'].append(cdm_feature.type)
            _table_feature_data['start'].append(cdm_feature.start)
            _table_feature_data['end'].append(cdm_feature.end)
            _table_feature_data['strand'].append(cdm_feature.strand)
            _table_feature_data['cds_phase'].append(cdm_feature.cds_phase)

        table_feature = pa.Table.from_arrays(
            [pa.array(_table_feature_data[k]) for k in _table_feature_names],
            names=_table_feature_names)

        return table_feature

    def export_feature_to_protein_table(self):
        _table_feature_x_protein_names = ["feature_id", "hash_protein_sequence", "has_stop_codon"]
        _table_feature_x_protein_data = {k: [] for k in _table_feature_x_protein_names}
        for feature_id, hash_protein_seq in tqdm(_data_feature_x_protein.items()):
            _table_feature_x_protein_data['feature_id'].append(feature_id)
            _table_feature_x_protein_data['hash_protein_sequence'].append(hash_protein_seq)
            _table_feature_x_protein_data['has_stop_codon'].append(True)
        table_feature_x_protein = pa.Table.from_arrays(
            [pa.array(_table_feature_x_protein_data[k]) for k in _table_feature_x_protein_names],
            names=_table_feature_x_protein_names)
        return table_feature_x_protein



class ETLLoadMysql:

    def __init__(self, engine):
        self.engine = engine
        self.df_stream = False
        self.load_config = {
            'contigset': ("SELECT hash_contigset FROM contigset;", lambda x: x[0]),
            'contig': ("SELECT hash_contigset, hash_contig FROM contig;", lambda x: (x[0], x[1])),
            'feature': ("SELECT feature_id FROM feature;", lambda x: x[0]),
            'protein': ("SELECT hash_protein_sequence FROM protein;", lambda x: x[0]),
            'encoded_feature': ("SELECT feature_id, hash_protein_sequence FROM encoded_feature;",
                                lambda x: (x[0], x[1])),
            'protocol': ("SELECT protocol_id FROM protocol;", lambda x: x[0]),
            'cluster': ("SELECT cluster_id FROM clusters;", lambda x: x[0]),
            'cluster_x_protein': ("SELECT cluster_id, hash_protein_sequence, feature_id FROM cluster_x_protein;",
                                  lambda x: (x[0], x[1], x[2])),

            'ani': ("SELECT hash_contigset_1, hash_contigset_2, protocol_id FROM ani;", lambda x: (x[0], x[1], x[2])),
            'rep': (),

            'name': (),
            'identifier': (),

            'y_rast': ("SELECT id_rast FROM y_rast;", lambda x: x[0]),
            'y_protein_x_rast': ("SELECT hash_protein_sequence, id_rast FROM y_protein_x_rast;", lambda x: (x[0], x[1])),

        }

    @staticmethod
    def get_str(s):
        if s is None:
            return 'NULL'
        else:
            return f"'{s}'"

    def fetch_keys(self, query, collect_key):
        def fn_fetch_keys():
            table_ids = set()
            with self.engine.connect() as con:
                rs = con.execute(text(query))
                for o in rs:
                    table_ids.add(collect_key(o))
            return table_ids

        return fn_fetch_keys

    def load_contigset(self, lazy_frame, pos=0, batch_size=20000):
        def build_query(row):
            val_hash_contigset = row[0]
            val_n_contigs = row[1]
            pk = val_hash_contigset
            sql_insert = f"""
            INSERT INTO contigset (hash_contigset, n_contigs) 
            VALUES ('{val_hash_contigset}', {val_n_contigs});
            """
            return pk, sql_insert

        query, collect_key = self.load_config['contigset']

        self.load_lazy_frame(lazy_frame, self.fetch_keys(query, collect_key), build_query, pos, batch_size)

    def load_contig(self, lazy_frame, pos=0, batch_size=20000):
        def build_query(row):
            hash_contigset, hash_contig, length, gc_content, base_count = row
            pk = (hash_contigset, hash_contig)
            json_str_base_count = json.dumps({k: v for k, v in base_count.items() if v is not None})
            sql_insert = f"""
            INSERT INTO `contig` (`hash_contigset`, `hash_contig`, `length`, `gc_content`, `base_count`) 
            VALUES ("{hash_contigset}", "{hash_contig}", {length}, {gc_content}, '{json_str_base_count}');
            """
            return pk, sql_insert

        query, collect_key = self.load_config['contig']

        self.load_lazy_frame(lazy_frame, self.fetch_keys(query, collect_key), build_query, pos, batch_size)

    def load_feature(self, lazy_frame, pos=0, batch_size=20000):
        def build_query(row):
            feature_id, hash_contig, hash_contigset, source_database, source_algorithm, \
            feature_type, start, end, strand, cds_phase, attributes = row
            pk = feature_id
            sql_insert = f"""
            INSERT 
                INTO `feature` (`feature_id`, `hash_contig`, `hash_contigset`, 
                `source_database`, `source_algorithm`, `type`, 
                `start`, `end`, `strand`, `cds_phase`)
                VALUES ("{feature_id}", "{hash_contig}", "{hash_contigset}", 
                "{source_database}", "{source_algorithm}", "{feature_type}", 
                {start}, {end}, '{strand}', NULL);
            """
            return pk, sql_insert

        query, collect_key = self.load_config['contigset']

        self.load_lazy_frame(lazy_frame, self.fetch_keys(query, collect_key), build_query, pos, batch_size)

    def load_protein(self, lazy_frame, pos=0, batch_size=20000):
        def build_query(row):
            hash_protein_sequence, description, length, sequence, evidence_for_existence = row
            pk = hash_protein_sequence
            sql_insert = f"""
            INSERT INTO protein (hash_protein_sequence, length, sequence) 
            VALUES ('{hash_protein_sequence}', {length}, '{sequence}');
            """
            return pk, sql_insert

        query, collect_key = self.load_config['protein']

        self.load_lazy_frame(lazy_frame, self.fetch_keys(query, collect_key), build_query, pos, batch_size)

    def load_encoded_feature(self, lazy_frame, pos=0, batch_size=20000):
        def build_query(row):
            feature_id, hash_protein_sequence, has_stop_codon = row
            pk = (feature_id, hash_protein_sequence)
            sql_insert = f"""
            INSERT INTO encoded_feature (feature_id, hash_protein_sequence, has_stop_codon) 
            VALUES ('{feature_id}', '{hash_protein_sequence}', {has_stop_codon});
            """
            return pk, sql_insert

        query, collect_key = self.load_config['encoded_feature']

        self.load_lazy_frame(lazy_frame, self.fetch_keys(query, collect_key), build_query, pos, batch_size)

    def load_protocol(self, lazy_frame, pos=0, batch_size=20000):
        def build_query(row):
            protocol_id = row[0]
            parent_protocol_id = self.get_str(row[1])
            name = self.get_str(row[2])
            pk = protocol_id
            sql_insert = f"""
            INSERT INTO protocol (protocol_id, parent_protocol_id, name) 
            VALUES ('{protocol_id}', {parent_protocol_id}, {name});
            """
            return pk, sql_insert

        query, collect_key = self.load_config['protocol']

        self.load_lazy_frame(lazy_frame, self.fetch_keys(query, collect_key), build_query, pos, batch_size)

    def load_cluster(self, lazy_frame, pos=0, batch_size=20000):
        def build_query(row):
            cluster_id, protocol_id, cluster_annotation_id, description, number_of_candidates, \
                cluster_annotation, uniprot_identity = row
            pk = cluster_id
            sql_insert = f"""
            INSERT INTO clusters (cluster_id, protocol_id , description, number_of_candidates) 
            VALUES ('{cluster_id}', '{protocol_id}', '{description}', {number_of_candidates});
            """
            return pk, sql_insert

        query, collect_key = self.load_config['cluster']

        self.load_lazy_frame(lazy_frame, self.fetch_keys(query, collect_key), build_query, pos, batch_size)

    def load_cluster_x_protein(self, lazy_frame, pos=0, batch_size=20000):
        def build_query(row):
            cluster_id, hash_protein_sequence, feature_id = row
            pk = (cluster_id, hash_protein_sequence, feature_id)
            sql_insert = f"""
            INSERT INTO cluster_x_protein (cluster_id, hash_protein_sequence, feature_id) 
            VALUES ('{cluster_id}', '{hash_protein_sequence}', '{feature_id}');
            """
            return pk, sql_insert

        query, collect_key = self.load_config['cluster_x_protein']

        self.load_lazy_frame(lazy_frame, self.fetch_keys(query, collect_key), build_query, pos, batch_size)

    def load_ani(self, lazy_frame, pos=0, batch_size=20000, protocol_id='default'):
        def build_query(row):
            hash_contigset_1, hash_contigset_2, ani, af, af_mapped, af_total = row
            pk = (hash_contigset_1, hash_contigset_2, protocol_id)
            sql_insert = f"""
            INSERT INTO ani (hash_contigset_1, hash_contigset_2, protocol_id, ani)
            VALUES ('{hash_contigset_1}', '{hash_contigset_2}', '{protocol_id}', {ani});
            """
            return pk, sql_insert

        query, collect_key = self.load_config['ani']

        self.load_lazy_frame(lazy_frame, self.fetch_keys(query, collect_key), build_query, pos, batch_size)

    def load_y_rast(self, lazy_frame, pos=0, batch_size=20000):
        def build_query(row):
            id_rast, annotation_rast = row
            annotation_rast = annotation_rast.replace("'", "\\'")
            pk = id_rast
            sql_insert = f"""
            INSERT INTO y_rast (id_rast, annotation_rast) 
            VALUES ('{id_rast}', '{annotation_rast}');
            """
            return pk, sql_insert

        query, collect_key = self.load_config['y_rast']

        self.load_lazy_frame(lazy_frame, self.fetch_keys(query, collect_key), build_query, pos, batch_size)

    def load_y_protein_x_rast(self, lazy_frame, pos=0, batch_size=20000):
        def build_query(row):
            hash_protein_sequence, id_rast = row
            pk = (hash_protein_sequence, id_rast)
            sql_insert = f"""
            INSERT INTO y_protein_x_rast (hash_protein_sequence, id_rast) 
            VALUES ('{hash_protein_sequence}', '{id_rast}');
            """
            return pk, sql_insert

        query, collect_key = self.load_config['y_protein_x_rast']

        self.load_lazy_frame(lazy_frame, self.fetch_keys(query, collect_key), build_query, pos, batch_size)

    def load_lazy_frame(self, lazy_frame, fn_fetch_keys, fn_build_query, pos=0, batch_size=20000):
        """

        :param lazy_frame: polars lazy dataframe
        :param fn_fetch_keys: fetch existing keys from table
        :param fn_build_query: convert slice row to primary key and sql insert query
        :param pos: starting slice position
        :param batch_size: slice size
        :return:
        """
        df_size = lazy_frame.select(pl.len()).collect()['len'][0]
        table_ids = fn_fetch_keys()
        print(df_size, len(table_ids))
        while pos < df_size:
            with Session(self.engine) as session:
                batch_ids = {}
                for row in lazy_frame.slice(pos, batch_size).collect(streaming=self.df_stream).iter_rows():
                    pk, sql_insert = fn_build_query(row)
                    if pk not in table_ids and pk not in batch_ids:
                        batch_ids[pk] = row
                        session.execute(text(sql_insert))
                    elif pk in batch_ids:
                        print(pk, row, batch_ids[pk])
                session.commit()
            pos += batch_size
            print(pos)
