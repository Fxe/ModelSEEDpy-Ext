import logging
from modelseedpy_ext.re.etl.etl_transform_graph import ETLTransformGraph

logger = logging.getLogger(__name__)


class ETLTransformUniref(ETLTransformGraph):
    """
    This transform turns a uniref record to a graph, because it requires extract/transform both uniprot and uniparc
    the method is super unstable and the code is rushed because I have better things to do rather than load 
    1000 database. Will improve in the far far future :)
    """

    def __init__(self, protein_store, uniref_collection_name,
                 etl_extract_uniparc, etl_transform_uniparc,
                 etl_extract_uniprot, etl_transform_uniprot,
                 etl_load, database):
        super().__init__()
        self.database = database
        self.uniref_collection = uniref_collection_name
        self.protein_store = protein_store
        self.etl_extract_uniparc = etl_extract_uniparc
        self.etl_extract_uniprot = etl_extract_uniprot
        self.etl_transform_uniparc = etl_transform_uniparc
        self.etl_transform_uniprot = etl_transform_uniprot
        self.etl_load = etl_load
        self.re_seq_protein_collection = 're_seq_protein'
        self.uniref_has_member = f'{self.uniref_collection}_has_member'
        self.aql_uniparc = """
        FOR edge IN uniparc_has_protein_sequence
            FILTER edge._from IN @l
            RETURN edge
        """
        self.aql_uniprot = """
        FOR acc IN @acc_list
            LET sprot = (
                LET f = (
                FOR q IN uniprotkb_sprot_has_accession
                    FILTER q._to == acc
                    RETURN q._from
                )
                FOR edge IN uniprotkb_sprot_has_protein_sequence
                    FILTER edge._from IN f
                    RETURN edge._to
            )

            LET trembl = (
                LET f = (
                FOR q IN uniprotkb_trembl_has_accession
                    FILTER q._to == acc
                    RETURN q._from
                )
                FOR edge IN uniprotkb_trembl_has_protein_sequence
                    FILTER edge._from IN f
                    RETURN edge._to
            )


            RETURN {'acc': acc, 'sprot': sprot, 'trembl': trembl}
        """

    def get_uniparc_hash(self, fetch_uniparc, nodes, fetch_missing=False):
        query = self.database.AQLQuery(self.aql_uniparc, bindVars={'l': list(fetch_uniparc)}, rawResults=True)
        for o in query:
            uniparc_id = o['_from']
            expected_len = fetch_uniparc[uniparc_id]['expected_len']
            if uniparc_id in fetch_uniparc:
                h = o['_to'].split('/')[-1]
                sequence = self.protein_store.get_sequence(h)
                if len(sequence) == expected_len:
                    node_sequence = self.add_node(h, self.re_seq_protein_collection, nodes, {'size': len(sequence)})
                    fetch_uniparc[uniparc_id]['hash'] = o['_to']
                    fetch_uniparc[uniparc_id]['node'] = node_sequence
                else:
                    logger.warning(f'expected length mismatch, expected: {expected_len}, actual: {len(sequence)}')
                # node_sequence = self.add_node(h, self.re_seq_protein_collection, nodes, {'size': len(sequence)})
        reload = False
        for uniparc_id in fetch_uniparc:
            data = fetch_uniparc[uniparc_id]
            if 'hash' not in data:
                i = uniparc_id.split('/')[-1]
                print('load uniparc', i)
                nodes, edges = self.etl_transform_uniparc.transform(self.etl_extract_uniparc.extract(i))
                self.etl_load.load({
                    'nodes': nodes,
                    'edges': edges
                })
                reload = True
        if reload:
            query = self.database.AQLQuery(self.aql_uniparc, bindVars={'l': list(fetch_uniparc)}, rawResults=True)
            for o in query:
                uniparc_id = o['_from']
                expected_len = fetch_uniparc[uniparc_id]['expected_len']
                if uniparc_id in fetch_uniparc:
                    h = o['_to'].split('/')[-1]
                    sequence = self.protein_store.get_sequence(h)
                    if len(sequence) == expected_len:
                        node_sequence = self.add_node(h, self.re_seq_protein_collection, nodes, {'size': len(sequence)})
                        fetch_uniparc[uniparc_id]['hash'] = o['_to']
                        fetch_uniparc[uniparc_id]['node'] = node_sequence
                    else:
                        logger.warning(f'expected length mismatch, expected: {expected_len}, actual: {len(sequence)}')

        return fetch_uniparc

    def get_uniprot_hash(self, fetch_uniprot, nodes, fetch_missing=False):
        query = self.database.AQLQuery(self.aql_uniprot, bindVars={'acc_list': list(fetch_uniprot)}, rawResults=True)
        prot_seq = {}
        reload = False
        for k in query:
            acc = k['acc']
            expected_len = fetch_uniprot[acc]['expected_len']
            # print(acc, k)
            sequences = []
            for seq in k['sprot']:
                h = seq.split('/')[-1]
                if h not in prot_seq:
                    sequence = self.protein_store.get_sequence(h)
                    prot_seq[h] = sequence
                else:
                    sequence = prot_seq[h]
                if len(sequence) == expected_len:
                    sequences.append(seq)
                    node_sequence = self.add_node(h, self.re_seq_protein_collection, nodes, {'size': len(sequence)})
                    fetch_uniprot[acc]['node'] = node_sequence
                    fetch_uniprot[acc]['hash'] = seq
                else:
                    logger.warning(f'expected length mismatch, expected: {expected_len}, actual: {len(sequence)}')
            for seq in k['trembl']:
                h = seq.split('/')[-1]
                if h not in prot_seq:
                    sequence = self.protein_store.get_sequence(h)
                    prot_seq[h] = sequence
                else:
                    sequence = prot_seq[h]
                if len(sequence) == fetch_uniprot[acc]['expected_len']:
                    sequences.append(seq)
                    node_sequence = self.add_node(h, self.re_seq_protein_collection, nodes, {'size': len(sequence)})
                    fetch_uniprot[acc]['node'] = node_sequence
                    fetch_uniprot[acc]['hash'] = seq
                else:
                    logger.warning(f'expected length mismatch, expected: {expected_len}, actual: {len(sequence)}')
            if len(sequences) == 1:
                fetch_uniprot[acc]['hash'] = sequences[0]
            else:
                nodes, edges = self.etl_transform_uniprot.transform(
                    self.etl_extract_uniprot.extract(acc.split('/')[-1]))
                self.etl_load.load({
                    'nodes': nodes,
                    'edges': edges
                })
                reload = True
        if reload:
            logger.warning('reload after new data')
            query = self.database.AQLQuery(self.aql_uniprot, bindVars={'acc_list': list(fetch_uniprot)},
                                           rawResults=True)
            for k in query:
                acc = k['acc']
                expected_len = fetch_uniprot[acc]['expected_len']

                sequences = []
                for seq in k['sprot']:
                    h = seq.split('/')[-1]
                    if h not in prot_seq:
                        sequence = self.protein_store.get_sequence(h)
                        prot_seq[h] = sequence
                    else:
                        sequence = prot_seq[h]
                    if len(sequence) == fetch_uniprot[acc]['expected_len']:
                        sequences.append(seq)
                        node_sequence = self.add_node(h, self.re_seq_protein_collection, nodes, {'size': len(sequence)})
                        fetch_uniprot[acc]['node'] = node_sequence
                        fetch_uniprot[acc]['hash'] = seq
                    else:
                        logger.warning(f'expected length mismatch, expected: {expected_len}, actual: {len(sequence)}')
                for seq in k['trembl']:
                    h = seq.split('/')[-1]
                    if h not in prot_seq:
                        sequence = self.protein_store.get_sequence(h)
                        prot_seq[h] = sequence
                    else:
                        sequence = prot_seq[h]
                    if len(sequence) == fetch_uniprot[acc]['expected_len']:
                        sequences.append(seq)
                        node_sequence = self.add_node(h, self.re_seq_protein_collection, nodes, {'size': len(sequence)})
                        fetch_uniprot[acc]['node'] = node_sequence
                        fetch_uniprot[acc]['hash'] = seq
                    else:
                        logger.warning(f'expected length mismatch, expected: {expected_len}, actual: {len(sequence)}')
                if len(sequences) == 1:
                    fetch_uniprot[acc]['hash'] = sequences[0]
        return fetch_uniprot

    def add_node(self, node_id, label, nodes, data=None):
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

        logger.error('add_node error')
        return None

    def add_edge(self, node_from, node_to, label, edges, data=None):
        if label not in edges:
            edges[label] = []
        if label in edges:
            _edge = self.transform_edge(node_from, node_to, data)
            edges[label].append(_edge)
            return _edge

        logger.error('add_edge error')
        return None

    def transform(self, r):
        nodes = {}
        edges = {}

        node_uniref = self.add_node(r['id'], self.uniref_collection, nodes, {
            'member count': r.get('member count', []),
            'common taxon': r.get('common taxon', []),
            'common taxon ID': r.get('common taxon ID', []),
            'name': r.get('name', []),
            'GO Molecular Function': r.get('GO Molecular Function', []),
            'GO Cellular Component': r.get('GO Cellular Component', [])
        })

        fetch_uniprot = {}
        fetch_uniparc = {}
        for o in r['members']:
            for db_ref in o['db_reference']:
                member_type = 'member'
                if 'isSeed' in db_ref and db_ref['isSeed'][0] == 'true':
                    member_type = 'seed'
                if db_ref['type'] == 'UniProtKB ID':
                    fetch_uniprot['uniprotkb_accession/' + db_ref['UniProtKB accession'][0]] = {
                        'expected_len': int(db_ref['length'][0]),
                        'member_type': member_type,
                        'db_reference': db_ref
                    }

                elif db_ref['type'] == 'UniParc ID':
                    fetch_uniparc['uniparc/' + db_ref['id']] = {
                        'expected_len': int(db_ref['length'][0]),
                        'member_type': member_type,
                        'db_reference': db_ref
                    }
                else:
                    logger.warning('not sure what to do:' + db_ref['type'])

        fetch_uniparc = self.get_uniparc_hash(fetch_uniparc, nodes)
        fetch_uniprot = self.get_uniprot_hash(fetch_uniprot, nodes)

        node_links = {}

        for representative_member in r['representative_members']:
            sequence = representative_member['sequence']['value']
            protein_hash = self.protein_store.store_sequence(sequence)
            node_sequence = self.add_node(protein_hash, self.re_seq_protein_collection, nodes, {'size': len(sequence)})
            if node_sequence.id not in node_links:
                node_links[node_sequence.id] = {'type': set(), 'reference': [], 'node': node_sequence}
            node_links[node_sequence.id]['type'].add('representative')
            for ref in representative_member['db_reference']:
                node_links[node_sequence.id]['reference'].append(ref)
                if 'isSeed' in ref and ref['isSeed'][0] == 'true':
                    node_links[node_sequence.id]['type'].add('seed')
            # add_edge(node_uniprotkb, node_sequence, self.uniprot_collection_has_protein_sequence)

        for uniprotkb in fetch_uniprot:
            o = fetch_uniprot[uniprotkb]
            if 'hash' in o:
                if o['hash'] not in node_links:
                    node_links[o['hash']] = {'type': set(), 'reference': []}
                node_links[o['hash']]['type'].add(o['member_type'])
                node_links[o['hash']]['reference'].append(o)
                node_links[o['hash']]['node'] = o['node']
            else:
                raise Exception('unable to match sequences for:', r['id'])
        for uniparc in fetch_uniparc:
            o = fetch_uniparc[uniparc]
            if 'hash' in o:
                if o['hash'] not in node_links:
                    node_links[o['hash']] = {'type': set(), 'reference': []}
                node_links[o['hash']]['type'].add(o['member_type'])
                node_links[o['hash']]['reference'].append(o)
                node_links[o['hash']]['node'] = o['node']
            else:
                raise Exception('unable to match sequences for:', r['id'])

        for link in node_links:
            data = node_links[link]
            self.add_edge(node_uniref, data['node'], self.uniref_has_member, edges, {
                'type': list(sorted(data['type'])),
                'reference': data['reference']
            })

        return nodes, edges
