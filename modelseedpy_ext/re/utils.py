import logging


logger = logging.getLogger(__name__)

def generate_ec_hierarchy(ec_lv1, ec_lv2, ec_lv3, ec_lv4):
    res = {}
    for ec in ec_lv1:
        res[ec] = {}
    for ec in ec_lv2:
        a, _, _, _ = ec.split('.')
        parent = f'{a}.-.-.-'
        if parent not in res:
            print(f'generate parent ec: {parent}')
            res[parent] = {}
        res[parent][ec] = {}
    for ec in ec_lv3:
        a, b, _, _ = ec.split('.')
        parent_lv1 = f'{a}.-.-.-'
        parent_lv2 = f'{a}.{b}.-.-'
        if parent_lv1 not in res:
            print(f'generate parent ec: {parent_lv1}')
            res[parent_lv1] = {}
        if parent_lv2 not in res[parent_lv1]:
            print(f'generate parent ec: {parent_lv2}')
            res[parent_lv1][parent_lv2] = {}
        res[parent_lv1][parent_lv2][ec] = set()
    for ec in ec_lv4:
        a, b, c, _ = ec.split('.')
        parent_lv1 = f'{a}.-.-.-'
        parent_lv2 = f'{a}.{b}.-.-'
        parent_lv3 = f'{a}.{b}.{c}.-'
        if parent_lv1 not in res:
            print(f'generate parent ec: {parent_lv1}')
            res[parent_lv1] = {}
        if parent_lv2 not in res[parent_lv1]:
            print(f'generate parent ec: {parent_lv2}')
            res[parent_lv1][parent_lv2] = {}
        if parent_lv3 not in res[parent_lv1][parent_lv2]:
            print(f'generate parent ec: {parent_lv3}')
            res[parent_lv1][parent_lv2][parent_lv3] = set()
        res[parent_lv1][parent_lv2][parent_lv3].add(ec)

    return res

def build_ec_tree(re, load_arango):
    col_ec_number = re.db['EC_number']
    ec_numbers = {doc['_key'] for doc in col_ec_number.fetchAll(rawResults=True)}

    print('total ec values found', len(ec_numbers))

    ec_lv1 = set()
    ec_lv2 = set()
    ec_lv3 = set()
    ec_lv4 = set()
    for ec in ec_numbers:
        if ec.count('.') == 3:
            a, b, c, d = ec.split('.')
            if b == '-' and c == '-' and d == '-':
                ec_lv1.add(ec)
            elif c == '-' and d == '-':
                ec_lv2.add(ec)
            elif d == '-':
                ec_lv3.add(ec)
            else:
                ec_lv4.add(ec)
        else:
            print('bad', ec)

    ec_hierarchy = generate_ec_hierarchy(ec_lv1, ec_lv2, ec_lv3, ec_lv4)

    new_ecs = set()
    for l1 in ec_hierarchy:
        if l1 not in ec_numbers:
            new_ecs.add(l1)
        for l2 in ec_hierarchy[l1]:
            if l2 not in ec_numbers:
                new_ecs.add(l2)
            for l3 in ec_hierarchy[l1][l2]:
                if l3 not in ec_numbers:
                    new_ecs.add(l3)
                for l4 in ec_hierarchy[l1][l2][l3]:
                    if l4 not in ec_numbers:
                        new_ecs.add(l4)

    def build_node(node_id, label, data=None):
        return Node(node_id, label, data)

    def add_node(nodes, node_id, label, data=None):
        if label not in nodes:
            nodes[label] = {}
        if label in nodes:
            _node = build_node(node_id, label, data)
            if _node.id not in nodes[label]:
                nodes[label][_node.id] = _node
            else:
                # print('dup', _node.id)
                pass

            return nodes[label][_node.id]

        logger.error('add_node error')
        return None

    nodes = {}
    for ec in new_ecs:
        add_node(nodes, ec, 'EC_number')

    load_arango.load( {
        'nodes': nodes,
        'edges': []
    })

    def transform_edge(src: Node, dst: Node, data=None):
        _edge = {
            '_key': '{}:{}'.format(src.key, dst.key),
            '_from': src,
            '_to': dst,
        }
        if data:
            _edge.update(data)
        return _edge

    def add_edge(edges, node_from, node_to, label, data=None):
        if label not in edges:
            edges[label] = []
        if label in edges:
            _edge = transform_edge(node_from, node_to, data)
            edges[label].append(_edge)
            return _edge

        logger.error('add_edge error')
        return None

    nodes = {}
    for ec in ec_numbers:
        add_node(nodes, ec, 'EC_number')
    edges = {}
    ec_edges = set()
    for l1 in ec_hierarchy:
        for l2 in ec_hierarchy[l1]:
            ec_edges.add((l1, l2))
            for l3 in ec_hierarchy[l1][l2]:
                ec_edges.add((l2, l3))
                for l4 in ec_hierarchy[l1][l2][l3]:
                    ec_edges.add((l3, l4))
    for src, dst in ec_edges:
        node_src = nodes['EC_number'][f'EC_number/{src}']
        node_dst = nodes['EC_number'][f'EC_number/{dst}']
        add_edge(edges, node_src, node_dst, 'ec_number_child', data=None)

    load_arango.load({
        'nodes': nodes,
        'edges': edges
    })

    print('done!')