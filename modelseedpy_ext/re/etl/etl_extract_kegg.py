import os


class ETLExtractKEGG:

    def __init__(self):
        pass

    def load_genomes(self, base_dir: str):
        genomes = {}
        for f in os.listdir(base_dir):
            if f[0] == 'T':
                filename = f'{base_dir}/{f}'
                kegg_genome = transform(read_kegg_flatfile(filename))
                if kegg_genome['_key'] not in genomes:
                    genomes[kegg_genome['_key']] = kegg_genome
        return genomes


def read_kegg_flatfile(filename):
    cursor = None
    with open(filename, 'r', encoding="utf8") as fh:
        data = {}
        l = fh.readline()
        while l:
            #print(l[:-1])
            current = l[:12].strip()
            #print(current)
            rest = l[12:]
            if current and current != cursor:
                cursor = current
                if cursor not in data:
                    data[cursor] = {'value': [], 'subfields': {}}
                data[cursor]['value'].append(rest)
            else:
                data[cursor]['value'].append(rest)
                #print(cursor)
            #print(l[:10], l[10:])
            if current == '///':
                return data
            l = fh.readline()


def transform(res):
    a, b = res['ENTRY']['value'][0].strip().split(' ', 1)
    a, b.strip()
    kegg_genome = {
        '_key': a,
        'name': res['NAME']['value'][0].strip(),
        'taxonomy': res['TAXONOMY']['value'][0].strip(),
        'lineage': res['LINEAGE']['value'][0].strip(),
    }
    if 'ORG_CODE' in res:
        kegg_genome['code'] = res['ORG_CODE']['value'][0].strip()
    if 'DATA_SOURCE' in res:
        kegg_genome['data_source'] = [x.strip() for x in res['DATA_SOURCE']['value']]
    if 'STATISTICS' in res:
        kegg_genome['statistics'] = [x.strip() for x in res['STATISTICS']['value']]
    if 'CHROMOSOME' in res:
        kegg_genome['chromosome'] = res['CHROMOSOME']['value'][0].strip()
    if 'COMMENT' in res:
        kegg_genome['comments'] = [x.strip() for x in res['COMMENT']['value']]
    if 'CATEGORY' in res:
        kegg_genome['category'] = res['CATEGORY']['value'][0].strip()
    if 'DISEASE' in res:
        kegg_genome['disease'] = res['DISEASE']['value'][0].strip()
    return kegg_genome