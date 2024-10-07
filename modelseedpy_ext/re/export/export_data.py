import pyarrow as pa
import pyarrow.parquet as pq


def export_ncbi_assembly(ncbi_assembly_data, filename):
    attrs_summary = {}

    for o in ncbi_assembly_data:
        for k in o:
            v = o[k]
            v_type = type(v)
            if k not in attrs_summary:
                attrs_summary[k] = {}
            if v_type not in attrs_summary[k]:
                attrs_summary[k][v_type] = []
            if len(attrs_summary[k][v_type]) < 4:
                attrs_summary[k][v_type].append(v)

    big_table = {
        'PropertyList': [],
        'Synonym_Genbank': [],
        'Synonym_RefSeq': [],
        'Synonym_Similarity': [],
    }
    ignore = {'_rev', 'Synonym', 'PropertyList'}
    for k in attrs_summary:
        if k not in ignore:
            big_table[k] = []
    for o in ncbi_assembly_data:
        for k in big_table:
            if k == 'PropertyList':
                big_table[k].append(';'.join(o.get(k, [])))
            elif k == 'Coverage':
                v = o.get('Coverage')
                if v:
                    v = float(v)
                big_table[k].append(v)
            elif k == 'ncbi_id':
                v = o.get('ncbi_id')
                if v:
                    v = str(v)
                big_table[k].append(v)
            elif k == 'Synonym_Genbank':
                big_table[k].append(o.get('Synonym', {}).get('Genbank'))
            elif k == 'Synonym_RefSeq':
                big_table[k].append(o.get('Synonym', {}).get('RefSeq'))
            elif k == 'Synonym_Similarity':
                big_table[k].append(o.get('Synonym', {}).get('Similarity'))
            else:
                big_table[k].append(o.get(k, None))

    pa_arrays = []
    names = []
    for k in big_table:
        names.append(k)
        pa_arrays.append(pa.array(big_table[k]))

    pa_table = pa.Table.from_arrays(pa_arrays, names=names)

    pq.write_table(pa_table, filename)
