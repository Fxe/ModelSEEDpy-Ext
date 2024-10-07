import subprocess
import pandas as pd
import os


class FastANIOutput:

    def __init__(self):
        self.d_ani = {}

    @staticmethod
    def from_tsv(path):
        d_ani = {}
        with open(path, 'r') as fh:
            line = fh.readline()
            while line:
                if line:
                    a, b, score, f1, f2 = line.split('\t')
                    a = a.split('/')[-1]
                    b = b.split('/')[-1]
                    if a not in d_ani:
                        d_ani[a] = {}
                    d_ani[a][b] = (score, f1, f2)
                line = fh.readline()
        res = FastANIOutput()
        res.d_ani = d_ani
        return res

    @staticmethod
    def to_parquet(path, filename, _fn_g_name=None):
        edges = {}
        with open(path, 'r') as fh:
            line = fh.readline()
            while line:
                if line:
                    a, b, score, f1, f2 = line.strip().split('\t')
                    if a != b:
                        _e = (a, b)
                        _e_rev = (b, a)
                        if _e not in edges and _e_rev not in edges:
                            edges[_e] = [float(score), None]
                        elif _e_rev in edges:
                            edges[_e_rev][1] = float(score)
                    # a = a.split('/')[-1]
                    # b = b.split('/')[-1]
                line = fh.readline()

        _data = {
            'h1': [],
            'h2': [],
            'ani': [],
        }
        for p in edges:
            _s1, _s2 = edges[p]
            if _s1 and _s2:
                _sq_err = (_s1 - _s2) ** 2
                if _sq_err > 1:
                    # print(p, _s1, _s2, _sq_err)
                    pass
                else:
                    _data['h1'].append(_fn_g_name(p[0]))
                    _data['h2'].append(_fn_g_name(p[1]))
                    _data['ani'].append(_s1)

        import pyarrow as pa
        ani_pq = pa.Table.from_arrays([
            pa.array(_data['h1']),
            pa.array(_data['h2']),
            pa.array(_data['ani'])], names=['h1', 'h2', 'ani'])
        import pyarrow.parquet as pq
        pq.write_table(ani_pq, filename)


class ProgramANI:

    def __init__(self, path_to_bin="fastANI"):
        self.bin = path_to_bin

    def run(self,
            query_file,
            query_library=None,
            library_file="/home/fliu/KE/data/fastani_library.txt",
            output_file="/home/fliu/KE/data/fastani.out",
            threads=1):

        cmd = [
            self.bin,
            "-t",
            str(threads),
            "--ql",
            query_file,
            "--rl",
            library_file,
            "-o",
            output_file,
        ]
        output = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )
        fast_ani_result = None
        file_size = os.stat(output_file)
        if file_size.st_size > 0:
            fast_ani_result = pd.read_csv(output_file, header=None, sep="\t")
        return fast_ani_result, output


def run_fastani(
    query_file,
    library_file="/home/fliu/KE/data/fastani_library.txt",
    output_file="/home/fliu/KE/data/fastani.out",
    threads=1,
    bin_fastani='/opt/FastANI/fastANI'
):
    cmd = [
        bin_fastani,
        "-t",
        str(threads),
        "--ql",
        query_file,
        "--rl",
        library_file,
        "-o",
        output_file,
    ]
    output = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )
    fast_ani_result = None
    if os.path.exists(output_file):
        file_size = os.stat(output_file)
        if file_size.st_size > 0:
            fast_ani_result = pd.read_csv(output_file, header=None, sep="\t")
        return fast_ani_result, output
    return None, output


def get_contig_from_genome(genome, kbase):
    assembly = kbase.get_from_ws(genome.assembly_ref)
    return get_contig_from_assembly(assembly)


def get_contig_from_assembly(assembly, token, kbase):
    file_id = assembly.fasta_handle_ref
    file_path = "/home/fliu/KE/data/"
    res = kbase.download_file_from_kbase(token, file_id, file_path)
    return res


def build_fast_ani_library(genomes, kbase, output_path='library_ani.txt',
                           kbase_cache='/scratch/fliu/data/kbase/cache/handle'):
    assemblies = {}

    for genome in genomes:
        if genome.id not in assemblies:
            assemblies[genome.id] = kbase.get_from_ws(genome.assembly_ref)
        else:
            print('dup', genome.id)

    with open(output_path, 'w') as fh:
        for genome_id, a in assemblies.items():
            fasta_handle_ref = a.fasta_handle_ref
            filename = f'{kbase_cache}/{fasta_handle_ref}'
            if not os.path.exists(filename):
                kbase.download_file_from_kbase2(fasta_handle_ref, filename)
            else:
                fh.write(filename + '\n')

    return assemblies
