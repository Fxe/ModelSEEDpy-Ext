import subprocess
import pandas as pd
import os


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
):
    cmd = [
        "/home/fliu/c++/FastANI/fastANI",
        "-t",
        str(threads),
        "-q",
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
