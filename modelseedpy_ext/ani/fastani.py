import subprocess
import pandas as pd
import os


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


def get_contig_from_genome(genome):
    assembly = kbase.get_from_ws(genome.assembly_ref)
    return get_contig_from_assembly(assembly)


def get_contig_from_assembly(assembly, token):
    file_id = assembly.fasta_handle_ref
    file_path = "/home/fliu/KE/data/"
    res = kbase.download_file_from_kbase(token, file_id, file_path)
    return res
