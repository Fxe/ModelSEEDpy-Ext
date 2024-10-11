import subprocess
import os
import polars as pl


def _read_search_output_as_parquet(filename, sep='\t'):
    import pyarrow as pa

    float_types = {'ANI', 'Align_fraction_ref', 'Align_fraction_query'}

    with open(filename, 'r') as fh:
        header = fh.readline()

        names = header.strip().split(sep)
        col_number = len(names)
        data = {n: [] for n in names}

        line = fh.readline()
        while line:
            row = line.strip().split(sep)
            for i in range(col_number):
                _d = row[i]
                if names[i] in float_types:
                    _d = float(_d)
                data[names[i]].append(_d)
            line = fh.readline()

        pa_table = pa.Table.from_arrays([pa.array(data[n]) for n in names], names=names)
        return pa_table


def _read_triangle_output_as_parquet(filename, sep='\t'):
    import pyarrow as pa

    data = {
        'Query_file1': [],
        'Query_file2': [],
        'ANI': [],
    }
    index = {}
    counter = 0
    with open(filename, 'r') as fh:
        line = fh.readline()
        sz = int(line.strip())
        print(sz)
        line = fh.readline()
        while line:
            _parts = line.strip().split(sep)
            h = _parts[0]
            index[counter] = h
            _scores = _parts[1:]
            for i in range(len(_scores)):
                _val = float(_scores[i])
                if _val != 0:
                    data['Query_file1'].append(h)
                    data['Query_file2'].append(index[i])
                    data['ANI'].append(_val)
            counter += 1
            line = fh.readline()

        names = ['Query_file1', 'Query_file2', 'ANI']
        pa_table = pa.Table.from_arrays([pa.array(data[n]) for n in names], names=names)
        return pa_table


class ANISkaniMatrix:

    def __init__(self, df_data):
        self.df_data = df_data
        #df_scores_alexey.filter(pl.col("ANI") >= 95).group_by("Query_file", maintain_order=True).max()
        pass

    def filter_ani(self, ani):
        return self.df_data.filter(pl.col("ANI") >= ani)

    @staticmethod
    def from_search_output(filename, sep='\t'):
        return ANISkaniMatrix(pl.from_arrow(_read_search_output_as_parquet(filename, sep=sep)))


class ProgramSkani:

    def __init__(self, path_to_bin="skani"):
        self.bin = path_to_bin

    def help(self):
        cmd = [self.bin, "--help"]
        output = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )
        return output

    def version(self):
        cmd = [self.bin, "--version"]
        output = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )
        return output

    def sketch(self, output_folder, library_file, compression_factor="--slow", threads=1):
        cmd = [
            self.bin, "sketch",
            "-t", str(threads),
            compression_factor,
            "-l",
            library_file,
            "-o",
            output_folder
        ]
        output = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )

        return output

    @staticmethod
    def write_library(library_folder, output):
        with open(output, 'w') as fh:
            for f in os.listdir(library_folder):
                if f.endswith('.sketch'):
                    fh.write(f'{library_folder}/{f}\n')

    def triangle_old(self, output_file, library_file, c=125):
        cmd = [
            self.bin,
            "-o",
            output_file,
            "-l",
            library_file,
            "-c",
            c
        ]
        output = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )
        file_size = os.stat(output_file)
        if file_size.st_size > 0:
            pass

        return output

    @staticmethod
    def read_lower_trig_mat(filename):
        index = {}
        counter = 0
        ani_scores_sk = {}
        with open(filename, 'r') as fh:
            line = fh.readline()
            sz = int(line.strip())
            print(sz)
            line = fh.readline()
            while line:
                _parts = line.strip().split('\t')
                h = _parts[0].split('/')[-1][:-4]
                index[counter] = h
                _scores = _parts[1:]
                for i in range(len(_scores)):
                    _val = float(_scores[i])
                    if _val != 0:
                        ani_scores_sk[(h, index[i])] = (_val, None, None)
                counter += 1
                line = fh.readline()
        return ani_scores_sk

    def run(self,
            query_library,
            library_folder="/scratch/fliu/data/ani/skani/slow/gtdb_r220_rep",
            output_file="./skani.out",
            threads=1):

        cmd = [
            self.bin,
            'search',
            "-t", str(threads),
            "--ql", query_library,
            "-d", library_folder,
            "-o", output_file,
        ]

        output = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )
        ani_matrix = None

        if os.path.exists(output_file) and os.stat(output_file).st_size > 0:
            ani_matrix = _read_search_output_as_parquet(output_file)
        return ani_matrix, output, cmd

    def triangle(self,
                 query_library,
                 compression_factor,
                 output_file="./skani.out",
                 threads=1):

        cmd = [
            self.bin,
            'triangle',
            "-t", str(threads),
            compression_factor,
            "-l", query_library,
            "-o", output_file,
        ]

        output = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )
        ani_matrix = None

        if os.path.exists(output_file) and os.stat(output_file).st_size > 0:
            ani_matrix = _read_triangle_output_as_parquet(output_file)
        return ani_matrix, output, cmd


class ANIMethodSkani:

    def __init__(self, program: ProgramSkani):
        self.program = program

    def build_library(self, genome_pointers, output_filename):
        if not output_filename.endswith('.txt'):
            raise ValueError('library filename must end with .txt')

        self.program.sketch
        # write
        pass

    def distance(self, library_query, library_reference, output_filename):

        pass