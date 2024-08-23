import subprocess
import os
import pyarrow as pa


class ProgramSkani:

    def __init__(self, path_to_bin="skani"):
        self.bin = path_to_bin

    def sketch(self, output_folder, library_file):
        cmd = [
            self.bin,
            "-o",
            output_folder,
            "-l",
            library_file
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

    def triangle(self, output_file, library_file, c=125):
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
