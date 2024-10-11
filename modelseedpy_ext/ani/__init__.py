import os
from pathlib import Path
from abc import ABC, abstractmethod
from modelseedpy_ext.ani.skani import ProgramSkani
from modelseedpy_ext.ani.fastani import ProgramFastANI


class ANIMethod(ABC):

    @abstractmethod
    def build_library(self, genomes, output_filename):
        pass

    @abstractmethod
    def distance(self, library_query, library_reference, output_filename):
        pass


class GenomeReference:

    def __init__(self):
        self.path = None

    @staticmethod
    def from_file(filename):
        ref = GenomeReference()
        ref.path = filename
        return ref

    @staticmethod
    def from_kbase(info):
        pass

    def get_path(self):
        return Path(self.path)

    def filename(self):
        return self.get_path().name

    def get_absolute_path(self):
        return str(self.get_path().absolute())


class ANIMethodSkani(ANIMethod):

    def __init__(self, program: ProgramSkani, threads=1):
        self.program = program
        self.threads = threads

    def build_library(self, genome_pointers, output_filename):
        if not output_filename.endswith('.txt'):
            raise ValueError('invalid library filename: must end with .txt')

        sketch_path = output_filename[:-4]
        if os.path.exists(sketch_path):
            raise ValueError(f'invalid library filename: {sketch_path} exists')

        with open(output_filename, 'w') as fh:
            for g in genome_pointers:
                fh.write(g.get_absolute_path() + '\n')

        self.program.sketch(sketch_path, output_filename, "--slow", threads=self.threads)

        sketch_library = sketch_path + '.sketch.txt'
        with open(sketch_library, 'w') as fh:
            for filename in os.listdir(sketch_path):
                if filename.endswith('.sketch'):
                    fh.write(sketch_path + '/' + filename + '\n')

        return sketch_library

    def distance(self, library_query, library_reference, output_filename):
        return self.program.run(library_query, library_reference, output_filename, threads=self.threads)


class ANIMethodFastANI(ANIMethod):

    def __init__(self, program: ProgramFastANI, threads=1):
        self.program = program
        self.threads = threads

    def build_library(self, genome_pointers, output_filename):
        if not output_filename.endswith('.txt'):
            raise ValueError('invalid library filename: must end with .txt')

        with open(output_filename, 'w') as fh:
            for g in genome_pointers:
                fh.write(g.get_absolute_path() + '\n')

        return output_filename

    def distance(self, library_query, library_reference, output_filename):
        return self.program.run(library_query,
                                library_file=library_reference,
                                output_file=output_filename,
                                threads=self.threads)
