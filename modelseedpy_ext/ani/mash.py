import subprocess

class ProgramMash:
    
    def __init__(self, path_to_bin='mash'):
        self.bin = path_to_bin
        
    def sketch(self, library_file, output_mash_library, sketch_size=1000, threads=1):
        cmd = [self.bin, 'sketch',
        '-p', str(threads),
        '-s', str(sketch_size),
        '-l', library_file, 
        '-o', output_mash_library]
        output = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        return output
    
    def dist(self, library_file, target_genome_file, threads=1):
        cmd = [self.bin, 'dist', library_file, target_genome_file]
        output = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        return output
