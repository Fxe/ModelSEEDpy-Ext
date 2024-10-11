import subprocess


class ProgramSkani:

    def __init__(self, path_to_bin="sourmash"):
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