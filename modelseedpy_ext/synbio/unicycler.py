import subprocess
import os


class Unicycler:

    def __init__(self):
        self.bin = '/opt/env/modelseed/bin/unicycler'
        self.spades_bin = '/opt/build/spades/bin'
        self.racon_bin = '/opt/build/racon.build/bin'

    def run(self, reads1, reads2, long_reads, output, threads=1):
        #if not os.path.exists(reads1):
        #    raise FileNotFoundError('reads1')
        #if not long and not os.path.exists(reads2):
        #    raise FileNotFoundError('reads2')

        cmd_args = [self.bin]
        if reads1 and reads2:
            cmd_args.extend([
                '-1', reads1,
                '-2', reads2,
            ])
        if long_reads:
            cmd_args.extend([
                '-l', long_reads
            ])

        cmd_args.extend([
            '-t', str(threads),
            '-o', output
        ])

        _env = os.environ.copy()
        _env["PATH"] = self.spades_bin + ":" + _env["PATH"]
        _env["PATH"] = self.racon_bin + ":" + _env["PATH"]

        try:
            # Run subprocess
            result = subprocess.run(
                cmd_args,
                env=_env,
                capture_output=True,
                text=True,
                check=True
            )

            return {
                'success': True,
                'output_files': os.listdir(output),
                'command': ' '.join(cmd_args),
                'stdout': result.stdout,
                'stderr': result.stderr,
                'returncode': result.returncode
            }

        except subprocess.CalledProcessError as e:
            raise Exception(f"failed with return code {e.returncode}: {e.stderr}")
        except FileNotFoundError:
            raise Exception(f"executable not found at {self.bin}")
