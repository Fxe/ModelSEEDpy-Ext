import subprocess
import os


class GeneMark2S:

    def __init__(self):
        self.bin = '/opt/build/genemark-2/gms2_linux_64/gms2.pl'

    def run(self, input_fna, output_faa, output_fna, work_dir='.', genome_type="bacteria"):
        if not os.path.exists(input_fna):
            raise FileNotFoundError('input_fna')

        if not os.path.exists(work_dir):
            os.makedirs(work_dir, exist_ok=True)

        cmd_args = [
            self.bin,
            '--seq', input_fna,
            '--genome-type', genome_type,
            '--faa', output_faa,
        ]

        if output_fna:
            cmd_args.extend(['--fnn', output_fna])

        _env = os.environ.copy()

        try:
            # Run subprocess
            result = subprocess.run(
                cmd_args,
                env=_env,
                cwd=work_dir,
                capture_output=True,
                text=True,
                check=True
            )

            return {
                'success': True,
                'output_files': os.listdir(work_dir),
                'command': ' '.join(cmd_args),
                'stdout': result.stdout,
                'stderr': result.stderr,
                'returncode': result.returncode
            }

        except subprocess.CalledProcessError as e:
            raise Exception(f"failed with return code {e.returncode}: {e.stderr}")
        except FileNotFoundError:
            raise Exception(f"executable not found at {self.bin}")
