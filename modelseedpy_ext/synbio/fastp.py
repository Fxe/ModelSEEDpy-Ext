import hashlib
import os
import subprocess
from io import BytesIO
from minio import Minio


def run_fastp(in1: str, in2:str, interleaved: bool,
              output_dir: str, threads: int, cmd:str='fastp', env:dict=None):
    """
    Run fastp for quality control and trimming of sequencing reads.
    
    Args:
        reads_dir (str): Directory containing input reads
        interleaved (bool): Whether reads are interleaved
        output_dir (str): Directory for output files
        threads (int): Number of threads to use
        
    Returns:
        dict: Results with output file paths and process info
    """
    _env = os.environ.copy()
    if env:
        _env.update(env)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Build command arguments
    cmd_args = [cmd]
    
    if interleaved:
        cmd_args.extend(['--interleaved_in', '--in1', in1])
    else:
        cmd_args.extend(['--in1', in1, '--in2', in2])
    
    cmd_args.extend([
        '-w', str(threads),
        '--out1', os.path.join(output_dir, '1.trim.fastq.gz'),
        '--out2', os.path.join(output_dir, '2.trim.fastq.gz'),
        '-j', os.path.join(output_dir, 'report.json'),
        '-h', os.path.join(output_dir, 'report.html')
    ])
    
    try:
        # Run fastp subprocess
        result = subprocess.run(
            cmd_args,
            env=_env,
            capture_output=True,
            text=True,
            check=True
        )
        
        return {
            'success': True,
            'output_files': os.listdir(output_dir),
            'command': ' '.join(cmd_args),
            'stdout': result.stdout,
            'stderr': result.stderr,
            'returncode': result.returncode
        }
        
    except subprocess.CalledProcessError as e:
        raise Exception(f"fastp failed with return code {e.returncode}: {e.stderr}")
    except FileNotFoundError:
        raise Exception(f"fastp executable not found at {cmd}")


def to_minio(local_path:str, minio_client: Minio, minio_bucket:str, minio_path:str):
    """
    Upload all files from a local directory to MinIO.
    
    Args:
        local_path (str): Local directory path
        minio_client (Minio): MinIO client instance
        minio_bucket (str): Target bucket name
        minio_path (str): Target path prefix in MinIO
    """
    for f in os.listdir(local_path):
        local_file_path = os.path.join(local_path, f)
        minio_object_name = os.path.join(minio_path, f)
        
        with open(local_file_path, 'rb') as fh:
            data = fh.read()
            minio_client.put_object(
                bucket_name=minio_bucket,
                object_name=minio_object_name,
                data=BytesIO(data),
                length=len(data),
                metadata={
                    'md5': hashlib.md5(data).hexdigest()
                }
            )
