import hashlib
import os
import subprocess


def ETL_reads(reads_dir:str, minio_client, bucket_name:str, object_name:str, type:str):
    # read reads_dir into DataBlob
    # save DataBlob to minio

    
    if type == "illumina": # if illumina, run fastp
        run_fastp(reads_dir, minio_client, bucket_name, object_name, type)
        # read fastp output into DataBlob
        # save fastp output to minio
    elif type == "nanopore": # if nanopore, run Nanoplot
        # save Nanoplot output to minio
        pass
    else:
        raise ValueError(f"Invalid reads type: {type}")
