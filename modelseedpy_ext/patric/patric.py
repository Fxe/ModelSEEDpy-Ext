import os
from ftplib import FTP


class PatricBulkDownloader:

    def __init__(self, base_path, max_block_size=5000):
        self.block_mapping = {0: set()}
        self.errors = {}
        self.current_block = 0
        self.max_block_size = max_block_size
        self.base_path = base_path
        self._PATRIC_FTP_URL = 'ftp.bvbrc.org'

    def get_next_block(self):
        next_block = 0
        for k in self.block_mapping:
            if len(self.block_mapping[k]) < self.max_block_size:
                next_block = int(k)
            else:
                next_block = int(k) + 1
        return next_block

    def _download_content(self, genome_summary, file_func):
        ftp = FTP(self._PATRIC_FTP_URL)
        ftp.login()
        num_records = len(genome_summary) - 1
        for i in range(num_records):
            m = genome_summary[i + 1].split('\t')[0]
            try:
                if len(self.block_mapping[self.current_block]) >= self.max_block_size:
                    self.current_block = self.get_next_block()
                if self.current_block not in self.block_mapping:
                    self.block_mapping[self.current_block] = set()
                _folder = f'{self.base_path}/{self.current_block}'
                f_name = file_func(m)
                if not os.path.exists(_folder):
                    os.mkdir(_folder)
                with open(f'{_folder}/{f_name}', 'wb') as fh:
                    ftp.retrbinary(f"RETR /genomes/{m}/{f_name}", fh.write, 8 * 1024)
                    self.block_mapping[self.current_block].add(m)
            except Exception as ex:
                self.errors[m] = str(ex)
        ftp.close()

    def download_genomes(self, genome_summary):
        def faa_file(i):
            return f'{i}.PATRIC.faa'
        self._download_content(genome_summary, faa_file)

    def download_contigs(self, genome_summary):
        def fna_file(i):
            return f'{i}.fna'
        self._download_content(genome_summary, fna_file)
