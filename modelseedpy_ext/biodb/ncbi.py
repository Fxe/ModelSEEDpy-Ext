import os
from modelseedpy_ext.re.core.genome import REAssembly
from modelseedpy.core import MSGenome
from modelseedpy_ext.re.core.genome import GffRecord


FTP_URL_PREFIX = "ftp://ftp.ncbi.nlm.nih.gov"


class NCBIAssembly:

    def __init__(self, data, cache=None):
        self.data = {}
        self.ftp_path_rs = data['FtpPath_RefSeq']
        self.ftp_path_gb = data['FtpPath_GenBank']
        self.cache_folder = cache

    @property
    def cwd_ftp_path_rs(self):
        if self.ftp_path_rs and self.ftp_path_rs.startswith(FTP_URL_PREFIX):
            _url_p = self.ftp_path_rs.split('/')
            return self.ftp_path_rs.split(FTP_URL_PREFIX)[1]

        return None

    @property
    def cwd_local_path_rs(self):
        if self.ftp_path_rs and self.ftp_path_rs.startswith(FTP_URL_PREFIX):
            _url_p = self.ftp_path_rs.split('/')
            return f"{self.cache_folder}/{'/'.join(_url_p[3:])}"

        return None

    @property
    def cwd_ftp_path_gb(self):
        if self.ftp_path_gb and self.ftp_path_gb.startswith(FTP_URL_PREFIX):
            _url_p = self.ftp_path_gb.split('/')
            return self.ftp_path_gb.split(FTP_URL_PREFIX)[1]

        return None

    @property
    def cwd_local_path_gb(self):
        if self.ftp_path_gb and self.ftp_path_gb.startswith(FTP_URL_PREFIX):
            _url_p = self.ftp_path_gb.split('/')
            return f"{self.cache_folder}/{'/'.join(_url_p[3:])}"

        return None

    def fetch_ncbi_ftp_data(self, ftp):
        """

        :param ftp: FTP client
        :return:
        """
        ftp_path = self.cwd_ftp_path_rs
        if ftp_path is None:
            ftp_path = self.cwd_ftp_path_gb

        if ftp_path:
            write_path = f'{self.cache_folder}/{ftp_path}'
            os.makedirs(write_path, exist_ok=True)

            ftp.cwd(ftp_path)
            files = ftp.nlst()
            for f in files:
                target_file = f'{write_path}/{f}'
                if f.endswith('_assembly_structure'):  # TODO: implement fetch _assembly_structure
                    os.makedirs(target_file, exist_ok=True)
                else:
                    with open(f'{write_path}/{f}', 'wb') as fh:
                        ftp.retrbinary(f"RETR {f}", fh.write)

    @property
    def local_path(self):
        local = self.cwd_local_path_rs
        if local is None:
            local = self.cwd_local_path_gb
        return local

    @property
    def local_genomic_fna_path(self):
        local = self.local_path
        if local is None:
            return None
        _p = local.split('/')
        file_fna = f'{_p[-1]}_genomic.fna.gz'
        return f'{local}/{file_fna}'

    def get_genomic_fna(self):
        local_genomic_fna_path = self.local_genomic_fna_path
        if os.path.exists(local_genomic_fna_path):
            return REAssembly.from_fasta(local_genomic_fna_path)

        raise ValueError('cache not found: ' + local_genomic_fna_path)

    def get_protein_faa(self):
        local = self.local_path
        if local and os.path.exists(local):
            _p = local.split('/')
            file_faa = f'{_p[-1]}_protein.faa.gz'
            if os.path.exists(f'{local}/{file_faa}'):
                return MSGenome.from_fasta(f'{local}/{file_faa}')

        raise ValueError('cache not found: ' + local)

    def get_gff(self):
        local = self.local_path
        if local and os.path.exists(local):
            _p = local.split('/')
            file_gff = f'{_p[-1]}_genomic.gff.gz'

            if os.path.exists(f'{local}/{file_gff}'):
                return _read_gff_features(f'{local}/{file_gff}')

        raise ValueError('cache not found: ' + local)
