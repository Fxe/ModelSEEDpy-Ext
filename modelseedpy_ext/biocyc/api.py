import os
import requests
from io import BytesIO
import lxml.etree as et
from enum import Enum


class BiocycQuery(Enum):
    Pathways = 'Pathway'
    Reactions = 'Reaction'
    Compounds = 'Compound'
    Proteins = 'Protein'
    Genes = 'Gene'
    Organisms = 'Organism'


def _detect_class(data):
    parser = et.iterparse(BytesIO(data), events=("end", "start"))
    for action, elem in parser:
        if action == "end" and elem.tag == "metadata":
            for _action, _elem in parser:
                if _action == 'start':
                    return _elem.tag

    raise Exception("unable to detect class")


def _to_list(data):
    res = {}

    parser = et.iterparse(BytesIO(data), events=("end", "start"))
    for _action, _elem in parser:
        if _action == "end" and _elem.tag == "metadata":
            for action, elem in parser:
                if action == 'start':
                    if elem.tag not in res:
                        res[elem.tag] = []
                    res[elem.tag].append([dict(elem.attrib), elem.text])
                elif action == "end" and elem.tag == "ptools-xml":
                    return res

    raise Exception("unexpected ending")


class BiocycAPI:

    def __init__(self, db='META', email=None, password=None):
        self.db = db
        self.end_point = 'https://websvc.biocyc.org'
        self.session = None
        if email and password:
            self.session = requests.Session()
            self.session.post(f'{self.end_point}/credentials/login/',
                              data={'email': email, 'password': password})

    def get_session(self):
        if self.session:
            return self.session
        else:
            return requests

    def fetch(self, frame_id, detail='high'):
        resp = self.get_session().get(f'{self.end_point}/getxml?id={self.db}:{frame_id}&detail={detail}')
        return resp.content

    def query(self, biocyc_query: BiocycQuery, detail='none'):
        """
        :param biocyc_query: query entity type
        :param detail: [none|low|full]
        """
        q = f'[x:x<-{self.db}^^{biocyc_query.name}]'

        resp = self.get_session().get(f'{self.end_point}/xmlquery?query={q}&detail={detail}')

        items = _to_list(resp.content)

        return {o[0]['frameid'] for o in items[biocyc_query.value]}


class BiocycAPICached(BiocycAPI):

    def __init__(self, cache, db='META', email=None, password=None):
        self.cache = cache
        super().__init__(db, email, password)

    @staticmethod
    def _to_frame_id_file(frame_id):
        return frame_id.replace('+', '__PLUS__')

    @staticmethod
    def _from_frame_id_file(frame_id):
        return frame_id.replace('__PLUS__', '+')

    def cache_fetch(self, frame_id):
        if not os.path.exists(f'{self.cache}/{self.db}'):
            return None

        frame_file = self._to_frame_id_file(frame_id) + '.xml'
        for d in os.listdir(f'{self.cache}/{self.db}'):
            path = f'{self.cache}/{self.db}/{d}'
            if os.path.isdir(path):
                files = os.listdir(path)
                if frame_file in files:
                    with open(f'{self.cache}/{self.db}/{d}/{frame_file}', 'rb') as fh:
                        return fh.read()
        return None

    def cache_store(self, frame_id, data, object_type):
        if not os.path.exists(f'{self.cache}/{self.db}'):
            os.mkdir(f'{self.cache}/{self.db}')
        if not os.path.exists(f'{self.cache}/{self.db}/{object_type}'):
            os.mkdir(f'{self.cache}/{self.db}/{object_type}')
        frame_file = self._to_frame_id_file(frame_id) + '.xml'
        with open(f'{self.cache}/{self.db}/{object_type}/{frame_file}', 'wb') as fh:
            fh.write(data)

    def cache_query(self, biocyc_query: BiocycQuery):
        if not os.path.exists(f'{self.cache}/{self.db}'):
            return None
        if not os.path.exists(f'{self.cache}/{self.db}/query'):
            return None
        if not os.path.exists(f'{self.cache}/{self.db}/query/{biocyc_query.name}.txt'):
            return None
        with open(f'{self.cache}/{self.db}/query/{biocyc_query.name}.txt', 'r') as fh:
            return fh.read().split()

    def save_query(self, biocyc_query: BiocycQuery, items):
        if not os.path.exists(f'{self.cache}/{self.db}'):
            os.mkdir(f'{self.cache}/{self.db}')
        if not os.path.exists(f'{self.cache}/{self.db}/query'):
            os.mkdir(f'{self.cache}/{self.db}/query')
        with open(f'{self.cache}/{self.db}/query/{biocyc_query.name}.txt', 'w') as fh:
            for i in items:
                fh.write(i + '\n')

    def fetch(self, frame_id, detail='high'):
        data = self.cache_fetch(frame_id)
        if data:
            return data
        else:
            data = super().fetch(frame_id, detail)
            object_type = _detect_class(data)
            print(object_type)
            self.cache_store(frame_id, data, object_type)

            return data

    def query(self, biocyc_query: BiocycQuery, detail='none'):
        items = self.cache_query(biocyc_query)
        if items is None:
            print('loading!!')
            items = super().query(biocyc_query, detail)
            self.save_query(biocyc_query, items)

        return items
