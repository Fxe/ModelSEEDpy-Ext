import requests
import json


class ETLExtractUniparc:
    def __init__(self, api):
        self.api = api

    def extract(self, entry):
        res = requests.get(f"https://rest.uniprot.org/uniparc/{entry}?format=json")
        if res.status_code == 200:
            return json.loads(res.content)


class ETLExtractUniprot:
    def __init__(self, api):
        self.api = api

    def extract(self, entry):
        res = requests.get(f"https://rest.uniprot.org/uniprotkb/{entry}?format=json")
        if res.status_code == 200:
            return json.loads(res.content)
