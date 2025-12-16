import json
import requests
import pandas as pd


ANL_ARGO_ENDPOINT_OLD = 'https://apps-test.inside.anl.gov/argoapi/api/v1/resource/chat/'
ANL_ARGO_ENDPOINT = 'https://apps.inside.anl.gov/argoapi/api/v1/resource/chat/'
ANL_ARGO_ENDPOINT_DEV = 'https://apps-dev.inside.anl.gov/argoapi/api/v1/resource/chat/'


class Argo:

    def __init__(self, engine, endpoint=ANL_ARGO_ENDPOINT_DEV):
        self.url = endpoint
        self.engine = engine
        self.user = 'fliu_bot'

    def gpt_query(self, q, model):
        _body = {
            "user": self.user,
            "system": "Answer the below question.",
            "model": model,
            "prompt": [q],
            "stop": []
        }
        response = requests.post(self.url, json=_body)
        d_response = json.loads(response.content)['response']
        return d_response

    def query(self, q, model='gpt4o'):
        prompt = f"""
        {q}

        Return the answer in the following json format:

        KEGG_ORTHOLOG_ID: explanations why it was chosen

        Example:

        K00001: reason because it is relevant

        If the question is not related to KEGG orthologs answer: NA
        """
        _body = {
            "user": self.user,
            "system": "Answer the below question in python json text.",
            "model": model,
            "prompt": [prompt],
            "stop": []
        }
        response = requests.post(self.url, json=_body)
        p_response, err = self.process_answer(response)

        if type(p_response) == str:
            return {
                'success': False,
                'response': p_response,
                'error': str(err)
            }
        else:
            if 'NA' in p_response:
                return {
                    'success': False,
                    'response': p_response,
                }
            return {
                'success': True,
                'response': self.process_ko_query(p_response),
            }

    def process_ko_query(self, p_response):
        response = {}
        for ko in p_response:
            query = f"SELECT * FROM `feature_x_kegg_ortholog` WHERE id_ko = 'KO:{ko}';"
            df = pd.read_sql(query, self.engine)
            response[ko] = {
                'query': query,
                'explain': p_response[ko],
                'data': df
            }
        return response

    @staticmethod
    def process_answer(response):
        # b_content = response.content
        d_response = json.loads(response.content)['response']
        clean_json_data = d_response.replace("```json", "").replace("```", "").strip()
        try:
            _json_data = json.loads(clean_json_data)
            return _json_data, None
        except Exception as ex:
            return d_response, ex
