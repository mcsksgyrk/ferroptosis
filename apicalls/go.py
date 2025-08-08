from .base import APIClient
from typing import List, Dict, TypedDict


class GOTerm(TypedDict):
    id: str
    name: str
    parents: List[str]
    childrens: List[str]


class GOClient(APIClient):
    """Client for interacting with Reactome API."""

    def __init__(self):
        super().__init__("https://www.ebi.ac.uk/QuickGO/services/ontology/go")

    def get_parent_terms(self, go_id: str) -> List[str]:
        headers = {'Accept': 'application/json'}
        response = self._make_request("GET",
                                      f"terms/{go_id}/ancestors",
                                      headers=headers)

        if response.status_code == 200:
            data = response.json()
            print("Full JSON response:")
            return data.get('results')[0].get('ancestors')
        else:
            print(f"Error: {response.status_code}")
            return []

    def get_relative_terms(self, go_id: str) -> List[GOTerm]:
        headers = {'Accept': 'application/json'}
        response = self._make_request("GET",
                                      f"terms/{go_id}/ancestors",
                                      headers=headers)

        if response.status_code == 200:
            data = response.json()
            go_id = data.get('results')[0].get('id')
            name = data.get('results')[0].get('name')
            parents = data.get('results')[0].get('ancestors')
            if data.get('results')[0].get('children'):
                childrens = [d['id'] for d in data.get('results')[0].get('children')]
            else:
                childrens = []
            return GOTerm(go_id=go_id, name=name, parents=parents, childrens=childrens)
        else:
            print(f"Error: {response.status_code}")
            return []

    def get_children_of_go_term(self, go_id: str) -> Dict:
        headers = {'Accept': 'application/json'}
        response = self._make_request("GET",
                                      f"terms/{go_id}/children",
                                      headers=headers)

        if response.status_code == 200:
            # Let's first see the actual structure
            data = response.json()
            results = data.get('results')
            if len(results) == 1:
                results = results[0].get('children', [])
                childrens = {}
                for result in results:
                    go_id = result['id']
                    name = result['name']
                    if go_id and name:
                        childrens[go_id] = name
                return childrens
        else:
            print(f"Error fetching data: {response.status_code}")
            return None

    def get_go_id_from_name(self, term_name: str) -> Dict:
        params = {
            'query': term_name,
            'page': 1,
            'limit': 1,
            'exact': True  # for exact match
        }
        headers = {'Accept': 'application/json'}

        response = self._make_request("GET",
                                      "search",
                                      params=params,
                                      headers=headers)

        go_dict = {}
        if response.status_code == 200:
            data = response.json()
            results = data.get('results', [])
            if results:
                # Extract both ID and name for verification
                go_dict[results[0].get('id')] = results[0].get('name')
                return go_dict
            return go_dict
        else:
            print(f"Error: {response.status_code}")
        return go_dict
