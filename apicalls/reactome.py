from apicalls.base import APIClient
import requests
from typing import List, Dict, Union


class ReactomeClient(APIClient):
    def __init__(self):
        super().__init__("https://reactome.org/ContentService")

    def map_protein_to_pathways(self, uniprot_id: str, species: str = "9606") -> List[str]:
        try:
            response = self._make_request(
                "GET",
                f"data/mapping/UniProt/{uniprot_id}/pathways",
                params={"species": species}
            )
            return self._parse_pathway_response(response.json())
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                return []
            raise requests.HTTPError(
                f"Reactome API request failed: {str(e)}"
            ) from e

    @staticmethod
    def _parse_pathway_response(resp: Union[Dict, List]) -> List[str]:
        if isinstance(resp, list):
            return [item['stId'] for item in resp]
        return [resp['stId']]
