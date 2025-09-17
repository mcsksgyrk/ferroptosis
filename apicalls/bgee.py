from apicalls.base import APIClient
from typing import List, Dict


class BGEEClient(APIClient):
    def __init__(self):
        super().__init__("https://www.bgee.org/api/")

    def _expression_data_call(self, gene_id: str, species_id=9606) -> Dict:
        params = {
            "display_type": "json",
            "page": "data",
            "action": "expr_calls",
            "gene_id": gene_id,
            "species_id": species_id,
            "cond_param": ["anat_entity", "cell_type"],
            "data_type": "all",
            "get_results": "true",
        }
        res = self._make_request("GET", "/gene/expression", params=params)

        return res.json()

    def _get_anatEntity(self, expressionCalls: List[Dict]) -> Dict:
        gene_entity_dict = dict()
        for cal in expressionCalls:
            gene = cal['gene']['geneId']
            anat_id = cal['condition']['anatEntity']['id']
            if gene in gene_entity_dict.keys():
                gene_entity_dict[gene].append(anat_id)
            else:
                gene_entity_dict.setdefault(gene, []).append(anat_id)
        return gene_entity_dict

    def get_expression_anat_entity(self, gene_id, species_id=9606):
        res = self._expression_data_call(gene_id, species_id)
        expressionCalls = res['data']['expressionData']['expressionCalls']
        return self._get_anatEntity(expressionCalls)
