from apicalls.base import APIClient
from typing import Dict, List


class MyGeneClient(APIClient):
    def __init__(self):
        super().__init__("https://mygene.info/v3")

    def query_gene(self, gene_symbol: str, species: str = "human",
                   fields: str = "symbol,name,ensembl.gene,uniprot,alias",
                   size: int = 1) -> Dict:
        params = {
            'q': gene_symbol,
            'species': species,
            'fields': fields,
            'size': size
        }

        response = self._make_request("GET", "query", params=params)
        return response.json()

    def batch_query_genes(self, gene_symbols: List[str],
                          species: str = "human",
                          fields: str = "symbol,name,ensembl.gene,uniprot,alias") -> List[Dict]:
        data = {
            'q': gene_symbols,
            'species': species,
            'fields': fields
        }
        response = self._make_request("POST", "query", json=data)
        return response.json()

    def get_gene_info(self, gene_id: str,
                      fields: str = "symbol,name,summary,ensembl") -> Dict:
        response = self._make_request("GET", f"gene/{gene_id}", params={"fields": fields})
        return response.json()
