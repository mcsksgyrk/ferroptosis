from apicalls.base import APIClient
from typing import List


class HPAClient(APIClient):
    def __init__(self):
        super().__init__("https://www.proteinatlas.org/api/search_download.php")

    def get_gene_info(self, id: str) -> List[str]:
        params = {
            'search': id,
            'format': 'json',
            'columns': 'g,up,eg,upbp,up_mf',
            'compress': 'no'
        }
        response = self._make_request(
            "GET",
            '',
            params=params
        )

        if response.status_code == 200:
            data = response.json()
            return data
        else:
            print(f"Error: {response.status_code}")
            return []
