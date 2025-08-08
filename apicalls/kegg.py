from .base import APIClient
import re
from typing import List


class KEGGClient(APIClient):
    def __init__(self):
        super().__init__("https://rest.kegg.jp")

    def get_molecule_info(self, mol_name: str) -> str:
        response = self._make_request("GET", f"get/{mol_name}")
        return response.text

    def get_pubchem_id(self, mol_name: str) -> int:
        response = self._make_request("GET", f"/conv/pubchem/{mol_name}")
        res = response.text.split('pubchem:')[1].strip()
        return int(res)

    def get_symbol(self, blob: str) -> List[str]:
        text = blob.splitlines()
        for line in text:
            if "SYMBOL" in line:
                print(line)  # logging
                stripped_line = self._strip_white_spaces(line)
                return stripped_line.removeprefix("SYMBOL").split(',')
        return []

    @staticmethod
    def _strip_white_spaces(line: str) -> str:
        pattern = re.compile(r'\s+')
        return re.sub(pattern, '', line)
