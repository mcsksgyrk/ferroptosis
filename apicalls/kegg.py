from apicalls.base import APIClient
import re
from typing import List, Dict


class KEGGClient(APIClient):
    def __init__(self):
        super().__init__("https://rest.kegg.jp")

    def get_molecule_info(self, hsa_ids: List[str]) -> Dict:
        if len(hsa_ids) == 1:
            param = hsa_ids[0]
        else:
            param = '+'.join(hsa_ids)
        response = self._make_request("GET", f"get/{param}")
        return self._parse_mol_resposne(response.text)

    def get_pubchem_id(self, mol_name: str) -> int:
        response = self._make_request("GET", f"/conv/pubchem/{mol_name}")
        res = response.text.split('pubchem:')[1].strip()
        return int(res)

    def get_symbol(self, blob: str) -> List[str]:
        text = blob.splitlines()
        for line in text:
            if "SYMBOL" in line:
                stripped_line = self._strip_white_spaces(line)
                return stripped_line.removeprefix("SYMBOL").split(',')
        return []

    def _parse_mol_resposne(self, text) -> List[Dict]:
        current_key = None
        current_values = []
        res_list = []
        res = dict()
        for line in text.split('\n'):
            if line == '///':
                res_list.append(res)
                res = dict()
                continue
            prev_key = current_key
            if line and line[0] != ' ':
                if current_key:
                    res[current_key] = current_values
                current_key = line.split()[0]
                current_values = []
                current_values.extend(self._line_parser(line, current_key))
            elif current_key == prev_key:
                current_values.extend(self._line_parser(line, current_key))
        return res_list

    @staticmethod
    def _strip_white_spaces(line: str) -> str:
        pattern = re.compile(r'\s+')
        return re.sub(pattern, '', line)

    @staticmethod
    def _line_parser(line, name):
        line_split = re.split(r'\s{2,}', line)
        if line_split[0] == '':
            line_split = line_split[1:]
        if name in line_split:
            name_index = line_split.index(name)
            values = line_split[name_index + 1:]
            return values
        else:
            return line_split
