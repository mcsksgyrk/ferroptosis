from .base import APIClient
import requests
from typing import Dict


class PubChemClient(APIClient):
    def __init__(self):
        super().__init__("https://pubchem.ncbi.nlm.nih.gov/rest/pug")

    def name_to_cid(self, name: str) -> tuple[int, int]:
        try:
            response = self._make_request("GET", f"compound/name/{name}/cids/TXT")
            if response.text.strip():
                return int(response.text.strip().split('\n')[0]), response.status_code
            return None, response.status_code
        except requests.exceptions.HTTPError as e:
            print(f"Error retrieving CID for {name}: {e}")
            return None, getattr(e.response, "status_code", None)

    def name_to_sid(self, name: str) -> tuple[int, int]:
        try:
            response = self._make_request("GET", f"substance/name/{name}/sids/TXT")
            if response.text.strip():
                return int(response.text.strip().split('\n')[0]), response.status_code
            return None, response.status_code
        except requests.exceptions.HTTPError as e:
            print(f"Error retrieving SID for {name}: {e}")
            return None, getattr(e.response, "status_code", None)

    def get_primary_cid_or_sid(self, name: str) -> Dict[str, int]:
        cid, cid_status = self.name_to_cid(name)
        if cid_status == 200 and cid:
            return {"cid": cid}

        sid, sid_status = self.name_to_sid(name)
        if sid_status == 200 and sid:
            return {"sid": sid}

        raise ValueError(f"Neither CID nor SID found for {name}")
