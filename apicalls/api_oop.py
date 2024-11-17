import requests
import time
import json
import re
from typing import List, Dict, Union, Optional


class APIClient:
    def __init__(self, base_url: str, polling_interval: int = 3):
        self.base_url = base_url
        self.polling_interval = polling_interval

    def _make_request(self, method: str, endpoint: str, **kwargs) -> requests.Response:
        """Make HTTP request and handle basic error checking."""
        url = f"{self.base_url}/{endpoint}"
        response = requests.request(method, url, **kwargs)
        response.raise_for_status()
        return response


class KEGGClient(APIClient):
    """Client for interacting with KEGG API."""

    def __init__(self):
        super().__init__("https://rest.kegg.jp")

    def get_molecule_info(self, mol_name: str) -> str:
        """Retrieve molecule information from KEGG."""
        response = self._make_request("GET", f"get/{mol_name}")
        return response.text

    def get_symbol(self, blob: str) -> List[str]:
        """Extract symbol information from KEGG response."""
        text = blob.splitlines()
        for line in text:
            if "SYMBOL" in line:
                print(line)  # Consider using logging instead of print
                stripped_line = self._strip_white_spaces(line)
                return stripped_line.removeprefix("SYMBOL").split(',')
        return []

    @staticmethod
    def _strip_white_spaces(line: str) -> str:
        """Remove extra whitespace from a string."""
        pattern = re.compile(r'\s+')
        return re.sub(pattern, '', line)


class UniProtClient(APIClient):
    """Client for interacting with UniProt API."""

    def __init__(self):
        super().__init__("https://rest.uniprot.org")

    def convert_to_uniprot_id(self, db: str, ids: List[str], human: bool = True) -> tuple[Dict, List]:
        """Convert IDs to UniProt IDs."""
        return self._execute_id_mapping(
            from_db=db,
            to_db="UniProtKB-Swiss-Prot",
            ids=ids,
            human=human
        )

    def convert_from_uniprot_id(self, db: str, ids: List[str], human: bool = True) -> tuple[Dict, List]:
        """Convert UniProt IDs to other database IDs."""
        return self._execute_id_mapping(
            from_db="UniProtKB_AC-ID",
            to_db=db,
            ids=ids,
            human=human
        )

    def _execute_id_mapping(self, from_db: str, to_db: str, ids: List[str], human: bool) -> tuple[Dict, List]:
        """Execute ID mapping operation for a list of IDs."""
        result_dict = {}
        failed_ids = []

        for id_value in ids:
            print(id_value)  # Consider using logging instead of print
            job_id = self._submit_id_mapping(from_db, to_db, [id_value], human)
            results = self._get_id_mapping_results(job_id)

            try:
                result_dict[id_value] = self._parse_results(results, from_db, to_db)
            except IndexError:
                print(f'Failed to map ID {id_value} with job ID {job_id}')  # Consider using logging
                failed_ids.append(id_value)

        return result_dict, failed_ids

    def _submit_id_mapping(self, from_db: str, to_db: str, ids: List[str], human: bool) -> str:
        """Submit an ID mapping job to UniProt."""
        data = {"from": from_db, "to": to_db, "ids": ids}
        if human:
            data["taxId"] = "9606"

        response = self._make_request("POST", "idmapping/run", data=data)
        return response.json()["jobId"]

    def _get_id_mapping_results(self, job_id: str) -> Dict:
        """Poll for ID mapping results."""
        while True:
            response = self._make_request("GET", f"idmapping/status/{job_id}")
            job = response.json()

            if "jobStatus" in job:
                if job["jobStatus"] == "RUNNING":
                    print(f"Retrying in {self.polling_interval}s")  # Consider using logging
                    time.sleep(self.polling_interval)
                else:
                    raise Exception(job["jobStatus"])
            else:
                return job

    def _parse_results(self, results: Dict, from_db: str, to_db: str) -> Union[str, Dict]:
        """Parse ID mapping results based on database types."""
        if "UniProt" in to_db:
            return results['results'][0]['to']['primaryAccession']
        elif "UniProt" in from_db:
            return results['results'][0]['to']
        raise ValueError("Unsupported database combination")


class ReactomeClient(APIClient):
    """Client for interacting with Reactome API."""

    def __init__(self):
        super().__init__("https://reactome.org/ContentService")

    def map_protein_to_pathways(self, uniprot_id: str, species: str = "9606") -> List[str]:
        """Map UniProt ID to Reactome pathways."""
        response = self._make_request(
            "GET",
            f"data/mapping/UniProt/{uniprot_id}/pathways",
            params={"species": species}
        )
        return self._parse_pathway_response(response.json())

    @staticmethod
    def _parse_pathway_response(resp: Union[Dict, List]) -> List[str]:
        """Parse pathway information from Reactome response."""
        if isinstance(resp, list):
            return [item['stId'] for item in resp]
        return [resp['stId']]


uniprot_client = UniProtClient()
res, er = uniprot_client.convert_from_uniprot_id("Gene_Name", ["P42345"], human=False)
