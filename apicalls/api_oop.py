import requests
import time
import re
from typing import List, Dict, Union, TypedDict


class GOTerm(TypedDict):
    id: str
    name: str
    parents: List[str]
    childrens: List[str]


class APIClient:
    def __init__(self, base_url: str, polling_interval: int = 1):
        self.base_url = base_url
        self.polling_interval = polling_interval

    def _make_request(self, method: str, endpoint: str, **kwargs) -> requests.Response:
        """Make HTTP request and handle basic error checking."""
        url = f"{self.base_url}/{endpoint}"
        response = requests.request(method, url, **kwargs)
        response.raise_for_status()
        return response


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


class KEGGClient(APIClient):
    """Client for interacting with KEGG API."""

    def __init__(self):
        super().__init__("https://rest.kegg.jp")

    def get_molecule_info(self, mol_name: str) -> str:
        """Retrieve molecule information from KEGG."""
        response = self._make_request("GET", f"get/{mol_name}")
        return response.text

    def get_pubchem_id(self, mol_name: str) -> int:
        """Retrieve molecule information from KEGG."""
        response = self._make_request("GET", f"/conv/pubchem/{mol_name}")
        res = response.text.split('pubchem:')[1].strip()
        return int(res)

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

    def get_pr_fnc(self, uniprot_id: str) -> List[str]:
        fncs = []
        go_ids = []
        response = self._make_request("GET", f"/uniprotkb/{uniprot_id}.json")
        data = response.json()
        if 'uniProtKBCrossReferences' in data:
            for ref in data['uniProtKBCrossReferences']:
                if ref['database'] == 'GO':
                    for prop in ref['properties']:
                        if 'F:' in prop.get('value', ''):
                            go_id = ref.get('id')
                            term = prop.get('value', '').lower()
                            fncs.append(term)
                            if go_id:
                                go_ids.append(go_id)
        return go_ids

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
        try:
            response = self._make_request(
                "GET",
                f"data/mapping/UniProt/{uniprot_id}/pathways",
                params={"species": species}
            )
            return self._parse_pathway_response(response.json())
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                return []  # Return empty list for unknown UniProt IDs
            raise requests.HTTPError(
                f"Reactome API request failed: {str(e)}"
            ) from e

    @staticmethod
    def _parse_pathway_response(resp: Union[Dict, List]) -> List[str]:
        """Parse pathway information from Reactome response."""
        if isinstance(resp, list):
            return [item['stId'] for item in resp]
        return [resp['stId']]


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
