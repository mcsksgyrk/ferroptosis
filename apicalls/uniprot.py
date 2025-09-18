from apicalls.base import APIClient
import time
import traceback
from typing import List, Dict, Union


class UniProtClient(APIClient):
    def __init__(self):
        super().__init__("https://rest.uniprot.org")

    def convert_to_uniprot_id(self, db: str, ids: List[str], human: bool = True) -> tuple[Dict, List]:
        return self._execute_id_mapping(
            from_db=db,
            to_db="UniProtKB-Swiss-Prot",
            ids=ids,
            human=human
        )

    def convert_from_uniprot_id(self, db: str, ids: List[str], human: bool = True) -> tuple[Dict, List]:
        return self._execute_id_mapping(
            from_db="UniProtKB_AC-ID",
            to_db=db,
            ids=ids,
            human=human
        )

    def batch_convert_to_uniprot_id(self, db: str, ids: List[str], batch_size: int = 25, human=False) -> tuple[Dict, List]:
        self.polling_interval = max(1.0, batch_size / 20)
        results_dict = {}
        failed_ids = []

        for i in range(0, len(ids), batch_size):
            batch = ids[i:i+batch_size]

            try:
                job_id = self._submit_id_mapping(db, "UniProtKB-Swiss-Prot", batch, human=human)
                results = self._get_id_mapping_results(job_id)
                for result in results.get('results', []):
                    from_id = result['from']
                    to_id = result['to']['primaryAccession']
                    results_dict[from_id] = to_id
                    print(f"from: {from_id}, to: {to_id}")
                if 'failedIds' in results:
                    failed_ids.extend(results['failedIds'])
            except Exception as e:
                print(f"ERROR: {e}")
                print(f"Error type: {type(e)}")
                traceback.print_exc()
                failed_ids.extend(batch)
        return results_dict, failed_ids

    def batch_convert_from_uniprot_id(self, db: str, ids: List[str], batch_size: int = 25) -> tuple[Dict, List]:
        self.polling_interval = max(1.0, batch_size / 20)
        results_dict = {}
        failed_ids = []

        for i in range(0, len(ids), batch_size):
            batch = ids[i:i+batch_size]

            try:
                job_id = self._submit_id_mapping("UniProtKB_AC-ID", db, batch, human=False)
                results = self._get_id_mapping_results(job_id)
                for result in results.get('results', []):
                    from_id = result['from']
                    if 'ensembl' in db.lower():
                        to_id = result['to'].split('.')[0]
                    else:
                        to_id = result['to']
                    results_dict[from_id] = to_id
                    print(f"from: {from_id}, to: {to_id}")
                if 'failedIds' in results:
                    failed_ids.extend(results['failedIds'])
            except Exception as e:
                print(f"ERROR: {e}")
                print(f"Error type: {type(e)}")
                traceback.print_exc()
                failed_ids.extend(batch)
        return results_dict, failed_ids

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
        result_dict = {}
        failed_ids = []

        for id_value in ids:
            print(id_value)
            job_id = self._submit_id_mapping(from_db, to_db, [id_value], human)
            results = self._get_id_mapping_results(job_id)

            try:
                result_dict[id_value] = self._parse_results(results, from_db, to_db)
            except IndexError:
                print(f'Failed to map ID {id_value} with job ID {job_id}')
                failed_ids.append(id_value)

        return result_dict, failed_ids

    def _submit_id_mapping(self, from_db: str, to_db: str, ids: List[str], human: bool) -> str:
        data = {"from": from_db, "to": to_db, "ids": ids}
        if human:
            data["taxId"] = "9606"

        response = self._make_request("POST", "idmapping/run", data=data)
        return response.json()["jobId"]

    def _get_id_mapping_results(self, job_id: str) -> Dict:
        while True:
            response = self._make_request("GET", f"idmapping/status/{job_id}")
            job = response.json()

            if "jobStatus" in job:
                if job["jobStatus"] in ("NEW", "RUNNING"):
                    print(f"Retrying in {self.polling_interval}s")
                    time.sleep(self.polling_interval)
                else:
                    raise Exception(job["jobStatus"])
            else:
                return job

    def _parse_results(self, results: Dict, from_db: str, to_db: str) -> Union[str, Dict]:
        if "UniProt" in to_db:
            return results['results'][0]['to']['primaryAccession']
        elif "UniProt" in from_db:
            return results['results'][0]['to']
        raise ValueError("Unsupported database combination")
