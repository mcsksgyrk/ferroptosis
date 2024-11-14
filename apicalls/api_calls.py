import requests
import time
import json
import re
import pickle

POLLING_INTERVAL = 3
UNIPROT_URL = "https://rest.uniprot.org"


def make_kegg_api_call(mol_name):
    url = f"https://rest.kegg.jp/get/{mol_name}"
    r = requests.get(url)
    return r.text


def get_kegg_symbol(blob):
    text = blob.splitlines()
    for line in text:
        if "SYMBOL" in line:
            print(line)
            line = strip_white_spaces(line)
            return line.removeprefix("SYMBOL").split(',')


def strip_white_spaces(line):
    pattern = re.compile(r'\s+')
    stripped_line = re.sub(pattern, '', line)
    return stripped_line


def submit_uni_id_mapping(fromDB, toDB, ids, human=True):
    if human:
        r = requests.post(
            f"{UNIPROT_URL}/idmapping/run",
            data={"from": fromDB, "to": toDB, "ids": ids, "taxId": "9606"},
            )
    else:
        r = requests.post(
            f"{UNIPROT_URL}/idmapping/run",
            data={"from": fromDB, "to": toDB, "ids": ids},
            )
    r.raise_for_status()
    return r.json()["jobId"]


def get_uni_id_mapping_results(job_id):
    while True:
        r = requests.get(f"{UNIPROT_URL}/idmapping/status/{job_id}")
        r.raise_for_status()
        job = r.json()
        if "jobStatus" in job:
            if job["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(job["jobStatus"])
        else:
            return job


def convert_to_uniprot_id(db, ids, human=True):
    fromDB = db
    toDB = "UniProtKB-Swiss-Prot"
    res, rip = execute_uniprot_api_calls(ids=ids, fromDB=fromDB,
                                         toDB=toDB, human=human)
    return res, rip


def convert_from_uniprot_id(db, ids, human=True):
    fromDB = "UniProtKB_AC-ID"
    toDB = db
    res, rip = execute_uniprot_api_calls(ids=ids, fromDB=fromDB,
                                         toDB=toDB, human=human)
    return res, rip


def parse_results(results, fromDB, toDB):
    if "UniProt" in toDB:
        return results['results'][0]['to']['primaryAccession']
    elif "UniProt" in fromDB:
        return results['results'][0]['to']


def execute_uniprot_api_calls(ids, fromDB, toDB, human=True):
    res_dict = dict()
    rippers = []
    for i in ids:
        print(i)
        job_id = submit_uni_id_mapping(
            fromDB=fromDB, toDB=toDB, ids=[i], human=human
        )
        results = get_uni_id_mapping_results(job_id)
        try:
            res_dict[i] = parse_results(results, fromDB, toDB)
        except IndexError:
            print(f'ripped id {i} with job id {job_id}')
            rippers.append(i)

    return res_dict, rippers


def get_from_reactome(uniprotID,species="Homo sapiens"):
    URL_REACTOM = "https://reactome.org/ContentService/data/mapping/UniProt/P42345/pathways?species=Homo%20sapiens"
        r = requests.post(
            f"{UNIPROT_URL}/idmapping/run",
            data={"from": fromDB, "to": toDB, "ids": ids},
            )
    r.raise_for_status()
    return r.json()["jobId"]

