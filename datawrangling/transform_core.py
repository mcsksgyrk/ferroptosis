import pandas as pd
from parsers.source_parsers import KEGGPathwayParser
from typing import Dict, List
from pathlib import Path
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from apicalls.api_oop import KEGGClient, UniProtClient
from database.sqlite_db_api import PsimiSQL


def make_tables(data: Dict[int, List[str]], source: str, core: int) -> pd.DataFrame:
    df = pd.DataFrame(list(data.items()), columns=["source_id", "uniprot_id"])
    df['source'] = source
    df['is_core'] = core
    return df


def convert_kegg(name: str):
    try:
        if "cpd" in name or "dr" in name:
            client = KEGGClient()
            try:
                res = client.get_pubchem_id(name)
                return [name, "", res]
            except Exception as e:
                print(f"Failed to convert compound {name}: {str(e)}")
                return [name, "", ""]

        elif "hsa" in name:
            client = UniProtClient()
            try:
                res = client.convert_to_uniprot_id("KEGG", [name], False)
                if res and len(res) > 0 and name in res[0]:
                    return [name, res[0][name], ""]
                else:
                    print(f"No UniProt mapping found for {name}")
                    return [name, "", ""]
            except Exception as e:
                print(f"Failed to convert gene {name}: {str(e)}")
                return [name, "", ""]
        else:
            return [name, "", ""]
    except Exception as e:
        print(f"Unexpected error converting {name}: {str(e)}")
        return [name, "", ""]


kegg_parser = KEGGPathwayParser()
kegg_src = kegg_parser.read_pathway("hsa04216.xml")
kegg_edges = kegg_parser.read_edges("hsa04216.xml")
kegg_id_list = kegg_parser.extract_gene_ids(kegg_src)

rows = []
for k, v in kegg_src.items():
    if len(v) > 1:
        print(f"{k}:{len(v)}")
        for i in v:
            if 'path' in i or i == "undefined":
                continue
            else:
                kegg_id, uniprot_id, pubchem_id = convert_kegg(i)
            rows.append([k, kegg_id, uniprot_id, pubchem_id])
    else:
        if 'path' in v[0] or v[0] == "undefined":
            continue
        else:
            kegg_id, uniprot_id, pubchem_id = convert_kegg(v[0])
            rows.append([k, kegg_id, uniprot_id, pubchem_id])

kegg_node_df = pd.DataFrame(data=rows, columns=['id', 'kegg_id', 'uniprot_id', 'pubchem_id'])

edge_rows = []
for edge in kegg_edges:
    source_rows = kegg_node_df[kegg_node_df.id == edge.source_id]
    target_rows = kegg_node_df[kegg_node_df.id == edge.target_id]

    if source_rows.empty or target_rows.empty:
        print(f"Warning: No nodes found for edge {edge.source_id} -> {edge.target_id}")
        continue

    for _, source_row in source_rows.iterrows():
        for _, target_row in target_rows.iterrows():
            if "hsa" in source_row['kegg_id']:
                source_id = source_row['uniprot_id']
                source_type = "gene"
            else:
                source_id = source_row['pubchem_id']
                source_type = "compound"

            if "hsa" in target_row['kegg_id']:
                target_id = target_row['uniprot_id']
                target_type = "gene"
            else:
                target_id = target_row['pubchem_id']
                target_type = "compound"

            edge_rows.append({
                'source_id': source_id,
                'source_type': source_type,
                'target_id': target_id,
                'target_type': target_type,
                'edge_type': edge.type
            })

            print(f"Interaction: {source_id} ({source_type}) -> {target_id} ({target_type}), Type: {edge.type}")
edge_df = pd.DataFrame(edge_rows)

kegg_node_df
SQL_SEED = PROJECT_ROOT / "database" / "network_db_seed.sql"
DB_DESTINATION = OUTPUTS_DIR / "kegg.db"
db_api = PsimiSQL(SQL_SEED)

for idx, row in kegg_node_df.iterrows():
    aux_dict = dict()
    aux_dict['name'] = row.uniprot_id
    aux_dict['gene_name'] = ""
    aux_dict['tax_id'] = 9606
    aux_dict['pathways'] = ""
    aux_dict['source'] = "KEGG"
    aux_dict['function'] = ""
    db_api.insert_node(aux_dict)

for idx, row in edge_df.iterrows():
    directed = 'false'
    direct = 'flase'
    print(row.to_dict())
    if row.is_directed == 1:
        directed = 'true'
        direct = 'true'
    elif row.is_directed == 0:
        directed = 'false'
    else:
        print("WARNING: unknown direction flag in line: " + idx)
    interaction_types = "is_directed:%s|is_direct:%s" % (directed, direct)
    edge_dict = {
        'source_db': 'OmniPath',
        'interaction_types': interaction_types,
        'layer': 1,
    }
    source_dict = db_api.get_node(row.source, 9606)
    target_dict = db_api.get_node(row.target, 9606)
    db_api.insert_edge(source_dict, target_dict, edge_dict)
db_api.save_db_to_file(str(DB_DESTINATION))
