import pandas as pd
from parsers.source_parsers import KEGGPathwayParser
from typing import Dict, List
from pathlib import Path
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from apicalls.api_oop import KEGGClient, UniProtClient
from database.sqlite_db_api3 import PsimiSQL


def make_tables(data: Dict[int, List[str]], source: str, core: int) -> pd.DataFrame:
    df = pd.DataFrame(list(data.items()), columns=["source_id", "uniprot_id"])
    df['source'] = source
    df['is_core'] = core
    return df


def convert_kegg(name: str):
    if "cpd" in name or "dr" in name:
        try:
            client = KEGGClient()
            res = client.get_pubchem_id(name)
            return [name, "", res]
        except Exception as e:
            print(f"Failed to convert compound {name}: {str(e)}")
            return [name, "", name]
    elif "hsa" in name:
        client = UniProtClient()
        res = client.convert_to_uniprot_id("KEGG", [name], False)
        if res and len(res) > 0 and name in res[0]:
            return [name, res[0][name], ""]
        else:
            print(f"No UniProt mapping found for {name}")
            return [name, name, ""]
    else:
        return [name, "", ""]


kegg_parser = KEGGPathwayParser()
kegg_src = kegg_parser.read_pathway("hsa04216.xml")
kegg_src
kegg_edges = kegg_parser.read_edges("hsa04216.xml")
kegg_id_list = kegg_parser.extract_gene_ids(kegg_src)
rows = []
for k, v in kegg_src.items():
    display_name = v['display_name']

    for kegg_item in v['kegg_id']:
        if 'path' in kegg_item or kegg_item == "undefined":
            continue

        kegg_id, uniprot_id, pubchem_id = convert_kegg(kegg_item)

        if uniprot_id:
            primary_id = "uniprot_id"
        elif pubchem_id:
            primary_id = "cid"
        elif kegg_id:
            primary_id = "kegg_id"
        else:
            primary_id = "unknown"

        rows.append([k, primary_id, kegg_id, uniprot_id, pubchem_id, display_name])

kegg_node_df = pd.DataFrame(data=rows, columns=['id', 'primary_id', 'kegg_id', 'uniprot_id', 'pubchem_id', 'display_name'])
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
edge_df['is_directed'] = 1
edge_df
kegg_node_df[kegg_node_df.kegg_id.str.contains('dr')]

SQL_SEED = PROJECT_ROOT / "database" / "network_db_seed3.sql"
DB_DESTINATION = OUTPUTS_DIR / "kegg_ext2.db"
db_api = PsimiSQL(SQL_SEED)

for idx, row in kegg_node_df.iterrows():
    node_dict = dict()
    node_id_dict = dict()
    if row.uniprot_id != "":
        node_dict['name'] = row.uniprot_id
        node_dict['type'] = 'protein'
    else:
        node_dict['name'] = row.pubchem_id
        node_dict['type'] = 'small_molecule'
    node_dict['display_name'] = ""
    node_dict['primary_id_type'] = row.primary_id
    node_dict['tax_id'] = 9606 if node_dict['type'] == 'protein' else None
    node_dict['pathways'] = ""
    node_dict['source'] = "KEGG"
    node_dict['display_name'] = '|'.join(row.display_name)
    node_dict['function'] = ""
    db_api.insert_node(node_dict)

    node_id = node_dict['id']
    if row.kegg_id != "":
        if primary_id == "kegg_id":
            db_api.insert_node_identifier(node_id, 'kegg_id', row.kegg_id, True)
        else:
            db_api.insert_node_identifier(node_id, 'kegg_id', row.kegg_id)
    if row.uniprot_id != "":
        if primary_id == "uniprot_id":
            db_api.insert_node_identifier(node_id, 'uniprot_id', row.uniprot_id, True)
        else:
            db_api.insert_node_identifier(node_id, 'uniprot_id', row.uniprot_id)
    if row.pubchem_id != "":
        if primary_id == "cid":
            db_api.insert_node_identifier(node_id, 'cid', row.pubchem_id, True)
        else:
            db_api.insert_node_identifier(node_id, 'cid', row.pubchem_id, True)

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
    interaction_types = "is_directed:%s|is_direct:%s|%s" % (directed, direct, row.edge_type)
    edge_dict = {
        'source_db': 'KEGG',
        'interaction_types': interaction_types,
        'layer': 0,
    }
    source_dict = db_api.get_node(row.source_id)
    target_dict = db_api.get_node(row.target_id)
    db_api.insert_edge(source_dict, target_dict, edge_dict)
db_api.save_db_to_file(str(DB_DESTINATION))
