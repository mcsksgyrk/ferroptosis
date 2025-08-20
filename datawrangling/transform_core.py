import pandas as pd
from parsers.source_parsers import KEGGPathwayParser
from typing import Dict, List
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from apicalls.kegg import KEGGClient
from database.sqlite_db_api3 import PsimiSQL


def make_tables(data: Dict[int, List[str]], source: str, core: int) -> pd.DataFrame:
    df = pd.DataFrame(list(data.items()), columns=["source_id", "uniprot_id"])
    df['source'] = source
    df['is_core'] = core
    return df


def extract_db_id(dblinks_list, prefix):
    for link in dblinks_list:
        if link.startswith(f'{prefix}: '):
            ids = link.replace(f'{prefix}: ', '').split()
            return ids[0] if ids else ''
    return ''


def convert_kegg(res: Dict, id: int, kegg_id: str) -> Dict:
    try:
        if 'Compound' in res['ENTRY']:
            pubchem_id = extract_db_id(res.get('DBLINKS', []), 'PubChem')
            primary_id_type = 'pubchem_id' if pubchem_id else 'kegg_id'
            cpd_dict = {
                'id': id,
                'primary_id_type': primary_id_type,
                'kegg_id': 'cpd:'+res['ENTRY'][0],
                'uniprot_id': '',
                'ensemble_id': '',
                'pubchem_id': pubchem_id,
                'type': 'compound',
                'display_name': res['NAME'][0].rstrip(';'),
                'function': '',
                'pathways': '',
                'role_in_ferroptosis': 'core',
                'name': res['NAME'][0].rstrip(';'),
                'primary_id_type': primary_id_type,
                'tax_id': ''
            }
        else:
            uniprot_id = extract_db_id(res.get('DBLINKS', []), 'UniProt')
            ensembl_id = extract_db_id(res.get('DBLINKS', []), 'Ensembl')
            cpd_dict = {
                'id': id,
                'primary_id_type': 'uniprot_id',
                'kegg_id': res['ORGANISM'][0]+':'+res['ENTRY'][0],
                'uniprot_id': uniprot_id,
                'ensemble_id': ensembl_id,
                'pubchem_id': '',
                'type': 'protein',
                'display_name': res['NAME'][0] if res.get('NAME') else '',
                'function': '',
                'pathways': '',
                'role_in_ferroptosis': 'core',
                'name': res['NAME'][0] if res.get('NAME') else '',
                'tax_id': 9606
            }
        return cpd_dict
    except Exception as e:
        print(f"Failed to convert entry {kegg_id}: {str(e)}")
        return {}


kegg_parser = KEGGPathwayParser()
kegg = KEGGClient()
kegg_src = kegg_parser.read_pathway("hsa04216.xml")
kegg_edges = kegg_parser.read_edges("hsa04216.xml")
kegg_id_list = kegg_parser.extract_gene_ids(kegg_src)

rows = []
for k, v in kegg_src.items():
    display_name = v['display_name'][0] if isinstance(v['display_name'], list) else v['display_name']
    display_name = display_name.rstrip('.')
    kegg_item = v['kegg_id']

    if 'path' in kegg_item or 'undefined' in kegg_item:
        continue

    kegg_res = kegg.get_molecule_info(kegg_item)

    for node_info in kegg_res:
        df_dict = convert_kegg(node_info, k, kegg_item)
        if df_dict:
            rows.append(df_dict)

kegg_node_df = pd.DataFrame(rows)
kegg_node_df.to_csv('wtf.csv')
# want to map to molecules, not to kegg nodes:
edge_rows = []
for edge in kegg_edges:
    # Get unique molecular identifiers, not all rows
    source_rows = kegg_node_df[kegg_node_df.id == edge.source_id]
    target_rows = kegg_node_df[kegg_node_df.id == edge.target_id]
    if source_rows.empty or target_rows.empty:
        print(f"Warning: No nodes found for edge {edge.source_id} -> {edge.target_id}")
        continue
    # Take first row (since duplicates should represent same molecule)
    source_row = source_rows.iloc[0]
    target_row = target_rows.iloc[0]
    # Determine molecular identifiers
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
    # Only add edge if both IDs exist
    if source_id and target_id:
        edge_rows.append({
            'source_id': source_id,
            'source_type': source_type,
            'target_id': target_id,
            'target_type': target_type,
            'edge_type': edge.type
        })

edge_df = pd.DataFrame(edge_rows).drop_duplicates()

edge_df['is_directed'] = 1
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
