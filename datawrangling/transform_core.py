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
            primary_id_value = pubchem_id if pubchem_id else 'cpd:'+res['ENTRY'][0]
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
                'name': primary_id_value,
                'tax_id': ''
            }
        elif 'Pathway' in res['ENTRY']:
            return {}
        elif 'Drug' in res['ENTRY']:
            pubchem_id = extract_db_id(res.get('DBLINKS', []), 'PubChem')
            primary_id_type = 'pubchem_id' if pubchem_id else 'kegg_id'
            primary_id_value = pubchem_id if pubchem_id else 'dr:'+res['ENTRY'][0]
            cpd_dict = {
                'id': id,
                'primary_id_type': primary_id_type,
                'kegg_id': 'dr:'+res['ENTRY'][0],
                'uniprot_id': '',
                'ensemble_id': '',
                'pubchem_id': pubchem_id,
                'type': 'compound',
                'display_name': res['NAME'][0].rstrip(';'),
                'function': '',
                'pathways': '',
                'role_in_ferroptosis': 'core',
                'name': primary_id_value,
                'tax_id': ''
            }
        else:
            uniprot_id = extract_db_id(res.get('DBLINKS', []), 'UniProt')
            ensembl_id = extract_db_id(res.get('DBLINKS', []), 'Ensembl')
            primary_id_value = uniprot_id if uniprot_id else res['ORGANISM'][0]+':'+res['ENTRY'][0]
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
                'name': primary_id_value,
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
# want to map to molecules, not to kegg nodes:
edge_rows = []
for edge in kegg_edges:
    source_rows = kegg_node_df[kegg_node_df.id == edge.source_id]
    target_rows = kegg_node_df[kegg_node_df.id == edge.target_id]
    if source_rows.empty or target_rows.empty:
        print(f"Warning: No nodes found for edge {edge.source_id} -> {edge.target_id}")
        continue
    source_row = source_rows.iloc[0]
    target_row = target_rows.iloc[0]

    source_primary_id_type = source_row['primary_id_type']
    source_id = source_row[source_primary_id_type]
    source_type = source_row['type']

    target_primary_id_type = target_row['primary_id_type']
    target_id = target_row[target_primary_id_type]
    target_type = target_row['type']

    is_directed_bool = 1 if edge.type != 'binding/association' else 0
    is_direct_bool = 0 if edge.type == 'repression' else 1
    directed_str = 'true' if is_directed_bool == 1 else 'false'
    direct_str = 'true' if is_direct_bool == 1 else 'false'

    if source_id and target_id:
        edge_rows.append({
            'interactor_a_node_name': source_id,
            'interactor_b_node_name': target_id,
            'source_type': source_type,
            'target_type': target_type,
            'is_directed': is_directed_bool,
            'is_direct': is_direct_bool,
            'edge_types': edge.type,
            'interaction_types': f"is_directed:{directed_str}|is_direct:{direct_str}|{edge.type}",
            'layer': 0,
            'source_db': 'KEGG'
        })
edge_df = pd.DataFrame(edge_rows).drop_duplicates()

SQL_SEED = PROJECT_ROOT / "database" / "network_db_seed3.sql"
DB_DESTINATION = OUTPUTS_DIR / "kegg.db"
db_api = PsimiSQL(SQL_SEED)

deduplicated_kegg_df = kegg_node_df.drop_duplicates()
deduplicated_kegg_df.drop('id', axis=1, inplace=True)
for idx, row in deduplicated_kegg_df.iterrows():
    node_dict = row.to_dict()
    node_dict['source_db'] = 'KEGG'
    db_api.insert_node(node_dict)

    node_id = node_dict['id']
    for key, value in row.items():
        if key.endswith('_id') and key != 'tax_id' and pd.notna(value) and value != '':
            is_primary = 1 if row['primary_id_type'] == key else 0
            db_api.insert_node_identifier(node_id, key, value, is_primary)

for idx, row in edge_df.iterrows():
    edge_dict = row.to_dict()
    source_dict = db_api.get_node(row.interactor_a_node_name)
    target_dict = db_api.get_node(row.interactor_b_node_name)
    db_api.insert_edge(source_dict, target_dict, edge_dict)
db_api.save_db_to_file(str(DB_DESTINATION))
