from apicalls.mygene import MyGeneClient
import pandas as pd
from database.external_db import DBconnector
from config import OUTPUTS_DIR
import unicodedata
import re


def create_node_row(row, table_name):
    return {
        'name': row.get('UniProtAC') or row.get('gene_name'),
        'primary_id_type': 'uniprot_id' if row.get('uniprot_id') else 'Symbol',
        'display_name': row.get('Symbol', ''),
        'uniprot_id': row.get('UniProtAC', ''),
        'hgnc_id': row.get('HGNC_ID', ''),
        'ensg_id': row.get('ENSG_stable', ''),
        'effect_on_ferroptosis': '',
        'source_table': 'ferrdb_'+table_name
    }


def has_unicode(text):
    if pd.isna(text):
        return False
    return any(ord(char) > 127 for char in str(text))


def extract_entities_from_pathway(pathway_str):
    if pd.isna(pathway_str):
        return []
    all_entities = []
    steps = [step.strip() for step in pathway_str.split(',')]
    for step in steps:
        entities = step.replace(':-:', '|').replace(':+:', '|').split('|')
        for entity in entities:
            entity = entity.strip()
            if entity:
                all_entities.append(entity)
    return all_entities


def extract_unicode_entities(pathway_str):
    entities = extract_entities_from_pathway(pathway_str)
    unicode_entities = []
    for entity in entities:
        if any(ord(char) > 127 for char in entity):
            unicode_entities.append(entity)
    return unicode_entities


def sort_by_pw(df):
    df['colon_count'] = df['Pathway'].fillna('').str.count(':')
    single_step_df = df[df['colon_count'] == 2]
    multi_step_df = df[df['colon_count'] > 2]
    no_pathway_df = df[df['colon_count'] == 0]
    return single_step_df, multi_step_df, no_pathway_df


def process_single_step_pathways(single_step_df, table_name):
    node_rows = []
    for idx, row in single_step_df.iterrows():
        node_row = create_node_row(row, table_name)
        pathway = row.get('Pathway', '')
        if 'ferroptosis' in pathway.lower():
            if ':+:' in pathway:
                node_row['effect_on_ferroptosis'] = 'promotes ferroptosis'
            elif ':-:' in pathway:
                node_row['effect_on_ferroptosis'] = 'suppresses ferroptosis'
        node_rows.append(node_row)
    return pd.DataFrame(node_rows)


def extract_ferroptosis_effect_from_pathway(pathway_str):
    if pd.isna(pathway_str) or pathway_str.strip() == '':
        return ''
    steps = [step.strip() for step in pathway_str.split(',')]
    effect = 1
    for step in steps:
        if ':+:' in step:
            effect *= 1
        elif ':-:' in step:
            effect *= -1
    if effect > 0:
        return 'promotes ferroptosis'
    elif effect < 0:
        return 'suppresses ferroptosis'
    else:
        return ''


def process_multi_step_pathways(multi_step_df, table_name):
    node_rows = []
    further_processing_rows = []
    for idx, row in multi_step_df.iterrows():
        node_row = create_node_row(row, table_name)
        pathway = row.get('Pathway', '')
        node_row['effect_on_ferroptosis'] = extract_ferroptosis_effect_from_pathway(pathway)
        node_rows.append(node_row)
        further_processing_rows.append({
            'gene_symbol': row.get('Symbol'),
            'uniprot_id': row.get('UniProtAC'),
            'pathway': pathway,
            'source_table': table_name,
            'original_row_idx': idx
        })
    return pd.DataFrame(node_rows), pd.DataFrame(further_processing_rows)


def clean_whitespace_only(pathway_str):
    if pd.isna(pathway_str):
        return pathway_str
    cleaned = pathway_str.replace('\xa0', ' ')
    cleaned = re.sub(r'\s+', ' ', cleaned).strip()
    return cleaned


def process_no_pathway(no_pathway_df, table_name):
    node_rows = []
    for idx, row in no_pathway_df.iterrows():
        node_row = create_node_row(row, table_name)
        node_rows.append(node_row)
    return pd.DataFrame(node_rows)


def parse_pw(df, table_name):
    single_step_df, multi_step_df, no_pathway_df = sort_by_pw(df)
    single_nodes = process_single_step_pathways(single_step_df, table_name)
    multi_nodes, further_processing = process_multi_step_pathways(multi_step_df, table_name)
    no_pathway_nodes = process_no_pathway(no_pathway_df, table_name)
    all_nodes = pd.concat([single_nodes, multi_nodes, no_pathway_nodes], ignore_index=True)
    return all_nodes, further_processing


def parse_mygene_response(response_data):
    if not response_data or not response_data.get('hits'):
        return None
    hit = response_data['hits'][0]
    gene_info = {
        'entrez_id': hit.get('_id'),
        'score': hit.get('_score'),
        'symbol': hit.get('symbol'),
        'name': hit.get('name'),
        'aliases': hit.get('alias', []),
    }
    ensembl_data = hit.get('ensembl', {})
    if isinstance(ensembl_data, dict):
        gene_info['ensembl_id'] = ensembl_data.get('gene')
    elif isinstance(ensembl_data, list) and ensembl_data:
        gene_info['ensembl_id'] = ensembl_data[0].get('gene')
    else:
        gene_info['ensembl_id'] = None
    uniprot_data = hit.get('uniprot', {})
    if isinstance(uniprot_data, dict):
        gene_info['uniprot_id'] = uniprot_data.get('Swiss-Prot')
        if not gene_info['uniprot_id'] and uniprot_data.get('TrEMBL'):
            gene_info['uniprot_id'] = uniprot_data['TrEMBL'][0]
    else:
        gene_info['uniprot_id'] = None
    return gene_info


mygene = MyGeneClient()
ferrdb_path = OUTPUTS_DIR / "ferrdb.db"
db = DBconnector(ferrdb_path)
query_suppressor = """
SELECT * FROM suppressor
WHERE LOWER(Exp_organism) LIKE '%human%'
AND Confidence = 'Validated'
AND Gene_type_hgnc_locus_type_or_other = 'gene with protein product'
"""
suppressor = db.query_to_dataframe(query_suppressor)
all_nodes, further_processing = parse_pw(suppressor, 'suppressor')
further_processing['pathway_clean'] = further_processing['pathway'].apply(clean_whitespace_only)
further_processing['pathway_clean']
unicode_containing_genes = []
for idx, row in further_processing[further_processing['pathway_clean'].apply(has_unicode)].iterrows():
    unicodes = extract_unicode_entities(row.pathway)
    unicode_containing_genes.extend(unicodes)
manual_parse = set(unicode_containing_genes)
unicode_to_standard = {
    'GSK-3β': 'GSK3B',
    'GSK3β': 'GSK3B',
    'HIF-1α': 'HIF1A',
    'HIF1α': 'HIF1A',
    'Ikβ-α': 'NFKBIA',
    'NF-κB': 'NFKB1',
    'NFκB': 'NFKB1',
    'PPARα': 'PPARA',
    'TGF-β2': 'TGFB2',
    'β-catenin': 'CTNNB1'
}
further_processing['pathway_clean'].apply(print)
mygene.query_gene('KMT5A')
suppressor
further_processing
