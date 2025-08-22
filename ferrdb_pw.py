from apicalls.mygene import MyGeneClient
import pandas as pd
from database.external_db import DBconnector
from config import OUTPUTS_DIR, SOURCES_DIR
import re


def make_compound_set(f_path):
    compounds = set()
    with open(f_path, 'r') as f:
        for line in f:
            kegg_entry = line.strip().split('\t', 1)[1]
            names = [name.strip().lower() for name in kegg_entry.split(';')]
            compounds.update(names)
    return compounds


def create_node_row(row, table_name):
    return {
        'name': row.get('UniProtAC') or row.get('uniprot_id') or row.get('gene_name') or row.get('symbol'),
        'primary_id_type': 'uniprot_id' if (row.get('UniProtAC') or row.get('uniprot_id')) else 'Symbol',
        'display_name': row.get('Symbol') or row.get('symbol'),
        'uniprot_id': row.get('UniProtAC') or row.get('uniprot_id'),
        'hgnc_id': row.get('HGNC_ID', ''),
        'ensg_id': row.get('ENSG_stable') or row.get('ensembl_id'),
        'effect_on_ferroptosis': '',
        'source_table': 'ferrdb_'+table_name,
        'type': row.get('node_type')
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
    if not response_data:
        return None
    gene_info = {
        'entrez_id': response_data.get('_id'),
        'score': response_data.get('_score'),
        'symbol': response_data.get('symbol'),
        'name': response_data.get('name'),
        'aliases': response_data.get('alias', []),
    }
    ensembl_data = response_data.get('ensembl', {})
    if isinstance(ensembl_data, dict):
        gene_info['ensembl_id'] = ensembl_data.get('gene')
    elif isinstance(ensembl_data, list) and ensembl_data:
        gene_info['ensembl_id'] = ensembl_data[0].get('gene')
    else:
        gene_info['ensembl_id'] = None
    uniprot_data = response_data.get('uniprot', {})
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
f_path = SOURCES_DIR / "kegg/kegg_compounds.txt"
compounds = make_compound_set(f_path)

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
flat_list = list(set(item.strip() for sublist in further_processing['pathway_clean'].str.split(r':-:|:\+:|,').to_list() for item in sublist if item.strip()))
cpd_nodes = set()
for i in flat_list:
    if i.strip().lower() in compounds:
        cpd_nodes.add(i)

cleaned_flat_list = [unicode_to_standard.get(term.strip(), term.strip()) for term in flat_list]
non_cpd = list(set(cleaned_flat_list)-cpd_nodes)
non = [x for x in non_cpd if not ('(' in x or ')' in x)]
res = mygene.batch_query_genes(non)

not_genes_not_cpds = dict()
genes_dicts = dict()
for r in res:
    if 'notfound' in r.keys():
        not_genes_not_cpds[r['query']] = r
    else:
        genes_dicts[r['query']] = r

new_edges = []
all_node_names = set(all_nodes['name'].tolist())

for idx, row in further_processing.iterrows():
    pw = row.pathway
    effect_on_ferroptosis = extract_ferroptosis_effect_from_pathway(pw)
    steps = [step.strip() for step in pw.split(',')]
    for step in steps:
        if ':+:' in step:
            parts = step.split(':+:')
            if len(parts) != 2:
                continue
            source, target = parts[0].strip(), parts[1].strip()
            interaction_type = 'positive'
        elif ':-:' in step:
            parts = step.split(':-:')
            if len(parts) != 2:
                continue
            source, target = parts[0].strip(), parts[1].strip()
            interaction_type = 'negative'
        else:
            continue
        problematic_entities = set(not_genes_not_cpds.keys())
        if source in problematic_entities or target in problematic_entities:
            continue
        if target.lower() == 'ferroptosis':
            continue
        for node in [source, target]:
            if node in all_node_names:
                continue
            if node in genes_dicts.keys():
                node_row = create_node_row(genes_dicts[node], 'suppressor')
                node_row['type'] = "protein"
                node_row['effect_on_ferroptosis'] = effect_on_ferroptosis
                all_nodes = pd.concat([all_nodes, pd.DataFrame([node_row])], ignore_index=True)
                all_node_names.add(node)
            elif node.lower() in cpd_nodes:
                node_dict = {
                    'type': "compound",
                    'name': node.lower(),
                    'primary_id_type': 'compound_name',
                    'display_name': node,
                    'effect_on_ferroptosis': effect_on_ferroptosis,
                    'source_table': 'ferrdb_suppressor'
                }
                all_nodes = pd.concat([all_nodes, pd.DataFrame([node_dict])], ignore_index=True)
                all_node_names.add(node)
            else:
                continue
        if (source in genes_dicts.keys() or source.lower() in cpd_nodes) and \
           (target in genes_dicts.keys() or target.lower() in cpd_nodes):
            edge = {
                'source': source,
                'target': target,
                'interaction_type': interaction_type,
                'pathway_effect': effect_on_ferroptosis,
                'source_pathway': pw,
                'source_table': 'ferrdb_suppressor'
            }
            new_edges.append(edge)
edges_df = pd.DataFrame(new_edges)
clean_edges = edges_df.drop_duplicates(['source', 'target', 'interaction_type']).to_dict('records')
final_nodes = all_nodes.drop_duplicates()
final_nodes
