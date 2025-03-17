import pandas as pd
from typing import Dict, List
from pathlib import Path
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from apicalls.api_oop import KEGGClient, UniProtClient
from database.sqlite_db_api2 import PsimiSQL
from database.external_db import DBconnector


def get_molecules_from_ferrdb_table(table_name: str,
                                    db: DBconnector,
                                    pr=True) -> pd.DataFrame:
    columns = [tup[1] for tup in db.get_columns(table_name)]
    sym_col = ""
    for col in columns:
        if "symbol" in col.lower():
            sym_col = col
            break
    cols = ["UniProtAC", "HGNC_ID", "ENSG_stable", "PMID", sym_col, "Exp_organism"]
    if pr:
        query = (f"SELECT {', '.join(cols)} "
                 f"FROM {table_name} "
                 f"WHERE Gene_type_hgnc_locus_type_or_other='gene with protein product'")
    else:
        query = (f"SELECT {', '.join(cols)} "
                 f"FROM {table_name} "
                 f"WHERE Gene_type_hgnc_locus_type_or_other!='gene with protein product'")
    res = db.custom_query(query)
    df_cols = ['unbiprot_id', "hgnc_id", "ensg_id", "pubmed_id", 'gene_name', 'topology']
    res_df = pd.DataFrame(data=res, columns=df_cols)
    return res_df


def get_ligands_from_ferrdb(table_name, db) -> pd.DataFrame:
    cols = ["Molecule", "Name", "PubChem_CID", "PMID"]
    query = (f"SELECT {', '.join(cols)} "
             f"FROM {table_name}")
    res = db.custom_query(query)
    df_cols = ["molecule", "name", "pubchem_id", "pmid"]
    res_df = pd.DataFrame(data=res, columns=df_cols)
    return res_df


def process_proteins(df, source_type):
    print(f"Processing {len(df)} {source_type} proteins...")
    for idx, row in df.iterrows():
        # Skip entries with no valid IDs
        if row.unbiprot_id == '_NA_' and (not row.gene_name or row.gene_name == '_NA_'):
            continue
        # Determine primary ID (UniProt if available, otherwise gene name)
        if row.unbiprot_id and row.unbiprot_id != '_NA_':
            name_id = row.unbiprot_id
            primary_id_type = 'uniprot_id'
        else:
            name_id = row.gene_name.replace(' ', '_')
            primary_id_type = 'gene_name'
        # Create node dictionary
        node_dict = {
            'name': name_id,
            'primary_id_type': primary_id_type,
            'display_name': row.gene_name if row.gene_name != '_NA_' else name_id,
            'tax_id': 9606 if 'Human' in str(row.topology) else None,
            'type': 'protein',
            'source': f'FerrDB_{source_type}',
            'function': '',
            'pathways': '',
        }
        # Insert the node into our database
        db_api.insert_node(node_dict)
        # Add identifiers to node_identifier table
        if node_dict.get('id'):
            # Add UniProt ID if available
            if row.unbiprot_id and row.unbiprot_id != '_NA_':
                db_api.insert_node_identifier(
                    node_dict['id'], 'uniprot_id', row.unbiprot_id,
                    is_primary=(primary_id_type == 'uniprot_id')
                )
            # Add gene name if available
            if row.gene_name and row.gene_name != '_NA_':
                db_api.insert_node_identifier(
                    node_dict['id'], 'gene_name', row.gene_name,
                    is_primary=(primary_id_type == 'gene_name')
                )

            # Add HGNC ID if available
            if row.hgnc_id and row.hgnc_id != '_NA_':
                db_api.insert_node_identifier(node_dict['id'], 'hgnc_id', row.hgnc_id, False)

            # Add Ensembl ID if available
            if row.ensg_id and row.ensg_id != '_NA_':
                db_api.insert_node_identifier(node_dict['id'], 'ensembl_id', row.ensg_id, False)


ferrdb = DBconnector(OUTPUTS_DIR / "ferrdb.db")
tables = ['driver', 'suppressor', 'unclassified']
results = {}

for table in tables:
    results[f"{table}_prs"] = get_molecules_from_ferrdb_table(table, ferrdb)
    results[f"{table}_non_prs"] = get_molecules_from_ferrdb_table(table, ferrdb, False)

for ligand in ["inducer", "inhibitor"]:
    results[f"{ligand}"] = get_ligands_from_ferrdb(ligand, ferrdb)

SQL_SEED = PROJECT_ROOT / "database" / "network_db_seed2.sql"
DB_DESTINATION = OUTPUTS_DIR / "ferrdb_network.db"
db_api = PsimiSQL(SQL_SEED)

process_proteins(results['driver_prs'], 'driver')
process_proteins(results['suppressor_prs'], 'suppressor')
process_proteins(results['unclassified_prs'], 'unclassified')
