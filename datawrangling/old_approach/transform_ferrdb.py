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
    df_cols = ['uniprot_id', "hgnc_id", "ensg_id", "pubmed_id", 'gene_name', 'topology']
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


def get_primary_id(row: pd.Series, cols: List[str], ligand: bool):
    if not ligand:
        if "uniprot_id" in cols and pd.notna(row.uniprot_id) and row.uniprot_id != '_NA_':
            return (row.uniprot_id, "uniprot_id")
        elif "ensg_id" in cols and pd.notna(row.ensg_id) and row.ensg_id != '_NA_':
            return (row.ensg_id, "ensg_id")
        elif "hgnc_id" in cols and pd.notna(row.hgnc_id) and row.hgnc_id != '_NA_':
            return (row.hgnc_id, "hgnc_id")
        elif "gene_name" in cols and pd.notna(row.gene_name) and row.gene_name != '_NA_':
            return (row.gene_name, "gene_name")
    else:
        if "pubchem_id" in cols and pd.notna(row.pubchem_id) and row.pubchem_id != '_NA_':
            return (row.pubchem_id, "pubchem_id")
        elif "molecule" in cols and pd.notna(row.molecule) and row.molecule != '_NA_':
            return (row.molecule, "molecule")
        elif "name" in cols and pd.notna(row.name) and row.name != '_NA_':
            return (row.name, "name")
    return (None, None)


def process_df(df, source_type, ligand=False):
    print(f"Processing {len(df)} {source_type}")
    df_cols = df.columns

    mol_type = 'protein'
    if 'non' in source_type.split('_'):
        mol_type = 'polynucleotide'

    for idx, row in df.iterrows():
        name_id, primary_id_type = get_primary_id(row, df_cols, ligand)
        if name_id is None:
            print(row)
        if not ligand:
            node_dict = {
                'name': name_id,
                'primary_id_type': primary_id_type,
                'display_name': row.gene_name if row.gene_name != '_NA_' else name_id,
                'tax_id': 9606 if 'Human' in str(row.topology) else None,
                'type': mol_type,
                'source': f'FerrDB_{source_type}',
                'function': '',
                'pathways': '',
            }
        else:
            node_dict = {
                'name': name_id,
                'primary_id_type': primary_id_type,
                'display_name': row.name if row.name != '_NA_' else name_id,
                'tax_id': '',
                'type': 'ligand',
                'source': f'FerrDB_{source_type}',
                'function': '',
                'pathways': '',
            }
        db_api.insert_node(node_dict)
        if not ligand:
            if node_dict.get('id'):
                if row.uniprot_id and pd.notna(row.uniprot_id):
                    db_api.insert_node_identifier(
                        node_dict['id'], 'uniprot_id',
                        row.uniprot_id,
                        (primary_id_type == 'uniprot_id')
                    )
                if row.gene_name and pd.notna(row.gene_name):
                    db_api.insert_node_identifier(
                        node_dict['id'], 'gene_name',
                        row.gene_name,
                        (primary_id_type == 'gene_name')
                    )

                if row.hgnc_id and pd.notna(row.hgnc_id):
                    db_api.insert_node_identifier(
                        node_dict['id'],
                        'hgnc_id',
                        row.hgnc_id,
                        (primary_id_type == 'hgnc_id')
                    )

                if row.ensg_id and pd.notna(row.ensg_id):
                    db_api.insert_node_identifier(
                        node_dict['id'],
                        'ensembl_id',
                        row.ensg_id,
                        (primary_id_type == 'ensg_id')
                    )
            else:
                if row.pubchem_id and pd.notna(row.pubchem_id):
                    db_api.insert_node_identifier(
                        node_dict['id'],
                        'pubchem_id',
                        row.pubchem_id,
                        (primary_id_type == 'pubchem_id')
                    )
                if row.molecule and pd.notna(row.molecule):
                    db_api.insert_node_identifier(
                        node_dict['id'],
                        'molecule',
                        row.molecule,
                        (primary_id_type == 'molecule')
                    )
                if row.name and pd.notna(row.name):
                    db_api.insert_node_identifier(
                        node_dict['id'],
                        'name',
                        row.name,
                        (primary_id_type == 'name')
                    )


ferrdb = DBconnector(OUTPUTS_DIR / "ferrdb.db")
tables = ['driver', 'suppressor', 'unclassified']
results = {}
for table in tables:
    df_prs = get_molecules_from_ferrdb_table(table, ferrdb)
    df_prs.replace('_NA_', None, inplace=True)
    results[f"{table}_prs"] = df_prs

    df_non_prs = get_molecules_from_ferrdb_table(table, ferrdb, False)
    df_non_prs.replace('_NA_', None, inplace=True)
    results[f"{table}_non_prs"] = df_non_prs

for ligand in ["inducer", "inhibitor"]:
    df_ligand = get_ligands_from_ferrdb(ligand, ferrdb)
    df_ligand.replace('_NA_', None, inplace=True)
    results[f"{ligand}"] = df_ligand

SQL_SEED = PROJECT_ROOT/"database"/"network_db_seed2.sql"
DB_DESTINATION = OUTPUTS_DIR/"ferrdb_test.db"
db_api = PsimiSQL(SQL_SEED)

for table, df in results.items():
    table_name = table.split("_")
    is_ligand = False
    if len(table_name) == 1:
        is_ligand = True
    process_df(df, table, is_ligand)

db_api.save_db_to_file(str(DB_DESTINATION))
