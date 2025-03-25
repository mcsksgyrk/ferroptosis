from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from typing import List, Union, Optional, Any
import pandas as pd
from database.external_db import DBconnector
from database.sqlite_db_api2 import PsimiSQL
from apicalls.api_oop import PubChemClient


def process_target_df(df, db_api):
    print(f"Processing {len(df)} targets")
    for idx, row in df.iterrows():
        if pd.isna(row.get('gene_name')) and pd.isna(row.get('uniprot_id')):
            print(f"Skipping row with no valid identifiers: {row}")
            continue

        primary_id = row.get('uniprot_id') if pd.notna(row.get('uniprot_id')) else row.get('gene_name')
        primary_id_type = 'uniprot_id' if pd.notna(row.get('uniprot_id')) else 'gene_name'

        node_dict = {
            'name': primary_id,
            'primary_id_type': primary_id_type,
            'display_name': row.get('gene_name') if pd.notna(row.get('gene_name')) else primary_id,
            'tax_id': 9606,
            'type': 'protein',
            'source': 'ferreg_target',
            'function': str(row.get('function', '')),
            'pathways': str(row.get('pathway', ''))
        }

        node_id = db_api.insert_node(node_dict)
        if node_id:
            if pd.notna(row.get('uniprot_id')):
                db_api.insert_node_identifier(
                    node_id, 'uniprot_id',
                    row.get('uniprot_id'),
                    (primary_id_type == 'uniprot_id')
                )

            if pd.notna(row.get('gene_name')):
                db_api.insert_node_identifier(
                    node_id, 'gene_name',
                    row.get('gene_name'),
                    (primary_id_type == 'gene_name')
                )

            if pd.notna(row.get('hgnc_id')):
                db_api.insert_node_identifier(
                    node_id, 'hgnc_id',
                    row.get('hgnc_id'),
                    (primary_id_type == 'hgnc_id')
                )

            if pd.notna(row.get('kegg_id')):
                db_api.insert_node_identifier(
                    node_id, 'kegg_id',
                    row.get('kegg_id'),
                    (primary_id_type == 'kegg_id')
                )


def process_reg_df(df, db_api):
    print(f"Processing {len(df)} regulators")
    for idx, row in df.iterrows():
        # Check if valid data exists
        if pd.isna(row.get('regulator_name')) and pd.isna(row.get('external_id')):
            print(f"Skipping row with no valid identifiers: {row}")
            continue

        primary_id_type, primary_id = get_external_id_type(row.get("external_id"))

        if primary_id_type is None and pd.notna(row.get('kegg_id')):
            primary_id = row.get('kegg_id')
            primary_id_type = 'kegg_id'
        elif primary_id_type is None and pd.notna(row.get('hgnc_id')):
            primary_id = row.get('hgnc_id')
            primary_id_type = 'hgnc_id'
        elif primary_id_type is None and pd.notna(row.get('regulator_name')):
            primary_id = row.get('regulator_name')
            primary_id_type = 'gene_name'

        node_dict = {
            'name': primary_id,
            'primary_id_type': primary_id_type,
            'display_name': row.get('regulator_name') if pd.notna(row.get('regulator_name')) else primary_id,
            'tax_id': 9606,  # Assuming human
            'type': row.get('type', 'protein'),
            'source': 'ferreg_regulator',
            'function': str(row.get('function', '')),
            'pathways': str(row.get('pathway', ''))
        }
        node_id = db_api.insert_node(node_dict)
        if node_id:
            if pd.notna(row.get('regulator_name')):
                db_api.insert_node_identifier(
                    node_id, 'gene_name',  # Use gene_name instead of regulator_name
                    row.get('regulator_name'),
                    (primary_id_type == 'gene_name')
                )
            if pd.notna(row.get('external_id')):
                id_type = primary_id_type
                db_api.insert_node_identifier(
                    node_id, id_type,
                    row.get('external_id'),
                    (primary_id_type == id_type)
                )
            if pd.notna(row.get('hgnc_id')):
                db_api.insert_node_identifier(
                    node_id, 'hgnc_id',
                    row.get('hgnc_id'),
                    (primary_id_type == 'hgnc_id')
                )
            if pd.notna(row.get('kegg_id')):
                db_api.insert_node_identifier(
                    node_id, 'kegg_id',
                    row.get('kegg_id'),
                    (primary_id_type == 'kegg_id')
                )


def process_ligand_df(df, db_api):
    print(f"Processing {len(df)} ligands")
    for idx, row in df.iterrows():
        if pd.notna(row.get('sid')):
            primary_id = row.get('sid')
            primary_id_type = ''
        elif pd.notna(row.get('inchikey')):
            primary_id = row.get('inchikey')
            primary_id_type = 'inchikey'
        elif pd.notna(row.get('drug_name')):
            primary_id = row.get('drug_name')
            primary_id_type = 'name'
        else:
            print(f"Skipping row with no valid identifiers: {row}")
            continue
        node_dict = {
            'name': primary_id,
            'primary_id_type': primary_id_type,
            'display_name': row.get('name') if pd.notna(row.get('name')) else primary_id,
            'tax_id': None,
            'type': 'small_molecule',
            'source': 'ferreg_ligand',
            'function': '',
            'pathways': ''
        }
        node_id = db_api.insert_node(node_dict)
        if node_id:
            if pd.notna(row.get('sid')):
                db_api.insert_node_identifier(
                    node_id, 'sid',
                    row.get('sid'),
                    (primary_id_type == 'sid')
                )
            if pd.notna(row.get('name')):
                db_api.insert_node_identifier(
                    node_id, 'name',
                    row.get('name'),
                    (primary_id_type == 'name')
                )
            if pd.notna(row.get('inchikey')):
                db_api.insert_node_identifier(
                    node_id, 'inchikey',
                    row.get('inchikey'),
                    False
                )


def get_external_id_type(id) -> str:
    if id is None or (hasattr(id, 'isna') and id.isna()) or pd.isna(id):
        return None, None
    possibilities = {'uniprot_id': ['p', 'q', 'o', 'a', 'b'], "mirbase_id": "MIMAT", "ensembl_id": "ENS"}
    if id[0].lower() in possibilities["uniprot_id"]:
        return "uniprot_id", id
    if id[:5].lower() == possibilities["mirbase_id"].lower():
        return "mirbase_id", id
    if id[:3].lower() == possibilities["ensembl_id"].lower():
        return "ensembl_id", id
    else:
        return None, id


def extractor(db: DBconnector,
              columns: str,
              table: str,
              filters: Optional[List[str]] = None) -> List[str]:
    query = db.compile_query_string(columns, table, filters)
    res = db.custom_query(query)
    return res


ferreg = DBconnector(OUTPUTS_DIR / "ferreg.db")
query_uniques = [
    {"db": ferreg,
     "table": "general_drug",
     "id_column": "drug_name,drug_type,ttd_id,drugmap_id,inchikey",
     "source": "ferreg",
     "is_core": 0
     },

    {"db": ferreg,
     "table": "general_target",
     "id_column": "gene_name, uniprot_id, pathway, type, function, `HGNC ID`, `KEGG ID`",
     "source": "ferreg",
     "is_core": 0
     },

    {"db": ferreg,
     "table": "general_regulator",
     "id_column": "regulator_id, type, regulator_name, External_id, Family, `GENE ID`, function, `HGNC ID`, `KEGG ID`",
     "source": "ferreg",
     "is_core": 0
     }
]
ligand_cols_name = query_uniques[0]['id_column'].split(',')
res = extractor(ferreg, query_uniques[0]["id_column"], query_uniques[0]['table'])
ligand_df = pd.DataFrame(data=res, columns=ligand_cols_name)
ligand_df = ligand_df.assign(cid=pd.NA, sid=pd.NA)

pubchem_client = PubChemClient()
for row_id, row in ligand_df.iterrows():
    if row["inchikey"] == ".":
        try:
            res = pubchem_client.get_primary_cid_or_sid(row["drug_name"])
            if isinstance(res, dict):
                if "cid" in res.keys():
                    ligand_df.loc[row_id, "cid"] = res["cid"]
                if "sid" in res.keys():
                    ligand_df.loc[row_id, "sid"] = res["sid"]
                break
        except ValueError as e:
            print(f"Warning: {e}")
            continue
ligand_df = ligand_df.replace([".", "NA", ""], pd.NA)
# target_df
target_cols_name = ["gene_name", "uniprot_id",
                    "pathway", "type", "function",
                    "hgnc_id", "kegg_id"]
res_target = extractor(ferreg, query_uniques[1]["id_column"],
                       query_uniques[1]['table'])
target_df = pd.DataFrame(data=res_target, columns=target_cols_name).replace([".", "NA", ""], pd.NA)
# regulator_df
reg_cols_name = ["regulator_id", "type",
                 "regulator_name", "external_id",
                 "family", "gene_id", "function",
                 "hgnc_id", "kegg_id"]
res_reg = extractor(ferreg, query_uniques[2]["id_column"],
                    query_uniques[2]['table'])
reg_df = pd.DataFrame(data=res_reg, columns=reg_cols_name).replace([".", "NA", ""], pd.NA)

# interaction parsing:
# regulation_information itt can minden info az interactionokról, unique_id köti össze táblákat
# target_regulator_drug_disease_pair ehhez kellenek más id-k is gecco

query_edge = [
    {"db": ferreg,
     "table": "regulation_information",
     "id_column": "*",
     "source": "ferreg",
     "is_core": 0
     },

    {"db": ferreg,
     "table": "target_regulator_drug_disease_pair",
     "id_column": "*",
     "source": "ferreg",
     "is_core": 0
     },

#    {"db": ferreg,
#     "table": "general_cellline",
#     "id_column": "*",
#     "source": "ferreg",
#     "is_core": 0
#     },

    {"db": ferreg,
     "table": "general_drug",
     "id_column": "drug_id, drug_name",
     "source": "ferreg",
     "is_core": 0
     },

    {"db": ferreg,
     "table": "general_regulator",
     "id_column": "regulator_id, External_id",
     "source": "ferreg",
     "is_core": 0
     },

    {"db": ferreg,
     "table": "general_target",
     "id_column": "target_id, uniprot_id",
     "source": "ferreg",
     "is_core": 0
     }
]


def find_ferreg_id(cols, table):
    if 'unique_id' in cols:
        return 'unique_id'
    elif len(table.split('_')):
        return table.split('_')[1]+"_id"


edge_dfs = dict()
for query in query_edge:
    key = query['table']
    res = extractor(query['db'], query['id_column'], query['table'])
    if query['id_column'] == "*":
        cols = [x[1].lower() for x in ferreg.get_columns(query['table'])]
        primary_key = find_ferreg_id(cols, key)
        val = pd.DataFrame(data=res, columns=cols).set_index(primary_key)
    else:
        cols = ['ferreg_id', 'ext_id']
        val = pd.DataFrame(data=res, columns=cols).set_index('ferreg_id')
    edge_dfs[key] = val


SQL_SEED = PROJECT_ROOT/"database"/"network_db_seed2.sql"
DB_DESTINATION = OUTPUTS_DIR/"ferreg_network.db"
db_api = PsimiSQL(SQL_SEED)
process_target_df(target_df, db_api)
process_reg_df(reg_df, db_api)
process_ligand_df(ligand_df, db_api)
db_api.save_db_to_file(str(DB_DESTINATION))
