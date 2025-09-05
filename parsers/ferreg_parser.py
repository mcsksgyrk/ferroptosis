import pandas as pd
from typing import List
from database.external_db import DBconnector


class FerregParser:
    def __init__(self, source_db_path):
        self.source = source_db_path
        self.db = DBconnector(self.source)
        self.interaction_df = self.get_complete_interactions()
        self.nodes = {}
        self.edges = []
        self.experiments = {}
        self.diseases = {}
        self.node_config = {
            'target': {
                'id_col': 'target_id',
                'columns': ['target_gene_name', 'target_uniprot', 'target_pathway',
                           'target_type', 'target_function', 'target_role4ferroptosis'],
                'schema_mapping': {
                    'name': ['target_uniprot', 'target_gene_name'],
                    'display_name': 'target_gene_name',
                    'tax_id': 9606,
                    'type': 'protein',
                    'pathways': 'target_pathway',
                    'role_in_ferroptosis': 'target_role4ferroptosis',
                    'function': 'target_function'
                }
            },
            'regulator': {
                'id_col': 'regulator_id',
                'columns': ['regulator_name', 'regulator_external_id', 'regulator_type',
                           'regulator_function', 'regulator2ferroptosis'],
                'schema_mapping': {
                    'name': ['regulator_external_id', 'regulator_name'],
                    'display_name': 'regulator_name',
                    'tax_id': 9606,
                    'type': 'regulator_type',
                    'pathways': '',
                    'role_in_ferroptosis': 'regulator2ferroptosis',
                    'function': 'regulator_function'
                }
            },
            'drug': {
                'id_col': 'drug_id',
                'columns': ['drug_name', 'drug_type', 'inchikey', 'ttd_id',
                           'drugmap_id', 'drug2ferroptosis'],
                'schema_mapping': {
                    'name': ['inchikey', 'drug_name'],
                    'display_name': 'drug_name',
                    'tax_id': None,
                    'type': 'small_molecule',
                    'pathways': '',
                    'role_in_ferroptosis': 'drug2ferroptosis',
                    'function': ''
                }
            }
        }

    def get_complete_interactions(self):
        query = """
        SELECT
            hub.unique_id,
            hub.target_id,
            hub.regulator_id,
            hub.drug_id,
            hub.disease_id,

            t.gene_name as target_gene_name,
            t.uniprot_id as target_uniprot,
            t.pathway as target_pathway,
            t.type as target_type,
            t.function as target_function,

            r.regulator_name,
            r.External_id as regulator_external_id,
            r.Type as regulator_type,
            r.Function as regulator_function,

            d.drug_name,
            d.drug_type,
            d.inchikey,
            d.ttd_id,
            d.drugmap_id,

            dis.Disease_name,
            dis."Disease ICD" as disease_icd,

            reg.Regulation,
            reg.drug2ferroptosis,
            reg.drug2regulator,
            reg.drug2target,
            reg."regulator to target gene",
            reg.regulator2ferroptosis,
            reg.target_role4ferroptosis,
            reg."Vivo model",
            reg.ReferenceID,
            reg.disease_detail_icd,
            reg.disease_detail,
            reg.in_vitro,
            reg."Cell Line",
            reg.pathway_id,
            reg.Cell_progress

        FROM target_regulator_drug_disease_pair hub
        LEFT JOIN general_target t ON hub.target_id = t.target_id
        LEFT JOIN general_regulator r ON hub.regulator_id = r.regulator_id
        LEFT JOIN general_drug d ON hub.drug_id = d.drug_id
        LEFT JOIN general_disease dis ON hub.disease_id = dis.disease_id
        LEFT JOIN regulation_information reg ON hub.unique_id = reg.unique_id
        WHERE hub.target_id != 'TAR99999'
        """
        return self.db.query_to_dataframe(query)

    def determine_id_type(self, node_type, col_name, value):
        if node_type == 'target':
            return 'uniprot_id'
        elif node_type == 'regulator':
            if 'name' in col_name.lower():
                return 'gene_name'
            else:
                value_str = str(value).upper()
                if (len(value_str) == 6 and value_str[0] in ['O', 'P', 'Q', 'A', 'B']):
                    return 'uniprot_id'
                elif value_str.startswith('ENSG'):
                    return 'ensembl_id'
                elif value_str.startswith('MIMAT'):
                    return 'mirbase_id'
                else:
                    return 'external_id'
        elif node_type == 'drug':
            if 'inchikey' in col_name.lower():
                return 'inchikey'
            else:
                return 'name'
        return 'name'

    def clean_value(self, value: str) -> str:
        if pd.isna(value) or value is None:
            return ''
        value_str = str(value).strip()
        empty_values = {'', '.', 'NA', 'null', 'none', 'undefined'}
        if value_str.lower() in empty_values:
            return ''
        return value_str

    def get_primary_identifier(self, node_type: str, row: pd.Series) -> List:
        config = self.node_config[node_type]
        for col in config['schema_mapping']['name']:
            value = row[col]
            if pd.notna(value) and value != '.':
                primary_id_type = self.determine_id_type(node_type, col, value)
                return value, primary_id_type
        return None, None

    def add_node(self, node_type: str, row: pd.Series):
        primary_id, primary_id_type = self.get_primary_identifier(node_type, row)
        if not primary_id:
            return
        config = self.node_config[node_type]
        node_id = row[config['id_col']]
        mapping = config['schema_mapping']
        node_dict = {
            'name': str(primary_id).strip(),
            'primary_id_type': primary_id_type,
            'display_name': self.clean_value(row[mapping['display_name']]) or str(primary_id).strip(),
            'tax_id': mapping['tax_id'],
            'type': self.get_type_value(node_type, row, mapping['type']),
            'pathways': self.clean_value(row[mapping['pathways']]) if mapping['pathways'] else '',
            'role_in_ferroptosis': self.clean_value(row[mapping['role_in_ferroptosis']]) if mapping['role_in_ferroptosis'] else '',
            'function': self.clean_value(row[mapping['function']]) if mapping['function'] else '',
            'source_db': 'ferreg',
            '_internal_id': node_id
        }
        self.nodes[node_id] = node_dict

    def get_type_value(self, node_type: str, row: pd.Series, type_field: str) -> str:
        if node_type == 'regulator':
            return self.map_regulator_type(row[type_field])
        else:
            return type_field

    def map_regulator_type(self, regulator_type: str) -> str:
        if pd.isna(regulator_type) or not regulator_type:
            return 'protein'
        regulator_type_str = str(regulator_type).lower()
        if 'protein' in regulator_type_str:
            return 'protein'
        elif 'lncrna' in regulator_type_str or 'lnc rna' in regulator_type_str:
            return 'lncRNA'
        elif 'mirna' in regulator_type_str or 'mir rna' in regulator_type_str:
            return 'miRNA'
        elif 'mrna' in regulator_type_str:
            return 'mRNA'
        else:
            return 'protein'

    def add_disease(self, row: pd.Series):
        unique_id = row['disease_id']
        if unique_id in self.diseases:
            return
        disease_icd = self.clean_value(row['disease_icd'])
        if disease_icd not in ['ICD-11: N.A.', 'N.A.', 'NA', '']:
            disease_id = disease_icd
        else:
            disease_id = row['Disease_name']
        disease_dict = {
            'disease_id': self.clean_value(disease_id),
            'disease_name': self.clean_value(row['Disease_name']),
            'description': self.clean_value(row['Regulation']),
            'unique_id': self.clean_value(row['unique_id'])
        }
        self.diseases[unique_id] = disease_dict

    def parse_edge_col_name(self, col_name: str) -> List:
        if "2" in col_name:
            return col_name.split("2")
        else:
            return ["regulator", "target"]

    def parse_edges(self, row: pd.Series):
        edges = [
            "regulator to target gene",
            "drug2target",
            "drug2regulator"
        ]
        for edge in edges:
            edge_value = self.clean_value(row[edge])
            if not edge_value:
                continue
            a_node, b_node = self.parse_edge_col_name(edge)
            a_node_id = row[a_node+"_id"]
            b_node_id = row[b_node+"_id"]
            if not self.clean_value(a_node_id) or not self.clean_value(b_node_id):
                continue
            edge_dict = {
                'interactor_a_node_name': a_node_id,
                'interactor_b_node_name': b_node_id,
                'layer': "ferreg",
                'interaction_types': row[edge],
                'effect_on_ferroptosis': "",
                'source_db': 'FerReg',
                '_unique_id_': row['unique_id']
            }
            self.edges.append(edge_dict)

    def parse_experiment(self, row: pd.Series):
        unique_id = row['unique_id']
        if unique_id not in self.experiments:
            exp_dict = {
                'cellline': row['Cell Line'],
                'in_vivo': self.clean_value(row['Vivo model']),
                'reference': 'FerReg'
            }
            self.experiments[unique_id] = exp_dict

    def parse_interactions(self):
        for idx, row in self.interaction_df.iterrows():
            nodes = [
                ('regulator', row['regulator_id']),
                ('drug', row['drug_id']),
                ('target', row['target_id'])
            ]
            for node_type, node_id in nodes:
                if pd.notna(node_id) and node_id != '.':
                    if node_id not in self.nodes:
                        self.add_node(node_type, row)

            self.add_disease(row)
            self.parse_edges(row)
            self.parse_experiment(row)
