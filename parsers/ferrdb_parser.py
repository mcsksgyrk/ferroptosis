
from apicalls.mygene import MyGeneClient
import pandas as pd
from database.external_db import DBconnector
from config import OUTPUTS_DIR, SOURCES_DIR
import re


class FerrdbParser():
    def __init__(self, df, compound_path):
        self.unicode_entities = {
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
        self.df = df
        self.compounds = self._make_compound_set(compound_path)
        self.all_entities = self._extract_entities_from_pathway()

    def _make_compound_set(self, f_path):
        compounds = set()
        with open(f_path, 'r') as f:
            for line in f:
                kegg_entry = line.strip().split('\t', 1)[1]
                names = [name.strip().lower() for name in kegg_entry.split(';')]
                compounds.update(names)
            return compounds

    def _extract_entities_from_pathway(self):
        all_entities = set()
        for idx, row in self.df.iterrows():
            pathway_str = row.Pathway
            if pd.isna(pathway_str):
                continue
            steps = [step.strip() for step in pathway_str.split(',')]
            for step in steps:
                entities = re.split(r':-:|:\+:', step)
                for entity in entities:
                    entity = entity.strip()
                    if entity:
                        all_entities.add(entity)
        return all_entities

    def extract_gene_products(self, client):
        entities = []
        for entity in self.all_entities:
            if entity[0] == '(' or entity.lower() == 'ferroptosis':
                continue
            elif self.unicode_entities.get(entity) is not None:
                entities.append(self.unicode_entities[entity])
            else:
                entities.append(entity.lower())

        res = client.batch_query_genes(entities)

        self.mygene = self._parse_mygene_response(res)

    def _parse_mygene_response(self, res):
        found_gene_symbols = dict()
        for response_data in res:
            if not response_data:
                continue
            if response_data.get('notfound'):
                continue
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
            found_gene_symbols[response_data['query']] = response_data
        return found_gene_symbols

    def pathway_to_edge(self):
        valid_nodes = self.compounds | set(self.mygene.keys())
        prev_target = None
        all_pws = []
        edge_df = pd.DataFrame(
            columns=['source', 'target', 'interaction_type']
        )
        current_pw = []
        for idx, row in self.df.iterrows():
            pw_str = row.Pathway
            if pd.isna(pw_str):
                continue
            interactions = re.split(r',|;', pw_str)
            for reaction in interactions:
                source_, target_ = [x.replace(' ', '') for x in re.split(r':\+:|:-:', reaction)]
                source = self.unicode_entities.get(source_, source_)
                target = self.unicode_entities.get(target_, target_)
                if source.lower() not in valid_nodes or target.lower() not in valid_nodes:
                    continue
                if '+' in reaction:
                    interaction_type = 1
                if '-' in reaction:
                    interaction_type = -1
                reaction_dict = {
                    'source': source,
                    'target': target,
                    'interaction_type': interaction_type
                }
                if current_pw == []:
                    prev_target = target
                if source == prev_target and reaction:
                    current_pw.append(reaction_dict)
                else:
                    if current_pw:
                        all_pws.append(current_pw)
                        pathway_df = pd.DataFrame(current_pw)
                        edge_df = pd.concat([edge_df, pathway_df], ignore_index=True)
                    current_pw = [reaction_dict]
                prev_target = target
            all_pws.append(current_pw)
        self.edges = edge_df.drop_duplicates().reset_index(drop=True)

    def parse_mygene_response(self, query):
        response_data = self.mygene.get(query)
        if not response_data:
            return None
        gene_info = {
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

    def create_node_row(self, row, table_name):
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

    def _process_rest_of_mygene(self, processed_nodes, table_name):
        nodes = []
        to_process = set(self.mygene.keys()) - processed_nodes
        print(len(to_process))
        for key in to_process:
            nodes.append(self.create_node_row(self.mygene[key], table_name))
        return nodes

    def make_nodes_df(self):
        if 'Symbol_or_reported_abbr' in self.df.columns:
            symbol_col = 'Symbol_or_reported_abbr'
        elif 'Symbol' in self.df.columns:
            symbol_col = 'Symbol'

        nodes = []
        processed_nodes = set()
        for idx, row in self.df.iterrows():
            if row[symbol_col].lower() in self.compounds:
                print(f"skipping {row[symbol_col]}")
                # kegg puchem id-khoz
                continue
            elif self.mygene.get(row[symbol_col].lower()):
                node_dict = self.create_node_row(row, "suppressor")
                node_dict['type'] = 'protein'
                nodes.append(node_dict)
                processed_nodes.add(row[symbol_col].lower())
            else:
                node_dict = self.create_node_row(row, "suppressor")
                node_dict['type'] = 'protein'
                nodes.append(node_dict)

        pw_nodes = self._process_rest_of_mygene(processed_nodes, "suppressor")
        self.nodes = pd.DataFrame(nodes+pw_nodes).drop_duplicates().reset_index(drop=True)


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
parser = FerrdbParser(df=suppressor, compound_path=f_path)
mygene = MyGeneClient()
parser.extract_gene_products(mygene)
parser.pathway_to_edge()
parser.make_nodes_df()

to_process = set(parser.mygene.keys()) - set(parser.nodes.display_name.str.lower())
print(to_process)
parser.df[parser.df.Pathway.str.contains('YY1')]
parser.nodes.to_csv('cleaner.csv')
parser.mygene.get('fas')
