from parsers.ferrdb_parser import FerrdbParser
from apicalls.mygene import MyGeneClient
from database.external_db import DBconnector
from config import OUTPUTS_DIR, SOURCES_DIR

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
parser = FerrdbParser(df=suppressor, compound_path=f_path, table_name="suppressor")
mygene = MyGeneClient()
parser.extract_gene_products(mygene)
parser.pathway_to_edge()
parser.make_nodes_df()
parser.parse_edge_nodes()
