from parsers.ferrdb_parser import FerrdbParser
from apicalls.mygene import MyGeneClient
from database.external_db import DBconnector
from apicalls.uniprot import UniProtClient
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
import pandas as pd
from database.sqlite_db_api3 import PsimiSQL

ferrdb_path = OUTPUTS_DIR / "ferrdb.db"
f_path = SOURCES_DIR / "kegg/kegg_compounds.txt"
db = DBconnector(ferrdb_path)
mygene = MyGeneClient()

query_dict = {
    'suppressor':  """
    SELECT * FROM suppressor
    WHERE LOWER(Exp_organism) LIKE '%human%'
    AND Confidence = 'Validated'
    AND Gene_type_hgnc_locus_type_or_other = 'gene with protein product'
    """,
    'driver': """
    SELECT * FROM suppressor
    WHERE LOWER(Exp_organism) LIKE '%human%'
    AND Confidence = 'Validated'
    AND Gene_type_hgnc_locus_type_or_other = 'gene with protein product'
    """
}

nodes_df_list = []
edges_df_list = []
for k, v in query_dict.items():
    print(f"parsing {k}")
    df = db.query_to_dataframe(v)
    parser = FerrdbParser(df=df, compound_path=f_path, table_name=k)
    parser.extract_gene_products(mygene)
    parser.make_nodes_df()
    if k != 'marker':
        parser.pathway_to_edge()
        parser.parse_edge_nodes()
    if getattr(parser, "nodes", None) is not None and not parser.nodes.empty:
        nodes_df_list.append(parser.nodes)
    if getattr(parser, "edges", None) is not None and not parser.edges.empty:
        edges_df_list.append(parser.edges)

final_edges = pd.concat(edges_df_list).drop_duplicates().reset_index(drop=True)
final_nodes = pd.concat(nodes_df_list).reset_index(drop=True)
parser.df
parser.nodes
edge_nodes = set(parser.edges['target'].to_list()).union(set(parser.edges['source'].to_list()))
node_nodes = set(parser.nodes.display_name.to_list())
yikes = edge_nodes - node_nodes
tp = parser


def testing_edge_nodes(k):
    unique_sources = k.edges.source.unique()
    unique_target = k.edges.target.unique()
    unique_edge_nodes = set(unique_sources) | set(unique_target)
    only_edge = unique_edge_nodes-set(k.nodes.display_name)
    nodes = []
    for oe in only_edge:
        if oe.lower() in k.compounds:
            node_dict = k.create_node_row({'symbol': oe})
            node_dict['type'] = "compound"
            nodes.append(node_dict)
        elif k.mygene.get(oe.lower()) or k.mygene.get(oe) or k.mygene.get(oe.upper()):
            gene_data = (k.mygene.get(oe.lower()) or
                         k.mygene.get(oe) or
                         k.mygene.get(oe.upper()))
            node_dict = k.create_node_row(gene_data)
            nodes.append(node_dict)
    reject_nodes = pd.DataFrame(nodes)
    if getattr(k, "nodes", None) is not None and not k.nodes.empty:
        k.nodes = pd.concat([reject_nodes, k.nodes]).drop_duplicates().reset_index(drop=True)
    else:
        k.nodes = reject_nodes
    return only_edge


gene_data = tp.mygene.get('yy1')
node_dict = tp.create_node_row(gene_data)
oe = testing_edge_nodes(tp)
tp_eo = set(tp.edges['target'].str.lower().to_list()).union(set(parser.edges['source'].str.lower().to_list()))
tp_nodes = set(tp.nodes.display_name.str.lower().to_list())
tp_yikes = tp_eo - tp_nodes
len(tp_yikes)
tp_yikes
tp.nodes[tp.nodes['type']=='compound']
'cbp'.lower() in tp.compounds
tp.mygene.get('yy1')
nombres = []
for oe in tp_yikes:
    gene_data = (tp.mygene.get(oe.lower()) or
                 tp.mygene.get(oe) or
                 tp.mygene.get(oe.upper()))
    node_dict = tp.create_node_row(gene_data)
    nombres.append(node_dict['display_name'])
