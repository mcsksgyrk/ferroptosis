from parsers.ferrdb_parser import FerrdbParser
from apicalls.mygene import MyGeneClient
from database.external_db import DBconnector
from apicalls.uniprot import UniProtClient
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
import pandas as pd
from database.sqlite_db_api3 import PsimiSQL


def get_node_dict(identifier, db_api, edge_node_dict, mygene):
    node_dict = db_api.get_node_by_any_identifier(identifier)
    if node_dict:
        return node_dict
    try:
        node_dict = db_api.get_node_by_any_identifier(edge_node_dict[identifier.lower()])
        if node_dict:
            return node_dict
    except:
        pass
    try:
        uniprot_id = mygene.query_gene(identifier)['hits'][0]['uniprot']['Swiss-Prot']
        node_dict = db_api.get_node_by_any_identifier(uniprot_id)
        if node_dict:
            return node_dict
    except:
        pass
    return None


def check_id_type(id):
    if id[0:4].lower() == "ensg":
        return 'ensambl_id'
    else:
        return 'entrez_id'


def convert_ferrdb_source():
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
        SELECT * FROM driver
        WHERE LOWER(Exp_organism) LIKE '%human%'
        AND Confidence = 'Validated'

        AND Gene_type_hgnc_locus_type_or_other = 'gene with protein product'
        """,
        'marker': """
        SELECT * FROM marker
        WHERE LOWER(Exp_organism) LIKE '%human%'
        AND Confidence = 'Validated'
        AND Gene_type_hgnc_locus_type_or_other = 'gene with protein product'
        """
    }
    nodes_df_list = []
    edges_df_list = []
    edge_node_dict = dict()
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
            for idx, row in parser.edges.iterrows():
                source = parser.get_display_name(row.source)
                source_key = row.source
                target = parser.get_display_name(row.target)
                target_key = row.target
                edge_node_dict[source_key.lower()] = source
                edge_node_dict[target_key.lower()] = target
            edges_df_list.append(parser.edges)

    final_edges = pd.concat(edges_df_list).drop_duplicates().reset_index(drop=True)
    final_nodes = pd.concat(nodes_df_list).reset_index(drop=True)

    final_nodes[final_nodes.duplicated(subset=['display_name'], keep=False)]

    source_agg = final_nodes.groupby('display_name')['source_table'].apply(lambda x: '|'.join(sorted(x.unique()))).reset_index()
    source_agg.columns = ['display_name', 'source_db']
    deduplicated_df = final_nodes.groupby('display_name').apply(lambda x: x.loc[x['uniprot_id'].notna().idxmax()], include_groups=False).reset_index()
    deduplicated_df = deduplicated_df.drop(columns=['level_1'], errors='ignore')
    final_nodes = deduplicated_df.merge(source_agg, on='display_name')

    uniprot = UniProtClient()
    symbol_nodes = final_nodes[
        (final_nodes.primary_id_type == 'Symbol') &
        (final_nodes.type != 'compound')
    ]




    mygene = MyGeneClient()
    genes_in_q = symbol_nodes.name.to_list()

    symbol_nodes.name.unique
    len(symbol_nodes)
    mygene_r = mygene.batch_query_genes(genes_in_q)
    nodes = {}
    for res in mygene_r:
        symbol = res.get('query')
        if nodes.get(symbol) and nodes.get(symbol).get('uniprot_id'):
            continue
        if nodes.get(symbol) and not res.get('uniprot'):
            continue
        node_dict = {}
        gene_id = check_id_type(res.get('_id'))
        node_dict['display_name'] = symbol
        node_dict[gene_id] = res.get('_id')
        if res.get('uniprot'):
            node_dict['uniprot_id'] = res.get('uniprot').get('Swiss-Prot')
            node_dict['primary_id_type'] = 'uniprot_id'
            node_dict['type'] = 'protein'
        else:
            node_dict['primary_id_type'] = gene_id
            node_dict['type'] = 'non_coding'
        nodes[symbol] = node_dict

    for symbol, node_dict in nodes.items():
        mask = final_nodes['display_name'] == symbol
        if node_dict.get('uniprot_id'):
            uniprot_id = node_dict['uniprot_id']
            if isinstance(uniprot_id, list):
                uniprot_id = uniprot_id[0]
            final_nodes.loc[mask, 'uniprot_id'] = uniprot_id
            final_nodes.loc[mask, 'primary_id_type'] = 'uniprot_id'
            final_nodes.loc[mask, 'type'] = 'protein'
        if node_dict.get('ensg_id'):
            final_nodes.loc[mask, 'ensg_id'] = node_dict['ensg_id']
        if node_dict.get('entrez_id'):
            final_nodes.loc[mask, 'entrez_id'] = node_dict['entrez_id']
    # final_nodes['uniprot_id'] = final_nodes.apply(lambda row: r.get(row['display_name'])
    #                                   if pd.isna(row['uniprot_id']) else row['uniprot_id'], axis=1)
    # final_nodes.loc[final_nodes['uniprot_id'].notna() & (final_nodes['uniprot_id'] != ''), 'primary_id_type'] = 'uniprot_id'
    # final_nodes.loc[final_nodes['uniprot_id'].notna() & (final_nodes['uniprot_id'] != ''), 'type'] = 'protein'
    # final_nodes['type'] = final_nodes['type'].fillna("nd")

    SQL_SEED = PROJECT_ROOT / "database" / "network_db_seed3.sql"
    DB_DESTINATION = OUTPUTS_DIR / "ferrdb_network.db"
    db_api = PsimiSQL(SQL_SEED)
    for idx, row in final_nodes.iterrows():
        node_dict = {k: (None if pd.isna(v) else v) for k, v in row.items()}
        node_dict['tax_id'] = 9606
        node_dict['type'] = node_dict.get('type') or 'nd'
        db_api.insert_node(node_dict)
        node_id = node_dict['id']
        for key, value in node_dict.items():
            if key.endswith('_id') and value is not None and value != '':
                is_primary = 1 if node_dict['primary_id_type'] == key else 0
                db_api.insert_node_identifier(node_id, key, value, is_primary)
        node_id = node_dict['id']
        for key, value in row.items():
            if key.endswith('_id') and pd.notna(value) and value != '':
                is_primary = 1 if row['primary_id_type'] == key else 0
                db_api.insert_node_identifier(node_id, key, value, is_primary)

    for idx, row in final_edges.iterrows():
        source_dict = get_node_dict(row.source, db_api, edge_node_dict, mygene)
        target_dict = get_node_dict(row.target, db_api, edge_node_dict, mygene)
        if not source_dict or not target_dict:
            print(f'Skipping {row.source} -> {row.target}: missing nodes')
            continue
        edge_type = 'activation' if row.interaction_type > 0 else 'inhibition'
        edge_dict = {
            'interactor_a_node_name': source_dict.get('id'),
            'interactor_b_node_name': target_dict.get('id'),
            'source_type': source_dict.get('type'),
            'target_type': target_dict.get('type'),
            'interaction_types': f"is_directed:true|is_direct:false|{edge_type}",
            'layer': "ferrdb_pw",
            'source_db': 'ferrdb'
        }
        db_api.insert_edge(source_dict, target_dict, edge_dict)
    db_api.save_db_to_file(str(DB_DESTINATION))
