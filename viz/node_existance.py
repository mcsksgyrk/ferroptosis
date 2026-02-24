from apicalls.mygene import MyGeneClient
from database.external_db import DBconnector
from config import *
import networkx as nx


def get_uniprot_ids(nodes_list):
    mygene = MyGeneClient()
    uniprot_ids = []
    for gene in nodes_list:
        mygene_res = mygene.query_gene(gene)
        uniprot_hits = mygene_res['hits'][0].get('uniprot')
        if uniprot_hits:
            swissprot = uniprot_hits.get('Swiss-Prot')
            if swissprot:
                if isinstance(swissprot, list):
                    uniprot_ids.append(swissprot[0])
                else:
                    uniprot_ids.append(swissprot)

    return uniprot_ids


nodes_to_check = 'GPX4, AKT1, MTOR, EIF4EBP1, ACSL4, SLC7A11, CDKN2A, DGAT1, DGAT2, SREBF1'
nodes_list = nodes_to_check.replace(' ', '').split(',')

uniprot_ids = get_uniprot_ids(nodes_list)
db = DBconnector(OUTPUTS_DIR / 'merged_ferroptosis_w_omnipath.db')
placeholders = ','.join([f"'{uid}'" for uid in uniprot_ids])
query = f"""
    SELECT *
    FROM node
    WHERE name IN ({placeholders})
"""
db.query_to_dataframe(query)
found = set(db.query_to_dataframe(query)['name'].tolist())
missing = [uid for uid in uniprot_ids if uid not in found]
print(missing)

core_query = """
    SELECT DISTINCT n.name, n.display_name
    FROM node n
    JOIN edge e ON n.name = e.interactor_a_node_name
        OR n.name = e.interactor_b_node_name
    WHERE e.layer = '0'
    AND n.name IN ({})
""".format(','.join([f"'{uid}'" for uid in uniprot_ids]))
db.query_to_dataframe(core_query)

core_edge_query = """
    SELECT n.name, n.display_name, COUNT(e.id) as edge_count_with_core
    FROM node n
    JOIN edge e ON n.name = e.interactor_a_node_name
        OR n.name = e.interactor_b_node_name
    JOIN edge e2 ON (e.interactor_a_node_name = e2.interactor_a_node_name
        OR e.interactor_b_node_name = e2.interactor_b_node_name)
    WHERE e2.layer = '0'
    AND n.name IN ({})
    GROUP BY n.name
""".format(','.join([f"'{uid}'" for uid in uniprot_ids]))
db.query_to_dataframe(core_edge_query)

disease_query = """
    SELECT n.display_name, d.disease_name, d.disease_id, COUNT(de.id) as association_count
    FROM node n
    JOIN edge e ON n.name = e.interactor_a_node_name
        OR n.name = e.interactor_b_node_name
    JOIN disease_edge de ON e.id = de.edge_id
    JOIN disease d ON de.disease_id = d.id
    WHERE n.name IN ({})
    AND d.disease_id = 'ICD-11: 2A00'
    GROUP BY n.name, d.disease_name
    ORDER BY association_count DESC
""".format(','.join([f"'{uid}'" for uid in uniprot_ids]))
db.query_to_dataframe(disease_query)


gbm_proteins = ['P36969', 'Q9UPY5', 'O60488', 'P31749', 'P42345']
placeholders = ','.join([f"'{uid}'" for uid in gbm_proteins])

edge_query = f"""
    SELECT e.interactor_a_node_name, n1.display_name as source_name,
           e.interactor_b_node_name, n2.display_name as target_name,
           e.interaction_types, e.layer, e.source_db
    FROM edge e
    JOIN node n1 ON e.interactor_a_node_name = n1.name
    JOIN node n2 ON e.interactor_b_node_name = n2.name
    WHERE e.interactor_a_node_name IN ({placeholders})
    AND e.interactor_b_node_name IN ({placeholders})
"""
db.query_to_dataframe(edge_query)

disease_edge_query = """
    SELECT n.display_name as node,
           e.interactor_a_node_name, n1.display_name as source_name,
           e.interactor_b_node_name, n2.display_name as target_name,
           e.interaction_types, e.layer, e.source_db,
           d.disease_name, d.disease_id
    FROM node n
    JOIN edge e ON n.name = e.interactor_a_node_name
        OR n.name = e.interactor_b_node_name
    JOIN node n1 ON e.interactor_a_node_name = n1.name
    JOIN node n2 ON e.interactor_b_node_name = n2.name
    JOIN disease_edge de ON e.id = de.edge_id
    JOIN disease d ON de.disease_id = d.id
    WHERE n.name IN ({})
    AND d.disease_id = 'ICD-11: 2A00'
    ORDER BY n.display_name
""".format(','.join([f"'{uid}'" for uid in uniprot_ids]))
db.query_to_dataframe(disease_edge_query).to_csv('uhoh.csv')


edge_query = """
    SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, layer, source_db
    FROM edge
"""
edges_df = db.query_to_dataframe(edge_query)

gpx4 = 'P36969'
slc7a11 = 'Q9UPY5'

upstream_gpx4 = db.query_to_dataframe(f"""
    SELECT n.display_name, n.type, e.interaction_types, e.source_db
    FROM edge e
    JOIN node n ON e.interactor_a_node_name = n.name
    WHERE e.interactor_b_node_name = '{gpx4}'
    AND n.type != 'compound'
""")

upstream_slc7a11 = db.query_to_dataframe(f"""
    SELECT n.display_name, n.type, e.interaction_types, e.source_db
    FROM edge e
    JOIN node n ON e.interactor_a_node_name = n.name
    WHERE e.interactor_b_node_name = '{slc7a11}'
    AND n.type != 'compound'
""")

downstream_gpx4 = db.query_to_dataframe(f"""
    SELECT n.display_name, n.type, e.interaction_types, e.source_db
    FROM edge e
    JOIN node n ON e.interactor_b_node_name = n.name
    WHERE e.interactor_a_node_name = '{gpx4}'
    AND n.type != 'compound'
""")

downstream_slc7a11 = db.query_to_dataframe(f"""
    SELECT n.display_name, n.type, e.interaction_types, e.source_db
    FROM edge e
    JOIN node n ON e.interactor_b_node_name = n.name
    WHERE e.interactor_a_node_name = '{slc7a11}'
    AND n.type != 'compound'
""")

print(f"Upstream GPX4: {len(upstream_gpx4)}")
print(f"Upstream SLC7A11: {len(upstream_slc7a11)}")
print(f"Downstream GPX4: {len(downstream_gpx4)}")
print(f"Downstream SLC7A11: {len(downstream_slc7a11)}")


upstream_gpx4_names = set(upstream_gpx4['display_name'])
upstream_slc7a11_names = set(upstream_slc7a11['display_name'])

shared = upstream_gpx4_names & upstream_slc7a11_names
only_gpx4 = upstream_gpx4_names - upstream_slc7a11_names
only_slc7a11 = upstream_slc7a11_names - upstream_gpx4_names

print(f"Shared upstream regulators: {len(shared)}")
print(f"Only regulating GPX4: {len(only_gpx4)}")
print(f"Only regulating SLC7A11: {len(only_slc7a11)}")


paper_proteins = {'AKT', 'MTOR', 'CDKN2A', 'ACSL4', 'SREBF1', 'EIF4EBP1'}
print("Paper proteins in upstream GPX4:", paper_proteins & upstream_gpx4_names)
print("Paper proteins in upstream SLC7A11:", paper_proteins & upstream_slc7a11_names)

shared_regulators = db.query_to_dataframe(f"""
    SELECT DISTINCT n.display_name, n.type, n.role_in_ferroptosis,
           e1.interaction_types as effect_on_gpx4,
           e2.interaction_types as effect_on_slc7a11
    FROM node n
    JOIN edge e1 ON n.name = e1.interactor_a_node_name
        AND e1.interactor_b_node_name = '{gpx4}'
    JOIN edge e2 ON n.name = e2.interactor_a_node_name
        AND e2.interactor_b_node_name = '{slc7a11}'
    WHERE n.type != 'compound'
""")

shared_regulators.to_csv('uhoh.csv')
