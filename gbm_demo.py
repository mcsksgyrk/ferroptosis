import sqlite3
import pandas as pd
from config import OUTPUTS_DIR


db_path = OUTPUTS_DIR / "merged_ferroptosis_w_omnipath.db"
conn = sqlite3.connect(db_path)


def get_neighbors(conn, node_names):
    if not node_names:
        return set()
    placeholders = ','.join(['?'] * len(node_names))
    query = f"""
        SELECT DISTINCT
            CASE
                WHEN interactor_a_node_name IN ({placeholders}) THEN interactor_b_node_name
                ELSE interactor_a_node_name
            END as neighbor
        FROM edge
        WHERE interactor_a_node_name IN ({placeholders})
           OR interactor_b_node_name IN ({placeholders})
    """
    names = list(node_names)
    df = pd.read_sql_query(query, conn, params=names * 3)
    return set(df.neighbor.tolist())


# 1. Nodes directly on glioblastoma-associated edges
direct_query = """
    SELECT DISTINCT n.name, n.display_name, n.type, n.source_db, e.layer
    FROM disease_edge de
    JOIN disease d ON de.disease_id = d.id
    JOIN edge e ON de.edge_id = e.id
    JOIN node n ON n.name IN (e.interactor_a_node_name, e.interactor_b_node_name)
    WHERE d.disease_id = 'ICD-11: 2A00'
    AND n.type = 'protein'
"""
direct_nodes = pd.read_sql_query(direct_query, conn)
direct_names = set(direct_nodes.name.tolist())
print(f"=== Direct glioblastoma-associated proteins: {len(direct_names)} ===")
print(direct_nodes.to_string(index=False))
print()


# 2. Indirect neighbors (1 intermediary)
hop1_neighbors = get_neighbors(conn, direct_names) - direct_names
hop2_neighbors = get_neighbors(conn, hop1_neighbors) - direct_names - hop1_neighbors

# filter to proteins only
if hop2_neighbors:
    placeholders = ','.join(['?'] * len(hop2_neighbors))
    hop2_query = f"""
        SELECT DISTINCT name, display_name, type, source_db
        FROM node
        WHERE name IN ({placeholders})
        AND type = 'protein'
    """
    indirect_df = pd.read_sql_query(hop2_query, conn, params=list(hop2_neighbors))
    print(f"=== Indirect GBM neighbors (via 1 intermediary): {len(indirect_df)} ===")
    print(indirect_df.head(20).to_string(index=False))
    print()
    all_gbm_candidates = direct_names | set(indirect_df.name.tolist())
else:
    all_gbm_candidates = direct_names


# 3. Distance from core ferroptosis nodes (KEGG source)
kegg_query = """
    SELECT DISTINCT name FROM node
    WHERE source_db LIKE '%KEGG%'
    AND primary_id_type = 'uniprot_id'
"""
kegg_core = set(pd.read_sql_query(kegg_query, conn).name.tolist())
print(f"=== KEGG core ferroptosis proteins: {len(kegg_core)} ===")

in_core = all_gbm_candidates & kegg_core
on_hop1 = (all_gbm_candidates - kegg_core) & get_neighbors(conn, kegg_core)
rest = all_gbm_candidates - kegg_core - on_hop1

print(f"GBM candidates IN ferroptosis core: {len(in_core)}")
for n in in_core:
    row = direct_nodes[direct_nodes.name == n]
    dn = row.iloc[0].display_name if not row.empty else n
    print(f"  {dn} ({n})")

print(f"GBM candidates 1 hop from core: {len(on_hop1)}")
print(f"GBM candidates >1 hop from core: {len(rest)}")
print()


# 4. Compounds affecting GBM candidates
if all_gbm_candidates:
    placeholders = ','.join(['?'] * len(all_gbm_candidates))
    compound_query = f"""
        SELECT DISTINCT
            n_cpd.display_name as compound,
            n_prot.display_name as target_protein,
            n_prot.name as target_id,
            e.interaction_types,
            e.layer
        FROM edge e
        JOIN node n_cpd ON n_cpd.name IN (e.interactor_a_node_name, e.interactor_b_node_name)
        JOIN node n_prot ON n_prot.name IN (e.interactor_a_node_name, e.interactor_b_node_name)
        WHERE n_cpd.type IN ('compound', 'small_molecule')
        AND n_prot.name IN ({placeholders})
        AND n_cpd.name != n_prot.name
    """
    compounds = pd.read_sql_query(compound_query, conn, params=list(all_gbm_candidates))
    print(f"=== Compounds affecting GBM candidates: {len(compounds)} ===")
    print(compounds.to_string(index=False))
    print()


conn.close()
