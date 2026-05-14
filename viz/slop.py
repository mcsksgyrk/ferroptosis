import sqlite3
import pandas as pd
from pathlib import Path


DB_PATH = "outputs/merged_ferroptosis_network.db"

conn = sqlite3.connect(DB_PATH)


def get_node_info(conn, uniprot_id):
    query = """
        SELECT n.*, ni.id_type, ni.id_value
        FROM node n
        LEFT JOIN node_identifier ni ON n.id = ni.node_id
        WHERE n.name = ?
        OR ni.id_value = ?
    """
    return pd.read_sql_query(query, conn, params=(uniprot_id, uniprot_id))


def get_direct_neighbors(conn, uniprot_id):
    query = """
        SELECT
            e.id as edge_id,
            e.interactor_a_node_name as source,
            e.interactor_b_node_name as target,
            n_a.display_name as source_display,
            n_b.display_name as target_display,
            n_a.type as source_type,
            n_b.type as target_type,
            e.interaction_types,
            e.layer,
            e.source_db,
            e.effect_on_ferroptosis
        FROM edge e
        JOIN node n_a ON e.interactor_a_node_id = n_a.id
        JOIN node n_b ON e.interactor_b_node_id = n_b.id
        WHERE e.interactor_a_node_name = ?
        OR e.interactor_b_node_name = ?
    """
    return pd.read_sql_query(query, conn, params=(uniprot_id, uniprot_id))


def get_disease_associations(conn, uniprot_id):
    query = """
        SELECT
            d.disease_name,
            d.disease_id,
            e.interactor_a_node_name as source,
            e.interactor_b_node_name as target,
            n_a.display_name as source_display,
            n_b.display_name as target_display,
            de.source_db
        FROM disease_edge de
        JOIN disease d ON de.disease_id = d.id
        JOIN edge e ON de.edge_id = e.id
        JOIN node n_a ON e.interactor_a_node_id = n_a.id
        JOIN node n_b ON e.interactor_b_node_id = n_b.id
        WHERE e.interactor_a_node_name = ?
        OR e.interactor_b_node_name = ?
    """
    return pd.read_sql_query(query, conn, params=(uniprot_id, uniprot_id))


def get_shared_neighbors(conn, id_1, id_2):
    query = """
        SELECT neighbor FROM (
            SELECT CASE
                WHEN e.interactor_a_node_name = ? THEN e.interactor_b_node_name
                ELSE e.interactor_a_node_name
            END as neighbor
            FROM edge e
            WHERE e.interactor_a_node_name = ?
            OR e.interactor_b_node_name = ?
        )
        INTERSECT
        SELECT neighbor FROM (
            SELECT CASE
                WHEN e.interactor_a_node_name = ? THEN e.interactor_b_node_name
                ELSE e.interactor_a_node_name
            END as neighbor
            FROM edge e
            WHERE e.interactor_a_node_name = ?
            OR e.interactor_b_node_name = ?
        )
    """
    return pd.read_sql_query(query, conn, params=(id_1, id_1, id_1, id_2, id_2, id_2))


def get_neighbor_details(conn, neighbor_names):
    if not neighbor_names:
        return pd.DataFrame()
    placeholders = ",".join("?" * len(neighbor_names))
    query = f"""
        SELECT n.name, n.display_name, n.type, n.role_in_ferroptosis, n.source_db
        FROM node n
        WHERE n.name IN ({placeholders})
    """
    return pd.read_sql_query(query, conn, params=neighbor_names)


def print_section(title):
    print("\n" + "=" * 60)
    print(title)
    print("=" * 60)


CDKN2A = "P42771"
GPX4 = "P36969"


print_section("CDKN2A (P42771) - Node Info")
cdkn2a_info = get_node_info(conn, CDKN2A)
print(cdkn2a_info.to_string())

print_section("GPX4 (P36969) - Node Info")
gpx4_info = get_node_info(conn, GPX4)
print(gpx4_info.to_string())

print_section("CDKN2A - Direct Neighbors")
cdkn2a_edges = get_direct_neighbors(conn, CDKN2A)
print(f"Total edges: {len(cdkn2a_edges)}")
if not cdkn2a_edges.empty:
    print("\nNeighbors as source (CDKN2A -> X):")
    outgoing = cdkn2a_edges[cdkn2a_edges.source == CDKN2A]
    for _, row in outgoing.iterrows():
        print(f"  -> {row.target_display} ({row.target_type}) | {row.interaction_types} | layer: {row.layer} | src: {row.source_db}")

    print("\nNeighbors as target (X -> CDKN2A):")
    incoming = cdkn2a_edges[cdkn2a_edges.target == CDKN2A]
    for _, row in incoming.iterrows():
        print(f"  <- {row.source_display} ({row.source_type}) | {row.interaction_types} | layer: {row.layer} | src: {row.source_db}")

print_section("GPX4 - Direct Neighbors")
gpx4_edges = get_direct_neighbors(conn, GPX4)
print(f"Total edges: {len(gpx4_edges)}")
if not gpx4_edges.empty:
    print("\nNeighbors as source (GPX4 -> X):")
    outgoing = gpx4_edges[gpx4_edges.source == GPX4]
    for _, row in outgoing.iterrows():
        print(f"  -> {row.target_display} ({row.target_type}) | {row.interaction_types} | layer: {row.layer} | src: {row.source_db}")

    print("\nNeighbors as target (X -> GPX4):")
    incoming = gpx4_edges[gpx4_edges.target == GPX4]
    for _, row in incoming.iterrows():
        print(f"  <- {row.source_display} ({row.source_type}) | {row.interaction_types} | layer: {row.layer} | src: {row.source_db}")

print_section("Shared Neighbors (CDKN2A and GPX4)")
shared = get_shared_neighbors(conn, CDKN2A, GPX4)
if not shared.empty:
    shared_names = shared.neighbor.tolist()
    print(f"Found {len(shared_names)} shared neighbors")
    details = get_neighbor_details(conn, shared_names)
    for _, row in details.iterrows():
        print(f"  {row.display_name} ({row['name']}) | type: {row.type} | ferroptosis role: {row.role_in_ferroptosis}")
else:
    print("No shared neighbors found")

print_section("CDKN2A - Disease Associations")
cdkn2a_diseases = get_disease_associations(conn, CDKN2A)
if not cdkn2a_diseases.empty:
    for _, row in cdkn2a_diseases.iterrows():
        print(f"  {row.disease_name} ({row.disease_id}) via {row.source_display} -> {row.target_display}")
else:
    print("No disease associations found")

print_section("GPX4 - Disease Associations")
gpx4_diseases = get_disease_associations(conn, GPX4)
if not gpx4_diseases.empty:
    for _, row in gpx4_diseases.iterrows():
        print(f"  {row.disease_name} ({row.disease_id}) via {row.source_display} -> {row.target_display}")
else:
    print("No disease associations found")

print_section("GPX4 Neighbors - Compounds Only")
if not gpx4_edges.empty:
    compound_neighbors = gpx4_edges[
        (gpx4_edges.source_type == "small_molecule") |
        (gpx4_edges.target_type == "small_molecule")
    ]
    if not compound_neighbors.empty:
        for _, row in compound_neighbors.iterrows():
            if row.source_type == "small_molecule":
                print(f"  {row.source_display} -> GPX4 | {row.interaction_types}")
            else:
                print(f"  GPX4 -> {row.target_display} | {row.interaction_types}")
    else:
        print("No compound neighbors found")

print_section("CDKN2A Neighbors - Compounds Only")
if not cdkn2a_edges.empty:
    compound_neighbors = cdkn2a_edges[
        (cdkn2a_edges.source_type == "small_molecule") |
        (cdkn2a_edges.target_type == "small_molecule")
    ]
    if not compound_neighbors.empty:
        for _, row in compound_neighbors.iterrows():
            if row.source_type == "small_molecule":
                print(f"  {row.source_display} -> CDKN2A | {row.interaction_types}")
            else:
                print(f"  CDKN2A -> {row.target_display} | {row.interaction_types}")
    else:
        print("No compound neighbors found")

print_section("Direct Edge Between CDKN2A and GPX4?")
direct = conn.execute("""
    SELECT * FROM edge
    WHERE (interactor_a_node_name = ? AND interactor_b_node_name = ?)
    OR (interactor_a_node_name = ? AND interactor_b_node_name = ?)
""", (CDKN2A, GPX4, GPX4, CDKN2A)).fetchall()
if direct:
    print("YES - direct edge exists")
    for row in direct:
        print(f"  {row}")
else:
    print("NO - no direct edge between CDKN2A and GPX4")


conn.close()
print("\nDone.")
