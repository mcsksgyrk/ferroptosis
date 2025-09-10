import sqlite3
from pathlib import Path
from config import OUTPUTS_DIR, PROJECT_ROOT, SOURCES_DIR
from database.sqlite_db_api3 import PsimiSQL


def merge_strings(string_1, string_2, separator="|"):
    if not string_1 and not string_2:
        return ""
    if not string_1:
        return string_2
    if not string_2:
        return string_1

    list_1 = [item for item in string_1.split(separator) if item and item != '-']
    list_2 = [item for item in string_2.split(separator) if item and item != '-']

    merged = list(set(list_1 + list_2))
    return separator.join(merged)


def main():
    SQL_SEED = PROJECT_ROOT / "database" / "network_db_seed3.sql"
    FERROPTOSIS_DB = OUTPUTS_DIR / "merged_ferroptosis_w_omnipat.db"
    ARN_DB = SOURCES_DIR / "arn" / "arn.db"
    OUTPUT_DB = OUTPUTS_DIR / "ferroptosis_autophagy.db"

    parser = PsimiSQL(SQL_SEED)
    parser.import_from_db_file(str(FERROPTOSIS_DB))

    print("Imported ferroptosis network database")

    # Build existing node lookup
    parser.cursor.execute("SELECT id, name FROM node")
    existing_nodes = {name: node_id for node_id, name in parser.cursor.fetchall()}
    print(f"Existing nodes: {len(existing_nodes)}")

    # Build existing edge lookup
    parser.cursor.execute("SELECT interactor_a_node_name, interactor_b_node_name, layer FROM edge")
    existing_edges = set()
    for a, b, layer in parser.cursor.fetchall():
        existing_edges.add((a, b, layer))
    print(f"Existing edges: {len(existing_edges)}")

    arn_db = sqlite3.connect(ARN_DB)
    arn_cursor = arn_db.cursor()

    # Process ARN nodes
    arn_cursor.execute("SELECT name, display_name, tax_id, type FROM node")
    arn_nodes = arn_cursor.fetchall()

    nodes_to_insert = []
    name_mapping = {}

    for name, display_name, tax_id, node_type in arn_nodes:
        if name in existing_nodes:
            name_mapping[name] = name
            # Update existing node to add ARN source
            parser.cursor.execute("UPDATE node SET source_db = source_db || '|ARN' WHERE name = ?", (name,))
        else:
            name_mapping[name] = name
            nodes_to_insert.append((
                name,
                'uniprot_id',
                display_name or name,
                tax_id or 9606,
                node_type or 'protein',
                '',
                '',
                '',
                'ARN'
            ))

    if nodes_to_insert:
        parser.cursor.executemany("""
            INSERT INTO node (name, primary_id_type, display_name, tax_id, type, pathways, role_in_ferroptosis, function, source_db)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, nodes_to_insert)
        parser.db.commit()

    print(f"ARN nodes: {len(nodes_to_insert)} added")

    # Rebuild node lookup after inserts
    parser.cursor.execute("SELECT id, name FROM node")
    all_nodes = {name: node_id for node_id, name in parser.cursor.fetchall()}

    # Process ARN edges in batches, only directed edges
    arn_cursor.execute("SELECT interactor_a_node_name, interactor_b_node_name, layer, interaction_types FROM edge WHERE interaction_types LIKE 'True|true%'")

    edges_to_insert = []
    edges_processed = 0
    edges_updated = 0
    batch_size = 10000

    for source_name, target_name, layer, interaction_types in arn_cursor:
        edges_processed += 1

        if edges_processed % 50000 == 0:
            print(f"Processed {edges_processed} edges...")

        if source_name not in all_nodes or target_name not in all_nodes:
            continue

        if (source_name, target_name, layer) in existing_edges:
            # Update existing edge to add ARN source
            parser.cursor.execute("""
                UPDATE edge
                SET source_db = source_db || '|ARN'
                WHERE interactor_a_node_name = ? AND interactor_b_node_name = ? AND layer = ?
            """, (source_name, target_name, layer))
            edges_updated += 1
            continue

        edges_to_insert.append((
            all_nodes[source_name],
            all_nodes[target_name],
            source_name,
            target_name,
            layer,
            interaction_types or '',
            '',
            'ARN'
        ))

        if len(edges_to_insert) >= batch_size:
            parser.cursor.executemany("""
                INSERT INTO edge (interactor_a_node_id, interactor_b_node_id, interactor_a_node_name, interactor_b_node_name, layer, interaction_types, effect_on_ferroptosis, source_db)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, edges_to_insert)
            parser.db.commit()
            print(f"Inserted batch of {len(edges_to_insert)} edges")
            edges_to_insert = []

    if edges_to_insert:
        parser.cursor.executemany("""
            INSERT INTO edge (interactor_a_node_id, interactor_b_node_id, interactor_a_node_name, interactor_b_node_name, layer, interaction_types, effect_on_ferroptosis, source_db)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, edges_to_insert)
        parser.db.commit()

    arn_db.close()

    print(f"Total ARN edges processed: {edges_processed}")
    print(f"Edges updated with ARN source: {edges_updated}")
    print(f"New edges added: {len(edges_to_insert) + (edges_processed // batch_size) * batch_size}")

    parser.save_db_to_file(str(OUTPUT_DB))
    print(f"Final merge complete: {OUTPUT_DB}")


if __name__ == "__main__":
    main()
