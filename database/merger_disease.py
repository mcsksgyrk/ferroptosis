import sqlite3
from config import OUTPUTS_DIR


def migrate_metadata():
    source_db_path = OUTPUTS_DIR / "ferreg_network.db"
    target_db_path = OUTPUTS_DIR / "merged_ferroptosis_network.db"

    source = sqlite3.connect(source_db_path)
    target = sqlite3.connect(target_db_path)
    src = source.cursor()
    tgt = target.cursor()

    # Build edge mapping: old edge id -> (source_name, target_name) from source db
    src.execute("SELECT id, interactor_a_node_name, interactor_b_node_name FROM edge")
    old_edge_names = {}
    for edge_id, a_name, b_name in src.fetchall():
        old_edge_names[edge_id] = (a_name, b_name)

    # Build disease mapping: old auto-increment id -> disease_id string from source db
    src.execute("SELECT id, disease_id FROM disease")
    old_disease_to_string = {}
    for auto_id, disease_id in src.fetchall():
        old_disease_to_string[auto_id] = disease_id

    # Build disease string -> new auto-increment id from target db
    tgt.execute("SELECT id, disease_id FROM disease")
    disease_string_to_new_id = {}
    for auto_id, disease_id in tgt.fetchall():
        disease_string_to_new_id[disease_id] = auto_id

    # Find new edge id by matching source/target names in target db
    def find_new_edge_id(a_name, b_name):
        tgt.execute(
            "SELECT id FROM edge WHERE interactor_a_node_name = ? AND interactor_b_node_name = ?",
            (a_name, b_name)
        )
        row = tgt.fetchone()
        if row:
            return row[0]
        return None

    # Migrate disease_edge
    src.execute("SELECT disease_id, edge_id, reference, source_db FROM disease_edge")
    disease_edge_rows = src.fetchall()

    disease_edge_count = 0
    disease_edge_pairs = set()

    for old_disease_auto_id, old_edge_id, reference, source_db in disease_edge_rows:
        a_name, b_name = old_edge_names.get(old_edge_id, (None, None))
        if not a_name or not b_name:
            continue

        new_edge_id = find_new_edge_id(a_name, b_name)
        if not new_edge_id:
            continue

        disease_string = old_disease_to_string.get(old_disease_auto_id)
        if not disease_string:
            continue

        new_disease_id = disease_string_to_new_id.get(disease_string)
        if not new_disease_id:
            continue

        pair = (new_disease_id, new_edge_id)
        if pair in disease_edge_pairs:
            continue
        disease_edge_pairs.add(pair)

        tgt.execute(
            "INSERT INTO disease_edge (disease_id, edge_id, reference, source_db) VALUES (?, ?, ?, ?)",
            (new_disease_id, new_edge_id, reference or '', source_db or '')
        )
        disease_edge_count += 1

    target.commit()
    print(f"Inserted {disease_edge_count} disease-edge associations")

    # Migrate experiment_model
    src.execute("SELECT edge_id, cellline, in_vivo, reference FROM experiment_model")
    experiment_rows = src.fetchall()

    experiment_count = 0

    for old_edge_id, cellline, in_vivo, reference in experiment_rows:
        a_name, b_name = old_edge_names.get(old_edge_id, (None, None))
        if not a_name or not b_name:
            continue

        new_edge_id = find_new_edge_id(a_name, b_name)
        if not new_edge_id:
            continue

        tgt.execute(
            "INSERT INTO experiment_model (edge_id, cellline, in_vivo, reference) VALUES (?, ?, ?, ?)",
            (new_edge_id, cellline or '', in_vivo or '', reference or '')
        )
        experiment_count += 1

    target.commit()
    print(f"Inserted {experiment_count} experiment models")

    source.close()
    target.close()
    print("Migration complete")
