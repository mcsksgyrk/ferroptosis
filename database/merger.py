import sqlite3
import os
from pathlib import Path
from config import OUTPUTS_DIR, PROJECT_ROOT
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


def get_union_of_nodes(node_1, node_2):
    # MODIFICATION: Updated for ferroptosis schema
    return {
        "name": node_1["name"],
        "primary_id_type": node_1["primary_id_type"],
        "display_name": node_1["display_name"],
        "tax_id": node_1["tax_id"],
        "type": node_1["type"],
        "pathways": merge_strings(node_1["pathways"], node_2["pathways"]),
        "role_in_ferroptosis": merge_strings(node_1["role_in_ferroptosis"], node_2["role_in_ferroptosis"]),
        "function": merge_strings(node_1["function"], node_2["function"]),
        "source_db": merge_strings(node_1["source_db"], node_2["source_db"])
    }


def main():
    SQL_SEED = PROJECT_ROOT / "database" / "network_db_seed3.sql"
    OUTPUT_DB = OUTPUTS_DIR / "merged_ferroptosis_network.db"

    input_dbs = [
        OUTPUTS_DIR / "kegg.db",
        OUTPUTS_DIR / "ferrdb_network.db",
        OUTPUTS_DIR / "ferreg_network.db"
    ]

    nodes = {}
    all_node_data = []  # Store all raw node data first
    all_identifiers = {}  # identifier_value -> list of node_data
    diseases = {}
    collected_edges_directed = {}
    collected_edges_undirected = {}
    disease_edges = []
    experiment_models = []

    # First pass: collect all raw data
    for db_file in input_dbs:
        if not os.path.exists(db_file):
            continue

        print(f"Collecting raw data from: {db_file}")
        db = sqlite3.connect(db_file)
        cursor = db.cursor()

        # Collect all nodes and their identifiers
        cursor.execute("SELECT * FROM node")
        node_id_mapping = {}

        for row in cursor.fetchall():
            id, name, primary_id_type, display_name, tax_id, type, pathways, role_in_ferroptosis, function, source_db = row

            node_data = {
                "db_id": id,
                "name": name,
                "primary_id_type": primary_id_type,
                "display_name": display_name,
                "tax_id": tax_id,
                "type": type,
                "pathways": pathways or "",
                "role_in_ferroptosis": role_in_ferroptosis or "",
                "function": function or "",
                "source_db": source_db or "",
                "db_file": str(db_file),
                "identifiers": {}
            }

            all_node_data.append(node_data)
            node_id_mapping[id] = len(all_node_data) - 1  # index in all_node_data

        # Collect identifiers for these nodes
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='node_identifier'")
        if cursor.fetchone():
            cursor.execute("SELECT * FROM node_identifier")
            for row in cursor.fetchall():
                node_db_id, id_type, is_primary, id_value = row
                node_index = node_id_mapping.get(node_db_id)
                if node_index is not None:
                    all_node_data[node_index]["identifiers"][id_type] = {
                        'value': id_value,
                        'is_primary': is_primary
                    }

        # Collect other data (diseases, edges, etc.) as before
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='disease'")
        if cursor.fetchone():
            cursor.execute("SELECT * FROM disease")
            for row in cursor.fetchall():
                id, disease_id, disease_name, description = row
                if disease_id not in diseases:
                    diseases[disease_id] = {
                        'disease_id': disease_id,
                        'disease_name': disease_name or "",
                        'description': description or ""
                    }
                else:
                    diseases[disease_id]['description'] = merge_strings(diseases[disease_id]['description'], description or "")

        # Store edge info for later relationship mapping
        edge_info = {}
        cursor.execute("SELECT * FROM edge")
        for row in cursor.fetchall():
            edge_db_id = row[0]
            source_name = row[3]
            target_name = row[4]
            layer = row[5]
            edge_info[edge_db_id] = {
                'source_name': source_name,
                'target_name': target_name,
                'layer': layer
            }

        # Collect disease edges and experiment models as before
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='disease_edge'")
        if cursor.fetchone():
            disease_mapping = {}
            cursor.execute("SELECT id, disease_id FROM disease")
            for row in cursor.fetchall():
                disease_mapping[row[0]] = row[1]

            cursor.execute("SELECT * FROM disease_edge")
            for row in cursor.fetchall():
                de_id, disease_db_id, edge_db_id, reference, source_db = row
                disease_id = disease_mapping.get(disease_db_id)
                edge_data = edge_info.get(edge_db_id)

                if disease_id and edge_data:
                    disease_edges.append({
                        'disease_id': disease_id,
                        'source_name': edge_data['source_name'],
                        'target_name': edge_data['target_name'],
                        'layer': edge_data['layer'],
                        'reference': reference,
                        'source_db': source_db
                    })

        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='experiment_model'")
        if cursor.fetchone():
            cursor.execute("SELECT * FROM experiment_model")
            for row in cursor.fetchall():
                edge_db_id, cellline, in_vivo, reference = row
                edge_data = edge_info.get(edge_db_id)

                if edge_data:
                    experiment_models.append({
                        'source_name': edge_data['source_name'],
                        'target_name': edge_data['target_name'],
                        'layer': edge_data['layer'],
                        'cellline': cellline,
                        'in_vivo': in_vivo,
                        'reference': reference
                    })

        db.close()

    # Second pass: build identifier mapping and detect duplicates
    print("Building identifier mapping and detecting duplicates...")

    # Create mapping: identifier_value -> list of node indices
    for i, node_data in enumerate(all_node_data):
        # Add primary name
        identifier_value = node_data["name"]
        if identifier_value not in all_identifiers:
            all_identifiers[identifier_value] = []
        all_identifiers[identifier_value].append(i)

        # Add all other identifiers
        for id_type, id_info in node_data["identifiers"].items():
            identifier_value = id_info["value"]
            if identifier_value and identifier_value not in all_identifiers:
                all_identifiers[identifier_value] = []
            if identifier_value:
                all_identifiers[identifier_value].append(i)

    # Group nodes that share any identifier
    node_groups = []
    processed = set()

    for i, node_data in enumerate(all_node_data):
        if i in processed:
            continue

        # Find all nodes that share identifiers with this one
        group = set([i])
        to_check = [i]

        while to_check:
            current_idx = to_check.pop()
            current_node = all_node_data[current_idx]

            # Check all identifiers of current node
            all_ids = [current_node["name"]] + [id_info["value"] for id_info in current_node["identifiers"].values() if id_info["value"]]

            for identifier_value in all_ids:
                if identifier_value in all_identifiers:
                    for related_idx in all_identifiers[identifier_value]:
                        if related_idx not in group:
                            group.add(related_idx)
                            to_check.append(related_idx)

        node_groups.append(list(group))
        processed.update(group)

    # Create merged nodes
    node_identifiers = {}
    for group in node_groups:
        # Pick the first node as base, merge others into it
        base_node = all_node_data[group[0]]
        merged_node = {
            "name": base_node["name"],
            "primary_id_type": base_node["primary_id_type"],
            "display_name": base_node["display_name"],
            "tax_id": base_node["tax_id"],
            "type": base_node["type"],
            "pathways": base_node["pathways"],
            "role_in_ferroptosis": base_node["role_in_ferroptosis"],
            "function": base_node["function"],
            "source_db": base_node["source_db"]
        }

        # Merge all nodes in group
        merged_identifiers = dict(base_node["identifiers"])

        for idx in group[1:]:
            other_node = all_node_data[idx]
            merged_node = get_union_of_nodes(merged_node, other_node)

            # Merge identifiers
            for id_type, id_info in other_node["identifiers"].items():
                if id_type not in merged_identifiers:
                    merged_identifiers[id_type] = id_info

        # Use the primary name as the key
        canonical_name = merged_node["name"]
        nodes[canonical_name] = merged_node
        node_identifiers[canonical_name] = merged_identifiers

    print(f"Total raw nodes collected: {len(all_node_data)}")
    print(f"Total unique nodes after deduplication: {len(nodes)}")
    print(f"Total unique diseases: {len(diseases)}")

    # Create output database and insert nodes
    parser = PsimiSQL(SQL_SEED)
    for node in nodes.values():
        parser.insert_unique_node(node)
        nodes[node['name']]['id'] = parser.cursor.lastrowid

    # Insert node identifiers
    for node_name, identifiers in node_identifiers.items():
        if node_name in nodes:
            node_id = nodes[node_name]['id']
            for id_type, id_data in identifiers.items():
                parser.insert_node_identifier(node_id, id_type, id_data['value'], id_data['is_primary'])

    # Insert diseases
    for disease in diseases.values():
        parser.insert_disease(disease)
        diseases[disease['disease_id']]['id'] = parser.cursor.lastrowid

    # Second pass: collect edges
    all_edge_counter = 0
    merged_edge_counter = 0

    for db_file in input_dbs:
        if not os.path.exists(db_file):
            continue

        print(f"Processing edges from: {db_file}")
        db = sqlite3.connect(db_file)
        cursor = db.cursor()
        cursor.execute("SELECT * FROM edge")

        for row in cursor.fetchall():
            all_edge_counter += 1
            # MODIFICATION: Updated for new schema
            edge_id, interactor_a_node_id, interactor_b_node_id, interactor_a_node_name, interactor_b_node_name, layer, interaction_types, effect_on_ferroptosis, source_db = row

            if interactor_a_node_name not in nodes or interactor_b_node_name not in nodes:
                continue

            edge_key = f"{interactor_a_node_name}@{interactor_b_node_name}@{layer}"
            edge_key_reverse = f"{interactor_b_node_name}@{interactor_a_node_name}@{layer}"

            current_edge = {
                'interaction_types': interaction_types or "",
                'effect_on_ferroptosis': effect_on_ferroptosis or "",
                'source_db': source_db or "",
                'layer': layer
            }

            # Check if directed
            directed = interaction_types and 'is_directed:true' in interaction_types.split('|')

            if directed:
                if edge_key not in collected_edges_directed:
                    collected_edges_directed[edge_key] = current_edge
                else:
                    existing = collected_edges_directed[edge_key]
                    existing['interaction_types'] = merge_strings(existing['interaction_types'], current_edge['interaction_types'])
                    existing['effect_on_ferroptosis'] = merge_strings(existing['effect_on_ferroptosis'], current_edge['effect_on_ferroptosis'])
                    existing['source_db'] = merge_strings(existing['source_db'], current_edge['source_db'])
                    merged_edge_counter += 1
            else:
                if edge_key not in collected_edges_undirected and edge_key_reverse not in collected_edges_undirected:
                    collected_edges_undirected[edge_key] = current_edge
                else:
                    existing_key = edge_key if edge_key in collected_edges_undirected else edge_key_reverse
                    existing = collected_edges_undirected[existing_key]
                    existing['interaction_types'] = merge_strings(existing['interaction_types'], current_edge['interaction_types'])
                    existing['effect_on_ferroptosis'] = merge_strings(existing['effect_on_ferroptosis'], current_edge['effect_on_ferroptosis'])
                    existing['source_db'] = merge_strings(existing['source_db'], current_edge['source_db'])
                    merged_edge_counter += 1

        db.close()

    print(f"Total edges read: {all_edge_counter}")
    print(f"Unique edges to write: {len(collected_edges_directed) + len(collected_edges_undirected)}")
    print(f"Merged edges: {merged_edge_counter}")

    # Insert edges
    edge_mapping = {}  # Track edge keys to new edge IDs
    for edge_key, edge_data in {**collected_edges_undirected, **collected_edges_directed}.items():
        node_a, node_b, layer = edge_key.split('@')
        parser.insert_edge(nodes[node_a], nodes[node_b], edge_data)
        edge_id = parser.cursor.lastrowid
        edge_mapping[edge_key] = edge_id

    # Insert experiment models
    print(f"Processing {len(experiment_models)} experiment models")
    for exp in experiment_models:
        edge_key = f"{exp['source_name']}@{exp['target_name']}@{exp['layer']}"
        edge_key_reverse = f"{exp['target_name']}@{exp['source_name']}@{exp['layer']}"

        edge_id = edge_mapping.get(edge_key) or edge_mapping.get(edge_key_reverse)
        if edge_id:
            exp_dict = {
                'edge_id': edge_id,
                'cellline': exp['cellline'],
                'in_vivo': exp['in_vivo'],
                'reference': exp['reference']
            }
            parser.insert_experiment_model(exp_dict)

    # Insert disease edges
    print(f"Processing {len(disease_edges)} disease edges")
    for de in disease_edges:
        edge_key = f"{de['source_name']}@{de['target_name']}@{de['layer']}"
        edge_key_reverse = f"{de['target_name']}@{de['source_name']}@{de['layer']}"

        edge_id = edge_mapping.get(edge_key) or edge_mapping.get(edge_key_reverse)
        disease_db_id = diseases.get(de['disease_id'], {}).get('id')

        if edge_id and disease_db_id:
            de_dict = {
                'disease_id': disease_db_id,
                'edge_id': edge_id,
                'reference': de['reference'],
                'source_db': de['source_db']
            }
            parser.insert_disease_edge(de_dict)

    parser.save_db_to_file(str(OUTPUT_DB))
    print(f"Merge complete: {OUTPUT_DB}")
    print(f"Merged {len(experiment_models)} experiment models")
    print(f"Merged {len(disease_edges)} disease edges")


if __name__ == "__main__":
    main()
