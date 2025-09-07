
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

    # Collect all nodes and identifiers
    all_nodes = []
    all_identifiers = {}
    diseases = {}
    collected_edges_directed = {}
    collected_edges_undirected = {}

    for db_file in input_dbs:
        if not os.path.exists(db_file):
            continue

        print(f"Reading: {db_file}")
        db = sqlite3.connect(db_file)
        cursor = db.cursor()

        # Get all columns for node table
        cursor.execute("PRAGMA table_info(node)")
        node_columns = [col[1] for col in cursor.fetchall()]

        # Read nodes
        cursor.execute("SELECT * FROM node")
        db_nodes = cursor.fetchall()

        # Read node identifiers if table exists
        node_identifiers = {}
        try:
            cursor.execute("SELECT * FROM node_identifier")
            for row in cursor.fetchall():
                node_db_id, id_type, is_primary, id_value = row
                if node_db_id not in node_identifiers:
                    node_identifiers[node_db_id] = {}
                node_identifiers[node_db_id][id_type] = {
                    'value': id_value,
                    'is_primary': is_primary
                }
        except:
            pass

        # Process each node
        for node_row in db_nodes:
            # Map row to columns
            node_dict = {}
            for i, col in enumerate(node_columns):
                node_dict[col] = node_row[i]

            node_data = {
                "name": node_dict['name'],
                "primary_id_type": node_dict['primary_id_type'],
                "display_name": node_dict['display_name'],
                "tax_id": node_dict['tax_id'],
                "type": node_dict['type'],
                "pathways": node_dict['pathways'] or "",
                "role_in_ferroptosis": node_dict.get('role_in_ferroptosis', '') or "",
                "function": node_dict['function'] or "",
                "source_db": node_dict.get('source_db', '') or "",
                "identifiers": node_identifiers.get(node_dict['id'], {})
            }

            node_index = len(all_nodes)
            all_nodes.append(node_data)
            if str(db_file).endswith('kegg.db'):
                print(f"KEGG node collected: {node_data['name']} - {node_data['type']}")

            # Map all identifiers
            all_ids = [node_data["name"]]
            for id_info in node_data["identifiers"].values():
                if id_info["value"]:
                    all_ids.append(id_info["value"])

            for identifier in all_ids:
                if identifier not in all_identifiers:
                    all_identifiers[identifier] = []
                all_identifiers[identifier].append(node_index)

        # Read diseases if table exists
        try:
            cursor.execute("SELECT * FROM disease")
            for row in cursor.fetchall():
                disease_id, disease_name, description = row[1], row[2], row[3]
                if disease_id not in diseases:
                    diseases[disease_id] = {
                        'disease_id': disease_id,
                        'disease_name': disease_name or "",
                        'description': description or ""
                    }
                else:
                    diseases[disease_id]['description'] = merge_strings(
                        diseases[disease_id]['description'],
                        description or ""
                    )
        except:
            pass

        db.close()

    print(f"Total nodes collected: {len(all_nodes)}")

    # Find connected components
    visited = set()
    node_groups = []

    def dfs(node_idx, group):
        if node_idx in visited:
            return
        visited.add(node_idx)
        group.append(node_idx)

        node = all_nodes[node_idx]
        all_ids = [node["name"]]
        for id_info in node["identifiers"].values():
            if id_info["value"]:
                all_ids.append(id_info["value"])

        for identifier in all_ids:
            if identifier in all_identifiers:
                for related_idx in all_identifiers[identifier]:
                    if related_idx not in visited:
                        dfs(related_idx, group)

    for i in range(len(all_nodes)):
        if i not in visited:
            group = []
            dfs(i, group)
            node_groups.append(group)

    print(f"Unique nodes after deduplication: {len(node_groups)}")
    for group_idx, group in enumerate(node_groups):
        for node_idx in group:
            node = all_nodes[node_idx]
            if node["name"] == "Q13772":
                print(f"Q13772 found in group {group_idx} with {len(group)} nodes:")
                print(f"  Canonical name will be: {all_nodes[group[0]]['name']}")
                for idx in group:
                    print(f"  - {all_nodes[idx]['name']} from {all_nodes[idx]['db_file']}")
                break
    # Merge nodes
    nodes = {}
    node_identifiers = {}
    name_mapping = {}

    for group in node_groups:
        base_node = all_nodes[group[0]]
        merged_node = dict(base_node)
        merged_identifiers = dict(base_node["identifiers"])

        for idx in group[1:]:
            other_node = all_nodes[idx]
            merged_node = get_union_of_nodes(merged_node, other_node)

            for id_type, id_info in other_node["identifiers"].items():
                if id_type not in merged_identifiers:
                    merged_identifiers[id_type] = id_info

        canonical_name = base_node["name"]
        nodes[canonical_name] = merged_node
        node_identifiers[canonical_name] = merged_identifiers

        # Map all names to canonical
        for idx in group:
            node = all_nodes[idx]
            name_mapping[node["name"]] = canonical_name
            for id_info in node["identifiers"].values():
                if id_info["value"]:
                    name_mapping[id_info["value"]] = canonical_name

    # Create database
    parser = PsimiSQL(SQL_SEED)

    for node in nodes.values():
        parser.insert_unique_node(node)
        nodes[node['name']]['id'] = parser.cursor.lastrowid

    for node_name, identifiers in node_identifiers.items():
        if node_name in nodes:
            node_id = nodes[node_name]['id']
            for id_type, id_data in identifiers.items():
                parser.insert_node_identifier(node_id, id_type, id_data['value'], id_data['is_primary'])

    for disease in diseases.values():
        parser.insert_disease(disease)
        diseases[disease['disease_id']]['id'] = parser.cursor.lastrowid

    print(f"Inserted {len(nodes)} nodes and {len(diseases)} diseases")

    # Process edges
    all_edge_counter = 0
    merged_edge_counter = 0

    for db_file in input_dbs:
        if not os.path.exists(db_file):
            continue

        print(f"Processing edges from: {db_file}")
        db = sqlite3.connect(db_file)
        cursor = db.cursor()

        # Get edge table columns
        cursor.execute("PRAGMA table_info(edge)")
        edge_columns = [col[1] for col in cursor.fetchall()]

        cursor.execute("SELECT * FROM edge")

        for row in cursor.fetchall():
            all_edge_counter += 1

            # Map row to columns
            edge_dict = {}
            for i, col in enumerate(edge_columns):
                edge_dict[col] = row[i]

            source_name = name_mapping.get(edge_dict['interactor_a_node_name'])
            target_name = name_mapping.get(edge_dict['interactor_b_node_name'])

            if not source_name or not target_name:
                continue

            edge_key = f"{source_name}@{target_name}@{edge_dict['layer']}"
            edge_key_reverse = f"{target_name}@{source_name}@{edge_dict['layer']}"

            current_edge = {
                'interaction_types': edge_dict['interaction_types'] or "",
                'effect_on_ferroptosis': edge_dict.get('effect_on_ferroptosis', '') or "",
                'source_db': edge_dict['source_db'] or "",
                'layer': edge_dict['layer']
            }

            directed = current_edge['interaction_types'] and 'is_directed:true' in current_edge['interaction_types'].split('|')

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
    for edge_key, edge_data in {**collected_edges_undirected, **collected_edges_directed}.items():
        node_a, node_b, layer = edge_key.split('@')
        parser.insert_edge(nodes[node_a], nodes[node_b], edge_data)

    parser.save_db_to_file(str(OUTPUT_DB))
    print(f"Merge complete: {OUTPUT_DB}")


if __name__ == "__main__":
    main()
