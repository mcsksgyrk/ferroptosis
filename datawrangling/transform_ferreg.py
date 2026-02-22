from database.sqlite_db_api3 import PsimiSQL
from config import PROJECT_ROOT, OUTPUTS_DIR
from parsers.ferreg_parser import FerregParser
from collections import defaultdict


def get_nodes_field(keys, parser, field='name'):
    true_names = []
    for key in keys:
        true_names.append(parser.nodes[key]['name'])
    return true_names


def determine_id_type(value):
    value_str = str(value).upper()
    if (len(value_str) == 6 and value_str[0] in ['O', 'P', 'Q', 'A', 'B']):
        return 'uniprot_id'
    elif value_str.startswith('ENSG'):
        return 'ensembl_id'
    elif value_str.startswith('MIMAT'):
        return 'mirbase_id'
    elif len(value_str) == 27:
        return 'inchikey'
    else:
        return 'external_id'


def save_to_database(output_path, parser):
    SQL_SEED = PROJECT_ROOT / "database" / "network_db_seed3.sql"
    db_api = PsimiSQL(SQL_SEED)

    unique_id_disease_map = {}
    for disease_key, disease_dict in parser.diseases.items():
        for unique_id in disease_dict['unique_ids']:
            if unique_id not in unique_id_disease_map:
                unique_id_disease_map[unique_id] = []
            unique_id_disease_map[unique_id].append(disease_key)

    # Check for unique_ids with multiple diseases
    multi_disease_uids = {uid: diseases for uid, diseases in unique_id_disease_map.items()
                          if len(diseases) > 1}
    if multi_disease_uids:
        print(f"Warning: {len(multi_disease_uids)} unique_ids associated with multiple diseases")
        for uid, diseases in list(multi_disease_uids.items())[:5]:
            print(f"  unique_id {uid}: {diseases}")

    print("Inserting nodes")
    for internal_id, node_dict in parser.nodes.items():
        db_api.insert_node(node_dict)

    print("Inserting identifiers")
    for internal_id, node_dict in parser.nodes.items():
        id_value = node_dict['name']
        id_type = determine_id_type(id_value)
        is_primary = True
        node_id = db_api.get_node_by_any_identifier(id_value)['id']
        db_api.insert_node_identifier(node_id, id_type, id_value, is_primary)

    # MODIFICATION: Insert all diseases (including duplicates for now to debug)
    print("Inserting diseases")
    disease_mapping = {}
    disease_dedup = {}  # Track unique diseases by their content

    for disease_key, disease_dict in parser.diseases.items():
        disease_name = disease_dict.get('disease_name', '')
        disease_id = disease_dict.get('disease_id', '')

        # Skip completely invalid diseases
        if (disease_id in ['ICD-11: N.A.', 'N.A.', 'NA', ''] and
            disease_name in ['N.A.', 'NA', '', 'Not Available']):
            continue

        # Create a unique key for deduplication
        dedup_key = f"{disease_id}|{disease_name}"

        if dedup_key in disease_dedup:
            # Disease already exists, map to existing ID
            disease_mapping[disease_key] = disease_dedup[dedup_key]
        else:
            # Insert new disease
            db_api.insert_disease(disease_dict)
            db_api.cursor.execute("SELECT last_insert_rowid()")
            db_disease_id = db_api.cursor.fetchone()[0]
            disease_mapping[disease_key] = db_disease_id
            disease_dedup[dedup_key] = db_disease_id

    print(f"Inserted {len(disease_dedup)} unique diseases")

    print("Inserting edges and experiments")
    unique_id_to_edge_ids = defaultdict(list)

    for edge_dict in parser.edges:
        a_ferreg_id = edge_dict['interactor_a_node_name']
        b_ferreg_id = edge_dict['interactor_b_node_name']
        a_name, b_name = get_nodes_field([a_ferreg_id, b_ferreg_id], parser)

        source_node = db_api.get_node_by_any_identifier(a_name)
        target_node = db_api.get_node_by_any_identifier(b_name)

        if source_node and target_node:
            db_edge_dict = {
                'layer': edge_dict['layer'],
                'interaction_types': edge_dict['interaction_types'],
                'effect_on_ferroptosis': edge_dict['effect_on_ferroptosis'],
                'source_db': edge_dict['source_db']
            }
            db_api.insert_edge(source_node, target_node, db_edge_dict)
            db_api.cursor.execute("SELECT last_insert_rowid()")
            edge_id = db_api.cursor.fetchone()[0]
            unique_id = edge_dict['_unique_id_']

            # Insert experiment
            if unique_id in parser.experiments:
                experiment_dict = parser.experiments[unique_id].copy()
                experiment_dict['edge_id'] = edge_id
                db_api.insert_experiment_model(experiment_dict)

            # Store mapping
            unique_id_to_edge_ids[unique_id].append(edge_id)

    # MODIFICATION: Correct disease-edge association
    print("Inserting disease-edge associations")
    disease_edge_count = 0
    edge_disease_pairs = set()  # Track unique pairs to avoid duplicates

    for unique_id, edge_ids in unique_id_to_edge_ids.items():
        # Find the disease for this unique_id (should be only one)
        disease_found = False

        for disease_key, disease_dict in parser.diseases.items():
            if unique_id in disease_dict['unique_ids']:
                if disease_key not in disease_mapping:
                    continue  # Skip invalid diseases

                db_disease_id = disease_mapping[disease_key]

                # Create entries for all edges with this unique_id
                for edge_id in edge_ids:
                    pair = (db_disease_id, edge_id)

                    # Only insert if this pair hasn't been inserted yet
                    if pair not in edge_disease_pairs:
                        disease_edge_dict = {
                            'disease_id': db_disease_id,
                            'edge_id': edge_id,
                            'reference': '',
                            'source_db': "FerReg"
                        }
                        db_api.insert_disease_edge(disease_edge_dict)
                        edge_disease_pairs.add(pair)
                        disease_edge_count += 1

                disease_found = True
                break  # Each unique_id should map to only one disease

        if not disease_found and unique_id_to_edge_ids[unique_id]:
            print(f"No disease found for unique_id {unique_id} with {len(edge_ids)} edges")

    print(f"Created {disease_edge_count} disease-edge associations")
    print(f"Unique disease-edge pairs: {len(edge_disease_pairs)}")

    db_api.save_db_to_file(str(output_path))


def convert_ferreg_source():
    db_path = OUTPUTS_DIR / "ferreg.db"
    parser = FerregParser(db_path)
    parser.parse_interactions()

    db_final = OUTPUTS_DIR / "ferreg_network.db"
    save_to_database(db_final, parser)
