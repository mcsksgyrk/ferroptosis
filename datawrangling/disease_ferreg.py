from database.sqlite_db_api3 import PsimiSQL
from config import PROJECT_ROOT, OUTPUTS_DIR
from parsers.ferreg_parser import FerregParser


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
    print("Inserting diseases")
    for disease_id, disease_dict in parser.diseases.items():
        db_api.insert_disease(disease_dict)
    print("Inserting eddges and experiments")
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
            print(edge_id)
            unique_id = edge_dict['_unique_id_']
            experiment_dict = parser.experiments[unique_id].copy()
            experiment_dict['edge_id'] = edge_id
            db_api.insert_experiment_model(experiment_dict)
        else:
            print(f"Warning: Could not find nodes for edge {edge_dict['interactor_a_node_name']} -> {edge_dict['interactor_b_node_name']}")

    db_api.save_db_to_file(str(output_path))


db_path = OUTPUTS_DIR / "ferreg.db"
parser = FerregParser(db_path)
parser.parse_interactions()

db_final = OUTPUTS_DIR / "tttttt.db"
save_to_database(db_final, parser)
