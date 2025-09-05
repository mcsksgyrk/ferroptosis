import pandas as pd
from database.sqlite_db_api3 import PsimiSQL
from config import PROJECT_ROOT, OUTPUTS_DIR
from parsers.ferreg_parser import FerregParser


def save_to_database(output_path, parser):
    SQL_SEED = PROJECT_ROOT / "database" / "network_db_seed3.sql"
    db_api = PsimiSQL(SQL_SEED)
    print("Inserting nodes...")
    for internal_id, node_dict in parser.nodes.items():
        db_api.insert_node(node_dict)
    print("Inserting diseases...")
    for disease_id, disease_dict in parser.diseases.items():
        db_api.insert_disease(disease_dict)
    print("Inserting edges...")
    for edge_dict in parser.edges:
        source_node = db_api.get_node_by_any_identifier(edge_dict['interactor_a_node_name'])
        target_node = db_api.get_node_by_any_identifier(edge_dict['interactor_b_node_name'])
        if source_node and target_node:
            db_edge_dict = {
                'layer': edge_dict['layer'],
                'interaction_types': edge_dict['interaction_types'],
                'effect_on_ferroptosis': edge_dict['effect_on_ferroptosis'],
                'source_db': edge_dict['source_db']
            }
            db_api.insert_edge(source_node, target_node, db_edge_dict)
        else:
            print(f"Warning: Could not find nodes for edge {edge_dict['interactor_a_node_name']} -> {edge_dict['interactor_b_node_name']}")

    db_api.save_db_to_file(output_path)


db_path = OUTPUTS_DIR / "ferreg.db"
parser = FerregParser(db_path)
parser.parse_interactions()

db_final = OUTPUTS_DIR / "tttttt.db"
save_to_database(db_final, parser)

len(parser.diseases.keys())
len(set(parser.diseases.keys()))
icds = set()
for k, v in parser.diseases.items():
    if v['disease_id'] in icds:
        print(v['disease_id'])
    icds.add(v['disease_id'])
len(set(icds))
