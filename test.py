from apicalls.mygene import MyGeneClient
from database.external_db import DBconnector
from config import OUTPUTS_DIR


def is_biomolecule_mygene(entity_name, mygene_client):
    non_biomolecules = {'FERROPTOSIS', 'APOPTOSIS', 'LIPID ROS', 'ROS', 'OXIDATIVE STRESS'}
    if entity_name.upper() in non_biomolecules:
        return False, "biological_process"
    try:
        result = mygene_client.query_gene(entity_name)
        if result.get('hits') and len(result['hits']) > 0:
            return True, "gene"
        return False, "unknown"
    except Exception as e:
        print(f"MyGene lookup failed for {entity_name}: {e}")
        return False, "unknown"


def trace_ferroptosis_effects(pathway_str, mygene_client):
    if not pathway_str:
        return [], ""
    steps = [step.strip() for step in pathway_str.split(',')]
    edges_with_effects = []
    direct_effect = ""
    ferroptosis_map = {}
    for step in steps:
        if 'ferroptosis' in step.lower():
            if ':-:' in step:
                source = step.split(':-:')[0].strip()
                ferroptosis_map[source] = "suppresses ferroptosis"
            elif ':+:' in step:
                source = step.split(':+:')[0].strip()
                ferroptosis_map[source] = "promotes ferroptosis"
    for step in steps:
        if 'ferroptosis' not in step.lower():
            if ':-:' in step or ':+:' in step:
                source, target = step.replace(':-:', '|').replace(':+:', '|').split('|')
                source, target = source.strip(), target.strip()
                source_is_bio, source_type = is_biomolecule_mygene(source, mygene_client)
                target_is_bio, target_type = is_biomolecule_mygene(target, mygene_client)
                if source_is_bio and target_is_bio:
                    effect = ferroptosis_map.get(target, "")
                    interaction_type = 'inhibition' if ':-:' in step else 'activation'
                    edges_with_effects.append((source, target, interaction_type, effect))
                else:
                    print(f"Skipping {source} -> {target}: source_is_bio={source_is_bio}, target_is_bio={target_is_bio}")
    return edges_with_effects, direct_effect


mygene = MyGeneClient()
ferrdb_path = OUTPUTS_DIR / "ferrdb.db"
db = DBconnector(ferrdb_path)
query = "SELECT * FROM suppressor"
suppressor = db.query_to_dataframe(query)

results = []
for idx, row in suppressor.iterrows():
    result = trace_ferroptosis_effects(row['Pathway'], mygene)
    results.append(results)

suppressor['edges'] = results

