import pandas as pd
from pathlib import Path
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from database.sqlite_db_api3 import PsimiSQL
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


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


def extend_arn_ferr_db_with_omnipath():
    arn_ferr_db_path = OUTPUTS_DIR / "ferroptosis_autophagy.db"
    output_db_path = OUTPUTS_DIR / "final.db"
    omnipath_file = SOURCES_DIR / "omnipath" / "omnipath_interactions.txt"

    if not arn_ferr_db_path.exists():
        raise FileNotFoundError(f"ARN-Ferroptosis database not found: {arn_ferr_db_path}")
    if not omnipath_file.exists():
        raise FileNotFoundError(f"OmniPath file not found: {omnipath_file}")

    logger.info(f"Loading existing database: {arn_ferr_db_path}")
    sql_seed = PROJECT_ROOT / "database" / "network_db_seed3.sql"
    db_api = PsimiSQL(sql_seed)
    db_api.import_from_db_file(str(arn_ferr_db_path))

    logger.info(f"Loading OmniPath interactions: {omnipath_file}")
    omnipath_df = pd.read_csv(omnipath_file, delimiter='\t')
    logger.info(f"Found {len(omnipath_df)} OmniPath interactions")

    logger.info("Step 1: Get KEGG core proteins (layer 0)")
    kegg_proteins = set()
    db_api.cursor.execute("SELECT name FROM node WHERE source_db LIKE '%KEGG%' AND primary_id_type = 'uniprot_id'")
    for name, in db_api.cursor.fetchall():
        kegg_proteins.add(name)
    logger.info(f"KEGG core proteins: {len(kegg_proteins)}")

    logger.info("Step 2: Get all ferroptosis proteins")
    ferroptosis_proteins = set()
    db_api.cursor.execute("SELECT name FROM node WHERE source_db LIKE '%ferr%' AND primary_id_type = 'uniprot_id'")
    for name, in db_api.cursor.fetchall():
        ferroptosis_proteins.add(name)
    logger.info(f"Total ferroptosis proteins: {len(ferroptosis_proteins)}")

    logger.info("Step 3: Get ARN-only proteins")
    arn_proteins = set()
    db_api.cursor.execute("SELECT name FROM node WHERE source_db = 'ARN'")
    for name, in db_api.cursor.fetchall():
        arn_proteins.add(name)
    logger.info(f"ARN-only proteins: {len(arn_proteins)}")

    logger.info("Step 4: Build interaction map")
    interactions = {}
    db_api.cursor.execute("SELECT interactor_a_node_name, interactor_b_node_name FROM edge")
    for a, b in db_api.cursor.fetchall():
        if a not in interactions:
            interactions[a] = set()
        if b not in interactions:
            interactions[b] = set()
        interactions[a].add(b)
        interactions[b].add(a)

    logger.info("Step 5: Classify layer 1 proteins")
    layer1_proteins = set()
    for protein in ferroptosis_proteins:
        if protein in kegg_proteins:
            continue
        if protein in interactions:
            for neighbor in interactions[protein]:
                if neighbor in kegg_proteins:
                    layer1_proteins.add(protein)
                    break
    logger.info(f"Layer 1 proteins: {len(layer1_proteins)}")

    logger.info("Step 6: Classify layer 2 proteins")
    layer2_proteins = set()
    for protein in ferroptosis_proteins:
        if protein in kegg_proteins or protein in layer1_proteins:
            continue
        if protein in interactions:
            for neighbor in interactions[protein]:
                if neighbor in layer1_proteins:
                    layer2_proteins.add(protein)
                    break
    logger.info(f"Layer 2 proteins: {len(layer2_proteins)}")

    logger.info("Step 7: Remaining are ferreg proteins")
    ferreg_proteins = ferroptosis_proteins - kegg_proteins - layer1_proteins - layer2_proteins
    logger.info(f"FerReg proteins: {len(ferreg_proteins)}")

    logger.info("Step 8: Build existing edge lookup")
    existing_edges = set()
    db_api.cursor.execute("SELECT interactor_a_node_name, interactor_b_node_name FROM edge")
    for a, b in db_api.cursor.fetchall():
        existing_edges.add((a, b))

    logger.info("Step 9: Process cross-network edges")
    edges_added = 0
    edges_skipped = 0
    layer_counts = {0: 0, 1: 0, 2: 0, 'arn_ferr_cross': 0}

    for _, interaction in omnipath_df.iterrows():
        try:
            source_id, target_id = interaction['source'], interaction['target']

            source_is_ferro = source_id in ferroptosis_proteins
            target_is_ferro = target_id in ferroptosis_proteins
            source_is_arn = source_id in arn_proteins
            target_is_arn = target_id in arn_proteins

            if not ((source_is_ferro and target_is_arn) or (source_is_arn and target_is_ferro)):
                continue

            if (source_id, target_id) in existing_edges:
                continue

            source_dict = db_api.get_node_by_any_identifier(source_id)
            target_dict = db_api.get_node_by_any_identifier(target_id)

            if not source_dict or not target_dict:
                edges_skipped += 1
                continue

            ferro_node = source_id if source_is_ferro else target_id

            if ferro_node in kegg_proteins:
                layer = 0
            elif ferro_node in layer1_proteins:
                layer = 1
            elif ferro_node in layer2_proteins:
                layer = 2
            elif ferro_node in ferreg_proteins:
                layer = 'arn_ferr_cross'
            else:
                layer = 'arn_ferr_cross'

            layer_counts[layer] += 1

            interaction_types = []
            if 'is_directed' in interaction and pd.notna(interaction['is_directed']):
                is_directed = 'true' if interaction['is_directed'] == 1 else 'false'
                interaction_types.extend([f"is_directed:{is_directed}", f"is_direct:{is_directed}"])
            else:
                interaction_types.extend(["is_directed:false", "is_direct:false"])

            if 'sources' in interaction and pd.notna(interaction['sources']):
                interaction_types.append(f"sources:{interaction['sources']}")

            edge_dict = {
                'source_db': 'Cross-Network',
                'interaction_types': '|'.join(interaction_types),
                'layer': str(layer),
                'effect_on_ferroptosis': ''
            }
            db_api.insert_edge(source_dict, target_dict, edge_dict)
            edges_added += 1

        except Exception as e:
            logger.warning(f"Failed to process edge {interaction['source']} -> {interaction['target']}: {e}")
            edges_skipped += 1

    logger.info(f"Added {edges_added} new cross-network edges, skipped {edges_skipped}")
    logger.info(f"Layer distribution: Layer 0: {layer_counts[0]}, Layer 1: {layer_counts[1]}, Layer 2: {layer_counts[2]}, Cross: {layer_counts['arn_ferr_cross']}")

    logger.info(f"Saving final database: {output_db_path}")
    db_api.save_db_to_file(str(output_db_path))

    logger.info("Cross-network integration complete!")


if __name__ == "__main__":
    extend_arn_ferr_db_with_omnipath()
