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


def extend_merged_db_with_omnipath():
    merged_db_path = OUTPUTS_DIR / "merged_ferroptosis_network.db"
    output_db_path = OUTPUTS_DIR / "merged_ferroptosis_w_omnipath.db"
    omnipath_file = SOURCES_DIR / "omnipath" / "omnipath_interactions.txt"

    if not merged_db_path.exists():
        raise FileNotFoundError(f"Merged database not found: {merged_db_path}")
    if not omnipath_file.exists():
        raise FileNotFoundError(f"OmniPath file not found: {omnipath_file}")

    logger.info(f"Loading existing database: {merged_db_path}")
    sql_seed = PROJECT_ROOT / "database" / "network_db_seed3.sql"
    db_api = PsimiSQL(sql_seed)
    db_api.import_from_db_file(str(merged_db_path))

    logger.info(f"Loading OmniPath interactions: {omnipath_file}")
    omnipath_df = pd.read_csv(omnipath_file, delimiter='\t')
    logger.info(f"Found {len(omnipath_df)} OmniPath interactions")

    logger.info("Building protein sets from existing database...")
    kegg_proteins = set()
    all_existing_proteins = set()

    db_api.cursor.execute("SELECT id, name, source_db FROM node WHERE type = 'protein'")
    for node_id, name, source in db_api.cursor.fetchall():
        is_kegg = ('KEGG' in source)

        db_api.cursor.execute("SELECT id_value FROM node_identifier WHERE node_id = ? AND id_type = 'uniprot_id'", (node_id,))
        uniprot_ids = [row[0] for row in db_api.cursor.fetchall()]

        if uniprot_ids:
            for uniprot_id in uniprot_ids:
                all_existing_proteins.add(uniprot_id)
                if is_kegg:
                    kegg_proteins.add(uniprot_id)
        else:
            logger.warning(f"Node '{name}' (source: {source}) has no UniProt ID")
            all_existing_proteins.add(name)
            if is_kegg:
                kegg_proteins.add(name)

    logger.info(f"Found {len(kegg_proteins)} KEGG proteins, {len(all_existing_proteins)} total existing proteins")

    logger.info("Building network layers...")
    layer1_proteins = set()
    layer2_proteins = set()
    new_nodes_added = set()

    # only considers directions pointing towards KEGG proteins
    #for _, row in omnipath_df.iterrows():
    #    source_id, target_id = row['source'], row['target']
    #    if target_id in kegg_proteins and source_id not in kegg_proteins:
    #        layer1_proteins.add(source_id)
    #        if source_id not in all_existing_proteins:
    #            new_nodes_added.add(source_id)

    #for _, row in omnipath_df.iterrows():
    #    source_id, target_id = row['source'], row['target']
    #    if target_id in layer1_proteins and source_id not in kegg_proteins:
    #        layer2_proteins.add(source_id)
    #        if source_id not in all_existing_proteins:
    #            new_nodes_added.add(source_id)


    # considers both directions, but same distance.
    for _, row in omnipath_df.iterrows():
        source_id, target_id = row['source'], row['target']
        if target_id in kegg_proteins and source_id not in kegg_proteins:
            layer1_proteins.add(source_id)
            if source_id not in all_existing_proteins:
                new_nodes_added.add(source_id)
        if source_id in kegg_proteins and target_id not in kegg_proteins:
            layer1_proteins.add(target_id)
            if target_id not in all_existing_proteins:
                new_nodes_added.add(target_id)

    for _, row in omnipath_df.iterrows():
        source_id, target_id = row['source'], row['target']
        if target_id in layer1_proteins and source_id not in kegg_proteins:
            layer2_proteins.add(source_id)
            if source_id not in all_existing_proteins:
                new_nodes_added.add(source_id)
        if source_id in layer1_proteins and target_id not in kegg_proteins:
            layer2_proteins.add(target_id)
            if target_id not in all_existing_proteins:
                new_nodes_added.add(target_id)

    logger.info(f"Layer 1: {len(layer1_proteins)} proteins, Layer 2: {len(layer2_proteins)} proteins")
    logger.info(f"New nodes to add: {len(new_nodes_added)}")

    if new_nodes_added:
        logger.info("Adding new nodes...")
        for protein_id in new_nodes_added:
            try:
                node_dict = {
                    'name': protein_id,
                    'primary_id_type': 'uniprot_id',
                    'display_name': protein_id,
                    'tax_id': 9606,
                    'type': 'protein',
                    'pathways': '',
                    'source_db': 'OmniPath',
                    'function': ''
                }
                db_api.insert_node(node_dict)
            except Exception as e:
                logger.warning(f"Failed to add node {protein_id}: {e}")
        logger.info(f"Added {len(new_nodes_added)} new nodes")

    logger.info("Processing edges...")
    edges_added = 0
    edges_updated = 0
    edges_skipped = 0
    layer_counts = {0: 0, 1: 0, 2: 0}

    for _, interaction in omnipath_df.iterrows():
        try:
            source_dict = db_api.get_node_by_any_identifier(interaction['source'])
            target_dict = db_api.get_node_by_any_identifier(interaction['target'])

            if not source_dict or not target_dict:
                edges_skipped += 1
                continue

            source_id, target_id = interaction['source'], interaction['target']

            if source_id in kegg_proteins and target_id in kegg_proteins:
                layer = 0
            elif (source_id in kegg_proteins and target_id in layer1_proteins) or \
                 (target_id in kegg_proteins and source_id in layer1_proteins):
                layer = 1
            elif source_id in layer1_proteins and target_id in layer1_proteins:
                layer = 2
            elif (source_id in layer1_proteins and target_id in layer2_proteins) or \
                 (target_id in layer1_proteins and source_id in layer2_proteins):
                layer = 2
            else:
                edges_skipped += 1
                continue

            layer_counts[layer] += 1

            interaction_types = []
            if 'is_directed' in interaction and pd.notna(interaction['is_directed']):
                is_directed = 'true' if interaction['is_directed'] == 1 else 'false'
                interaction_types.extend([f"is_directed:{is_directed}", f"is_direct:{is_directed}"])
            else:
                interaction_types.extend(["is_directed:false", "is_direct:false"])

            if 'sources' in interaction and pd.notna(interaction['sources']):
                interaction_types.append(f"sources:{interaction['sources']}")

            # Check for existing edge
            existing_query = """
                SELECT id, layer, source_db, interaction_types FROM edge
                WHERE interactor_a_node_id = ? AND interactor_b_node_id = ?
            """
            db_api.cursor.execute(existing_query, (source_dict['id'], target_dict['id']))
            existing_edge = db_api.cursor.fetchone()

            if existing_edge:
                # Update existing edge
                edge_id, old_layer, old_source_db, old_interaction_types = existing_edge

                # Merge source databases
                new_source_db = merge_strings(old_source_db, 'OmniPath')

                # Merge interaction types
                new_interaction_types = merge_strings(old_interaction_types, '|'.join(interaction_types))

                # Update with new layer
                update_query = """
                    UPDATE edge
                    SET layer = ?, source_db = ?, interaction_types = ?
                    WHERE id = ?
                """
                db_api.cursor.execute(update_query, (str(layer), new_source_db, new_interaction_types, edge_id))
                db_api.db.commit()
                edges_updated += 1
            else:
                # Insert new edge
                edge_dict = {
                    'source_db': 'OmniPath',
                    'interaction_types': '|'.join(interaction_types),
                    'layer': str(layer)
                }
                db_api.insert_edge(source_dict, target_dict, edge_dict)
                edges_added += 1

        except Exception as e:
            logger.warning(f"Failed to process edge {interaction['source']} -> {interaction['target']}: {e}")
            edges_skipped += 1

    logger.info(f"Added {edges_added} new edges, updated {edges_updated} existing edges, skipped {edges_skipped}")
    logger.info(f"Layer distribution: Layer 0: {layer_counts[0]}, Layer 1: {layer_counts[1]}, Layer 2: {layer_counts[2]}")

    logger.info(f"Saving extended database: {output_db_path}")
    db_api.save_db_to_file(str(output_db_path))

    logger.info("OmniPath integration complete!")
