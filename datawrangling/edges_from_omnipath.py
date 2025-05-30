import pandas as pd
from pathlib import Path
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from database.sqlite_db_api2 import PsimiSQL
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def extend_merged_db_with_omnipath():
    """Add OmniPath interactions to the existing merged database with layered expansion."""

    # Paths
    merged_db_path = OUTPUTS_DIR / "merged_network.db"
    test_db_path = OUTPUTS_DIR / "test_omnipath.db"
    omnipath_file = SOURCES_DIR / "omnipath" / "omnipath_interactions.txt"

    if not merged_db_path.exists():
        raise FileNotFoundError(f"Merged database not found: {merged_db_path}")
    if not omnipath_file.exists():
        raise FileNotFoundError(f"OmniPath file not found: {omnipath_file}")

    # Load existing database
    logger.info(f"Loading existing database: {merged_db_path}")
    sql_seed = PROJECT_ROOT / "database" / "network_db_seed2.sql"
    db_api = PsimiSQL(sql_seed)
    db_api.import_from_db_file(str(merged_db_path))

    # Load OmniPath data
    logger.info(f"Loading OmniPath interactions: {omnipath_file}")
    omnipath_df = pd.read_csv(omnipath_file, delimiter='\t')
    logger.info(f"Found {len(omnipath_df)} OmniPath interactions")

    # Build existing protein sets (UniProt IDs only, since OmniPath uses UniProt)
    logger.info("Building protein sets from existing database...")
    kegg_proteins = set()
    all_existing_proteins = set()

    db_api.cursor.execute("SELECT id, name, source FROM node WHERE type = 'protein'")
    for node_id, name, source in db_api.cursor.fetchall():
        is_kegg = (source == 'KEGG')

        # Get UniProt IDs for this node
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

    # Layered network expansion
    logger.info("Building network layers...")
    layer1_proteins = set()
    layer2_proteins = set()
    new_nodes_added = set()

    # Find Layer 1: proteins that influence KEGG
    for _, row in omnipath_df.iterrows():
        source_id, target_id = row['source'], row['target']
        if target_id in kegg_proteins and source_id not in kegg_proteins:
            layer1_proteins.add(source_id)
            if source_id not in all_existing_proteins:
                new_nodes_added.add(source_id)

    # Find Layer 2: proteins that influence Layer 1
    for _, row in omnipath_df.iterrows():
        source_id, target_id = row['source'], row['target']
        if target_id in layer1_proteins and source_id not in kegg_proteins:
            layer2_proteins.add(source_id)
            if source_id not in all_existing_proteins:
                new_nodes_added.add(source_id)

    logger.info(f"Layer 1: {len(layer1_proteins)} proteins, Layer 2: {len(layer2_proteins)} proteins")
    logger.info(f"New nodes to add: {len(new_nodes_added)}")

    # Add new nodes to database
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
                    'source': 'OmniPath',
                    'function': ''
                }
                db_api.insert_node(node_dict)
            except Exception as e:
                logger.warning(f"Failed to add node {protein_id}: {e}")
        logger.info(f"Added {len(new_nodes_added)} new nodes")

    # Process edges with layer assignment
    logger.info("Processing edges...")
    edges_added = 0
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

            # Determine layer
            if source_id in kegg_proteins and target_id in kegg_proteins:
                layer = 0  # KEGG ↔ KEGG
            elif (source_id in kegg_proteins and target_id in layer1_proteins) or \
                 (target_id in kegg_proteins and source_id in layer1_proteins):
                layer = 1  # KEGG ↔ Layer 1
            elif source_id in layer1_proteins and target_id in layer1_proteins:
                layer = 2  # Layer 1 ↔ Layer 1
            elif (source_id in layer1_proteins and target_id in layer2_proteins) or \
                 (target_id in layer1_proteins and source_id in layer2_proteins):
                layer = 2  # Layer 1 ↔ Layer 2
            else:
                edges_skipped += 1
                continue

            layer_counts[layer] += 1

            # Build interaction types
            interaction_types = []
            if 'is_directed' in interaction and pd.notna(interaction['is_directed']):
                is_directed = 'true' if interaction['is_directed'] == 1 else 'false'
                interaction_types.extend([f"is_directed:{is_directed}", f"is_direct:{is_directed}"])
            else:
                interaction_types.extend(["is_directed:false", "is_direct:false"])

            if 'sources' in interaction and pd.notna(interaction['sources']):
                interaction_types.append(f"sources:{interaction['sources']}")

            # Check if edge already exists
            existing_query = """
                SELECT id FROM edge
                WHERE interactor_a_node_id = ? AND interactor_b_node_id = ?
            """
            db_api.cursor.execute(existing_query, (source_dict['id'], target_dict['id']))
            if db_api.cursor.fetchone():
                edges_skipped += 1
                continue

            # Create edge
            edge_dict = {
                'source_db': 'OmniPath',
                'interaction_types': '|'.join(interaction_types),
                'layer': str(layer)
            }

            db_api.insert_edge(source_dict, target_dict, edge_dict)
            edges_added += 1

        except Exception as e:
            logger.warning(f"Failed to add edge {interaction['source']} -> {interaction['target']}: {e}")
            edges_skipped += 1

    logger.info(f"Added {edges_added} edges, skipped {edges_skipped}")
    logger.info(f"Layer distribution: Layer 0: {layer_counts[0]}, Layer 1: {layer_counts[1]}, Layer 2: {layer_counts[2]}")

    # Save database
    backup_path = OUTPUTS_DIR / "merged_network_backup.db"
    logger.info(f"Creating backup: {backup_path}")
    merged_db_path.rename(backup_path)

    logger.info(f"Saving updated database: {test_db_path}")
    db_api.save_db_to_file(str(test_db_path))

    logger.info("OmniPath integration complete!")


if __name__ == "__main__":
    extend_merged_db_with_omnipath()
