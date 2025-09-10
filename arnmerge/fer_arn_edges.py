import pandas as pd
from pathlib import Path
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from database.sqlite_db_api3 import PsimiSQL
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def extend_arn_ferr_with_cross_edges():
    arn_ferr_db_path = OUTPUTS_DIR / "ferroptosis_autophagy.db"
    output_db_path = OUTPUTS_DIR / "final.db"
    omnipath_file = SOURCES_DIR / "omnipath" / "omnipath_interactions.txt"

    logger.info(f"Loading existing database: {arn_ferr_db_path}")
    sql_seed = PROJECT_ROOT / "database" / "network_db_seed3.sql"
    db_api = PsimiSQL(sql_seed)
    db_api.import_from_db_file(str(arn_ferr_db_path))

    logger.info(f"Loading OmniPath interactions: {omnipath_file}")
    omnipath_df = pd.read_csv(omnipath_file, delimiter='\t')
    omnipath_df = omnipath_df[omnipath_df['is_directed'] == 1]
    logger.info(f"Found {len(omnipath_df)} directed OmniPath interactions")

    # Update ARN edge layers to distinguish them from ferroptosis layers
    logger.info("Updating ARN edge layers (0->10, 1->11, 2->12)...")
    db_api.cursor.execute("""
        UPDATE edge
        SET layer = CASE
            WHEN layer = '0' THEN '10'
            WHEN layer = '1' THEN '11'
            WHEN layer = '2' THEN '12'
            ELSE layer
        END
        WHERE source_db = 'ARN'
    """)
    db_api.db.commit()

    # Get KEGG proteins (ferroptosis layer 0 core)
    kegg_proteins = set()
    db_api.cursor.execute("""
        SELECT DISTINCT n.name
        FROM node n
        WHERE n.source_db LIKE '%KEGG%'
        AND n.primary_id_type = 'uniprot_id'
    """)
    for name, in db_api.cursor.fetchall():
        kegg_proteins.add(name)
    logger.info(f"KEGG core proteins: {len(kegg_proteins)}")

    # Build ferroptosis layer 0 edges and proteins
    layer0_proteins = set(kegg_proteins)
    db_api.cursor.execute("""
        SELECT DISTINCT e.interactor_b_node_name
        FROM edge e
        JOIN node n ON e.interactor_b_node_name = n.name
        WHERE e.interactor_a_node_name IN ({})
        AND n.primary_id_type = 'uniprot_id'
        AND e.source_db NOT LIKE '%ARN%'
    """.format(','.join(['?' for _ in kegg_proteins])), list(kegg_proteins))
    for name, in db_api.cursor.fetchall():
        layer0_proteins.add(name)

    # Build ferroptosis layer 1 proteins
    layer1_proteins = set()
    if layer0_proteins:
        db_api.cursor.execute("""
            SELECT DISTINCT e.interactor_b_node_name
            FROM edge e
            JOIN node n ON e.interactor_b_node_name = n.name
            WHERE e.interactor_a_node_name IN ({})
            AND e.interactor_b_node_name NOT IN ({})
            AND n.primary_id_type = 'uniprot_id'
            AND e.source_db NOT LIKE '%ARN%'
        """.format(','.join(['?' for _ in layer0_proteins]), ','.join(['?' for _ in layer0_proteins])),
           list(layer0_proteins) + list(layer0_proteins))
        for name, in db_api.cursor.fetchall():
            layer1_proteins.add(name)

    logger.info(f"Ferroptosis layer 0 proteins: {len(layer0_proteins)}")
    logger.info(f"Ferroptosis layer 1 proteins: {len(layer1_proteins)}")

    # Get all ARN proteins
    arn_proteins = set()
    db_api.cursor.execute("""
        SELECT DISTINCT name
        FROM node
        WHERE source_db LIKE '%ARN%'
        AND primary_id_type = 'uniprot_id'
    """)
    for name, in db_api.cursor.fetchall():
        arn_proteins.add(name)
    logger.info(f"ARN proteins: {len(arn_proteins)}")

    # Get all ferroptosis proteins (non-ARN)
    ferroptosis_proteins = set()
    db_api.cursor.execute("""
        SELECT DISTINCT name
        FROM node
        WHERE source_db NOT LIKE '%ARN%'
        AND primary_id_type = 'uniprot_id'
    """)
    for name, in db_api.cursor.fetchall():
        ferroptosis_proteins.add(name)
    logger.info(f"Ferroptosis proteins: {len(ferroptosis_proteins)}")

    # Get existing edges to avoid duplicates
    existing_edges = set()
    db_api.cursor.execute("SELECT interactor_a_node_name, interactor_b_node_name FROM edge")
    for a, b in db_api.cursor.fetchall():
        existing_edges.add((a, b))

    # Process cross-network edges from OmniPath
    logger.info("Processing cross-network edges...")
    cross_edges_added = 0
    layer_counts = {'cross_0': 0, 'cross_1': 0, 'cross_other': 0}

    for _, interaction in omnipath_df.iterrows():
        source_id = interaction['source']
        target_id = interaction['target']

        # Check if this is a cross-network edge
        source_in_ferro = source_id in ferroptosis_proteins
        target_in_ferro = target_id in ferroptosis_proteins
        source_in_arn = source_id in arn_proteins
        target_in_arn = target_id in arn_proteins

        is_cross_network = (source_in_ferro and target_in_arn) or (source_in_arn and target_in_ferro)

        if not is_cross_network:
            continue

        if (source_id, target_id) in existing_edges:
            continue

        source_dict = db_api.get_node_by_any_identifier(source_id)
        target_dict = db_api.get_node_by_any_identifier(target_id)

        if not source_dict or not target_dict:
            continue

        # Determine cross-network layer based on ferroptosis protein involved
        ferro_protein = source_id if source_in_ferro else target_id

        if ferro_protein in layer0_proteins:
            layer = 'cross_0'
            layer_counts['cross_0'] += 1
        elif ferro_protein in layer1_proteins:
            layer = 'cross_1'
            layer_counts['cross_1'] += 1
        else:
            layer = 'cross_other'
            layer_counts['cross_other'] += 1

        # Only add cross_0 and cross_1 edges as requested
        if layer in ['cross_0', 'cross_1']:
            interaction_types = ["is_directed:true", "is_direct:true"]
            if 'sources' in interaction and pd.notna(interaction['sources']):
                interaction_types.append(f"sources:{interaction['sources']}")

            edge_dict = {
                'source_db': 'omnipath_cross',
                'interaction_types': '|'.join(interaction_types),
                'layer': layer,
                'effect_on_ferroptosis': ''
            }

            db_api.insert_edge(source_dict, target_dict, edge_dict)
            cross_edges_added += 1
            existing_edges.add((source_id, target_id))

    logger.info(f"Added {cross_edges_added} cross-network edges")
    logger.info(f"Cross-network distribution: {layer_counts}")

    # Final summary
    logger.info("\nFinal edge summary:")
    db_api.cursor.execute("""
        SELECT source_db, layer, COUNT(*) as count
        FROM edge
        GROUP BY source_db, layer
        ORDER BY source_db, layer
    """)

    for source_db, layer, count in db_api.cursor.fetchall():
        logger.info(f"  {source_db} - Layer {layer}: {count} edges")

    logger.info(f"\nSaving final database: {output_db_path}")
    db_api.save_db_to_file(str(output_db_path))
    logger.info("Cross-network integration complete!")


if __name__ == "__main__":
    extend_arn_ferr_with_cross_edges()
