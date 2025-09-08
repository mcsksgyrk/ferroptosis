import sqlite3
import networkx as nx
import pandas as pd
from pathlib import Path
from config import OUTPUTS_DIR, SOURCES_DIR
from collections import defaultdict

def analyze_edge_distances():
    db_path = OUTPUTS_DIR / "extended_omnipath_network.db"
    omnipath_file = SOURCES_DIR / "omnipath" / "omnipath_interactions.txt"

    # Build OmniPath network graph from source file
    print("Building OmniPath network graph from source file...")
    omnipath_df = pd.read_csv(omnipath_file, sep='\t')

    G = nx.Graph()
    for _, row in omnipath_df.iterrows():
        G.add_edge(row['source'], row['target'])

    print(f"OmniPath graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Get FerReg/ferrdb edges from database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    print("\nGetting FerReg/ferrdb edges...")
    query_target = """
    SELECT DISTINCT n1.name, n2.name, e.source_db, e.interaction_types
    FROM edge e
    JOIN node n1 ON e.interactor_a_node_id = n1.id
    JOIN node n2 ON e.interactor_b_node_id = n2.id
    WHERE (e.source_db = 'FerReg' OR e.source_db = 'ferrdb')
    AND n1.primary_id_type = 'uniprot_id'
    AND n2.primary_id_type = 'uniprot_id'
    """

    cursor.execute(query_target)
    target_edges = cursor.fetchall()
    print(f"Found {len(target_edges)} edges from FerReg/ferrdb with uniprot IDs")

    # Analyze distances
    distance_stats = defaultdict(list)
    unreachable_pairs = []
    direct_connections = []
    not_in_omnipath = []

    for source, target, source_db, interaction_type in target_edges:
        # Check if both nodes exist in OmniPath network
        if source not in G and target not in G:
            not_in_omnipath.append((source, target, source_db, 'both'))
            continue
        elif source not in G:
            not_in_omnipath.append((source, target, source_db, 'source'))
            continue
        elif target not in G:
            not_in_omnipath.append((source, target, source_db, 'target'))
            continue

        # Check if edge already exists in OmniPath
        if G.has_edge(source, target):
            direct_connections.append((source, target, source_db))
            distance_stats[source_db].append(1)
        else:
            # Calculate shortest path
            try:
                distance = nx.shortest_path_length(G, source, target)
                distance_stats[source_db].append(distance)
            except nx.NetworkXNoPath:
                unreachable_pairs.append((source, target, source_db))

    # Print results
    print("\n=== DISTANCE ANALYSIS RESULTS ===")

    for db in ['FerReg', 'ferrdb']:
        if db in distance_stats:
            distances = distance_stats[db]
            if distances:
                print(f"\n{db}:")
                print(f"  Total edges analyzed: {len(distances)}")
                print(f"  Average distance: {sum(distances)/len(distances):.2f}")
                print(f"  Min distance: {min(distances)}")
                print(f"  Max distance: {max(distances)}")

                # Distance distribution
                dist_count = defaultdict(int)
                for d in distances:
                    dist_count[d] += 1

                print(f"  Distance distribution:")
                for dist in sorted(dist_count.keys()):
                    print(f"    Distance {dist}: {dist_count[dist]} edges")

    print(f"\nDirect connections (already in OmniPath): {len(direct_connections)}")
    print(f"Unreachable pairs (both nodes in OmniPath but not connected): {len(unreachable_pairs)}")
    print(f"Nodes not in OmniPath: {len(not_in_omnipath)}")

    # Analyze nodes not in OmniPath
    missing_count = defaultdict(int)
    for _, _, db, missing_type in not_in_omnipath:
        missing_count[missing_type] += 1

    print(f"\nNodes not in OmniPath breakdown:")
    for missing_type, count in missing_count.items():
        print(f"  {missing_type} missing: {count}")

    # Show some examples of different distances
    print("\n=== EXAMPLE EDGES BY DISTANCE ===")
    examples_by_distance = defaultdict(list)

    for source, target, source_db, interaction_type in target_edges:
        if source in G and target in G:
            if not G.has_edge(source, target):
                try:
                    distance = nx.shortest_path_length(G, source, target)
                    if len(examples_by_distance[distance]) < 3:
                        examples_by_distance[distance].append((source, target, source_db))
                except:
                    pass

    for dist in sorted(examples_by_distance.keys())[:5]:
        print(f"\nDistance {dist}:")
        for source, target, db in examples_by_distance[dist][:3]:
            print(f"  {source} -> {target} ({db})")

    conn.close()
    return distance_stats, unreachable_pairs, direct_connections, not_in_omnipath

if __name__ == "__main__":
    distance_stats, unreachable_pairs, direct_connections, not_in_omnipath = analyze_edge_distances()
