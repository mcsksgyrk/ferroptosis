import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import sqlite3
import pandas as pd
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from typing import List, Dict, Set
from pathlib import Path
import json


def get_core_arn_proteins(merged_db_path: Path) -> List[str]:
    conn = sqlite3.connect(merged_db_path)
    query = """
        SELECT DISTINCT n.name
        FROM node n
        WHERE n.source_db LIKE '%ARN%'
        AND n.type = 'protein'
        AND (n.pathways LIKE '%autophagy%' OR n.function LIKE '%autophagy%')
    """
    res = conn.execute(query).fetchall()
    conn.close()
    return [val[0] for val in res]


def get_core_fer_proteins(merged_db_path: Path) -> List[str]:
    conn = sqlite3.connect(merged_db_path)
    query = """
        SELECT DISTINCT n.name
        FROM node n
        WHERE (n.source_db LIKE '%KEGG%' OR n.source_db LIKE '%ferrdb%' OR n.source_db LIKE '%FerReg%')
        AND n.type = 'protein'
        AND n.role_in_ferroptosis != ''
    """
    res = conn.execute(query).fetchall()
    conn.close()
    return [val[0] for val in res]


def get_all_arn_proteins(merged_db_path: Path) -> List[str]:
    conn = sqlite3.connect(merged_db_path)
    query = """
        SELECT DISTINCT n.name
        FROM node n
        WHERE n.source_db LIKE '%ARN%'
        AND n.type = 'protein'
    """
    res = conn.execute(query).fetchall()
    conn.close()
    return [val[0] for val in res]


def get_all_fer_proteins(merged_db_path: Path) -> List[str]:
    conn = sqlite3.connect(merged_db_path)
    query = """
        SELECT DISTINCT n.name
        FROM node n
        WHERE (n.source_db LIKE '%KEGG%' OR n.source_db LIKE '%ferrdb%' OR n.source_db LIKE '%FerReg%')
        AND n.type = 'protein'
        AND n.primary_id_type = 'uniprot_id'
    """
    res = conn.execute(query).fetchall()
    conn.close()
    return [val[0] for val in res]


def find_node_edges(node: str, edges_df: pd.DataFrame, cond_set: Set = None) -> pd.DataFrame:
    if cond_set is not None:
        res = edges_df[
            (edges_df.interactor_a_node_name == node) &
            (edges_df.interactor_b_node_name.isin(cond_set))
        ]
    else:
        res = edges_df[
            (edges_df.interactor_a_node_name == node)
        ]
    return res


def analyze_downstream_connectivity(node_set: Set, edges_df: pd.DataFrame,
                                   arn_core: Set, fer_core: Set) -> Dict:
    shared_core = arn_core & fer_core
    arn_only_core = arn_core - shared_core
    fer_only_core = fer_core - shared_core

    res = {}
    for node in list(node_set):
        arn_conn = len(find_node_edges(node, edges_df, arn_only_core))
        fer_conn = len(find_node_edges(node, edges_df, fer_only_core))
        shared_conn = len(find_node_edges(node, edges_df, shared_core))
        all_conn = len(find_node_edges(node, edges_df))
        score = arn_conn + fer_conn + shared_conn
        crosstalk_score = min(arn_conn + shared_conn, fer_conn + shared_conn)

        res[node] = {
            'fer_conn': fer_conn,
            'arn_conn': arn_conn,
            'shared_conn': shared_conn,
            'all_conn': all_conn,
            'score': score,
            'crosstalk_score': crosstalk_score
        }
    return res


def main():
    merged_db_path = OUTPUTS_DIR / 'merged_ferroptosis_network.db'

    extended_db_path = OUTPUTS_DIR / 'extended_omnipath_network.db'
    if extended_db_path.exists():
        merged_db_path = extended_db_path

    print(f"Using database: {merged_db_path}")

    arn_core = set(get_core_arn_proteins(merged_db_path))
    fer_core = set(get_core_fer_proteins(merged_db_path))
    common_core = arn_core & fer_core

    print(f"Core ARN proteins: {len(arn_core)}")
    print(f"Core ferroptosis proteins: {len(fer_core)}")
    print(f"Shared core proteins: {len(common_core)}")

    all_arn = set(get_all_arn_proteins(merged_db_path))
    all_fer = set(get_all_fer_proteins(merged_db_path))

    plt.figure(figsize=(8, 6), dpi=150)
    venn2(subsets=(len(all_arn - all_fer), len(all_fer - all_arn), len(all_arn & all_fer)),
          set_labels=('Autophagy (ARN)', 'Ferroptosis'))
    plt.title('Protein Overlap Between Autophagy and Ferroptosis Networks')
    plt.savefig(OUTPUTS_DIR / 'arn_ferroptosis_venn.png', dpi=150, bbox_inches='tight')
    plt.show()

    conn = sqlite3.connect(merged_db_path)
    edges_df = pd.read_sql_query("SELECT * FROM edge", conn)
    conn.close()

    print(f"Total edges: {len(edges_df)}")

    # 1. Direct crosstalk: ARN_core → FER_core or FER_core → ARN_core
    crosstalk = edges_df[
        ((edges_df.interactor_a_node_name.isin(arn_core)) &
         (edges_df.interactor_b_node_name.isin(fer_core))) |
        ((edges_df.interactor_a_node_name.isin(fer_core)) &
         (edges_df.interactor_b_node_name.isin(arn_core)))
    ]
    print(f"Direct crosstalk edges: {len(crosstalk)}")

    # 2. Shared targets: ARN_core → X ← FER_core
    arn_targets = set(edges_df[
        edges_df.interactor_a_node_name.isin(arn_core)
    ].interactor_b_node_name)

    fer_targets = set(edges_df[
        edges_df.interactor_a_node_name.isin(fer_core)
    ].interactor_b_node_name)

    shared_targets = arn_targets & fer_targets
    print(f"Shared target proteins: {len(shared_targets)}")

    # 3. Shared regulators: X → ARN_core and X → FER_core
    arn_regulators = set(edges_df[
        edges_df.interactor_b_node_name.isin(arn_core)
    ].interactor_a_node_name)

    fer_regulators = set(edges_df[
        edges_df.interactor_b_node_name.isin(fer_core)
    ].interactor_a_node_name)

    shared_regulators = arn_regulators & fer_regulators
    print(f"Shared regulator proteins: {len(shared_regulators)}")

    res = analyze_downstream_connectivity(shared_regulators, edges_df, arn_core, fer_core)
    df_upstream_reg = pd.DataFrame.from_dict(res, orient='index')

    # Filter for significant crosstalk regulators
    significant_regulators = df_upstream_reg[
        (df_upstream_reg.crosstalk_score > 2) &
        (df_upstream_reg.fer_conn > 0) &
        (df_upstream_reg.arn_conn > 0)
    ].index.tolist()

    print(f"Significant crosstalk regulators: {len(significant_regulators)}")

    # Detailed analysis of significant regulators
    shared_core = arn_core & fer_core
    arn_only_core = arn_core - shared_core
    fer_only_core = fer_core - shared_core

    regulators = {}
    for node in significant_regulators:
        edges = find_node_edges(node, edges_df, fer_core | arn_core)

        if node not in regulators:
            regulators[node] = {
                'fer': [],
                'arn': [],
                'shared': []
            }

        for idx, edge in edges.iterrows():
            if edge.interactor_b_node_name in shared_core:
                regulators[node]['shared'].append(edge.interactor_b_node_name)
            elif edge.interactor_b_node_name in arn_only_core:
                regulators[node]['arn'].append(edge.interactor_b_node_name)
            elif edge.interactor_b_node_name in fer_only_core:
                regulators[node]['fer'].append(edge.interactor_b_node_name)

    # Remove duplicates
    for node in regulators:
        regulators[node]['fer'] = list(set(regulators[node]['fer']))
        regulators[node]['arn'] = list(set(regulators[node]['arn']))
        regulators[node]['shared'] = list(set(regulators[node]['shared']))

    # Save results
    output_path = OUTPUTS_DIR / "crosstalk_analysis.json"
    with open(output_path, "w") as f:
        json.dump({
            'statistics': {
                'arn_core_proteins': len(arn_core),
                'fer_core_proteins': len(fer_core),
                'shared_core_proteins': len(common_core),
                'all_arn_proteins': len(all_arn),
                'all_fer_proteins': len(all_fer),
                'direct_crosstalk_edges': len(crosstalk),
                'shared_targets': len(shared_targets),
                'shared_regulators': len(shared_regulators),
                'significant_crosstalk_regulators': len(significant_regulators)
            },
            'regulators': regulators
        }, f, indent=2)

    print(f"\nResults saved to {output_path}")

    # Save upstream regulator analysis
    df_upstream_reg.to_csv(OUTPUTS_DIR / "upstream_regulators_analysis.csv")

    return {
        'arn_core': arn_core,
        'fer_core': fer_core,
        'shared_targets': shared_targets,
        'shared_regulators': shared_regulators,
        'significant_regulators': significant_regulators,
        'crosstalk_edges': crosstalk
    }


if __name__ == "__main__":
    results = main()
