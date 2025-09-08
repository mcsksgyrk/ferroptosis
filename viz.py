import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from matplotlib_venn import venn3, venn3_circles
import networkx as nx

OUTPUTS_DIR = Path('outputs')
KEGG_DB = OUTPUTS_DIR / 'kegg.db'
FERRDB_DB = OUTPUTS_DIR / 'ferrdb_network.db'
FERREG_DB = OUTPUTS_DIR / 'ferreg_network.db'
MERGED_DB = OUTPUTS_DIR / 'merged_ferroptosis_network.db'
EXTENDED_OMNIPATH_DB = OUTPUTS_DIR / 'extended_omnipath_network.db'

plt.style.use('seaborn-whitegrid')
sns.set_palette("husl")

def get_db_stats(db_path):
    """Extract basic statistics from a database - UniProt proteins only"""
    if not db_path.exists():
        return {'nodes': 0, 'edges': 0, 'diseases': 0}

    conn = sqlite3.connect(db_path)
    stats = {}

    try:
        stats['nodes'] = conn.execute(
            "SELECT COUNT(*) FROM node WHERE type = 'protein' AND primary_id_type = 'uniprot_id'"
        ).fetchone()[0]
    except:
        stats['nodes'] = 0

    try:
        query = """
        SELECT COUNT(*)
        FROM edge e
        JOIN node n1 ON e.interactor_a_node_id = n1.id
        JOIN node n2 ON e.interactor_b_node_id = n2.id
        WHERE n1.type = 'protein' AND n2.type = 'protein'
        AND n1.primary_id_type = 'uniprot_id' AND n2.primary_id_type = 'uniprot_id'
        """
        stats['edges'] = conn.execute(query).fetchone()[0]
    except:
        stats['edges'] = 0

    try:
        stats['diseases'] = conn.execute("SELECT COUNT(*) FROM disease").fetchone()[0]
    except:
        stats['diseases'] = 0

    conn.close()
    return stats

def get_protein_names(db_path):
    """Get all UniProt protein names from database"""
    if not db_path.exists():
        return set()

    conn = sqlite3.connect(db_path)
    try:
        nodes = conn.execute(
            "SELECT DISTINCT name FROM node WHERE type = 'protein' AND primary_id_type = 'uniprot_id'"
        ).fetchall()
        node_set = {node[0] for node in nodes}
    except:
        node_set = set()
    conn.close()
    return node_set

def create_figure_1():
    """Create Figure 1: Database Integration Overview including OmniPath Extension"""
    fig = plt.figure(figsize=(20, 6))

    kegg_stats = get_db_stats(KEGG_DB)
    ferrdb_stats = get_db_stats(FERRDB_DB)
    ferreg_stats = get_db_stats(FERREG_DB)
    merged_stats = get_db_stats(MERGED_DB)
    extended_stats = get_db_stats(EXTENDED_OMNIPATH_DB)

    # A. Integration Flow
    ax1 = plt.subplot(141)

    sources = ['KEGG', 'FerrDB', 'FerReg']
    source_nodes = [kegg_stats['nodes'], ferrdb_stats['nodes'], ferreg_stats['nodes']]
    source_edges = [kegg_stats['edges'], ferrdb_stats['edges'], ferreg_stats['edges']]

    y_positions = [0.8, 0.6, 0.4]
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']

    for i, (source, nodes, edges, y, color) in enumerate(zip(sources, source_nodes, source_edges, y_positions, colors)):
        rect = Rectangle((0.05, y-0.05), 0.15, 0.08, facecolor=color, alpha=0.7)
        ax1.add_patch(rect)
        ax1.text(0.125, y, f'{source}\n{nodes}\n{edges}',
                ha='center', va='center', fontsize=8, fontweight='bold')

        ax1.arrow(0.22, y, 0.13, 0.6-y, head_width=0.015, head_length=0.015,
                 fc=color, ec=color, alpha=0.5, width=0.003)

    rect = Rectangle((0.38, 0.55), 0.12, 0.1, facecolor='#95E77E', alpha=0.7)
    ax1.add_patch(rect)
    ax1.text(0.44, 0.6, f'Merged\n{merged_stats["nodes"]}\n{merged_stats["edges"]}',
            ha='center', va='center', fontsize=8, fontweight='bold')

    ax1.arrow(0.52, 0.6, 0.18, 0, head_width=0.02, head_length=0.02,
             fc='#95E77E', ec='#95E77E', alpha=0.7, width=0.005)

    rect = Rectangle((0.72, 0.55), 0.22, 0.1, facecolor='#FFD93D', alpha=0.7)
    ax1.add_patch(rect)
    ax1.text(0.83, 0.6, f'Extended\n+ OmniPath\n{extended_stats["nodes"]} | {extended_stats["edges"]}',
            ha='center', va='center', fontsize=8, fontweight='bold')

    ax1.set_xlim(0, 1)
    ax1.set_ylim(0.2, 1)
    ax1.axis('off')
    ax1.set_title('A. Integration Flow', fontsize=12, fontweight='bold')

    # B. Source Database Overlap
    ax2 = plt.subplot(142)

    kegg_proteins = get_protein_names(KEGG_DB)
    ferrdb_proteins = get_protein_names(FERRDB_DB)
    ferreg_proteins = get_protein_names(FERREG_DB)

    only_kegg = len(kegg_proteins - ferrdb_proteins - ferreg_proteins)
    only_ferrdb = len(ferrdb_proteins - kegg_proteins - ferreg_proteins)
    only_ferreg = len(ferreg_proteins - kegg_proteins - ferrdb_proteins)
    kegg_ferrdb = len((kegg_proteins & ferrdb_proteins) - ferreg_proteins)
    kegg_ferreg = len((kegg_proteins & ferreg_proteins) - ferrdb_proteins)
    ferrdb_ferreg = len((ferrdb_proteins & ferreg_proteins) - kegg_proteins)
    all_three = len(kegg_proteins & ferrdb_proteins & ferreg_proteins)

    venn = venn3(subsets=(only_kegg, only_ferrdb, kegg_ferrdb, only_ferreg,
                          kegg_ferreg, ferrdb_ferreg, all_three),
                 set_labels=('KEGG', 'FerrDB', 'FerReg'),
                 ax=ax2)

    ax2.set_title('B. Source Database Overlap', fontsize=12, fontweight='bold')

    # C. Database Statistics Comparison
    ax3 = plt.subplot(143)

    databases = ['KEGG', 'FerrDB', 'FerReg', 'Merged', 'Extended\n+OmniPath']
    all_stats = [kegg_stats, ferrdb_stats, ferreg_stats, merged_stats, extended_stats]

    x = np.arange(len(databases))
    width = 0.35

    nodes_counts = [s['nodes'] for s in all_stats]
    edges_counts = [s['edges'] for s in all_stats]

    bars1 = ax3.bar(x - width/2, nodes_counts, width, label='UniProt Proteins', color='#FF6B6B')
    bars2 = ax3.bar(x + width/2, edges_counts, width, label='PPI Edges', color='#4ECDC4')

    ax3.set_xlabel('Database', fontsize=11)
    ax3.set_ylabel('Count', fontsize=11)
    ax3.set_title('C. Database Statistics', fontsize=12, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels(databases, fontsize=9)
    ax3.legend()
    ax3.set_yscale('log')

    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax3.text(bar.get_x() + bar.get_width()/2., height,
                        f'{int(height)}', ha='center', va='bottom', fontsize=7, rotation=45)

    # D. Network Growth
    ax4 = plt.subplot(144)

    steps = ['KEGG', 'FerrDB', 'FerReg', 'Merged', 'Extended']
    protein_growth = [kegg_stats['nodes'], ferrdb_stats['nodes'], ferreg_stats['nodes'],
                     merged_stats['nodes'], extended_stats['nodes']]
    edge_growth = [kegg_stats['edges'], ferrdb_stats['edges'], ferreg_stats['edges'],
                  merged_stats['edges'], extended_stats['edges']]

    ax4.plot(steps, protein_growth, 'o-', color='#FF6B6B', linewidth=2, markersize=6, label='Proteins')
    ax4.plot(steps, edge_growth, 's-', color='#4ECDC4', linewidth=2, markersize=6, label='PPI Edges')

    ax4.set_ylabel('Count', fontsize=11)
    ax4.set_title('D. Network Growth', fontsize=12, fontweight='bold')
    ax4.legend()
    ax4.set_yscale('log')
    ax4.tick_params(axis='x', rotation=45)

    plt.suptitle('Figure 1: Protein-Protein Interaction Database Integration with OmniPath Extension',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    return fig

def create_figure_4():
    """Create Figure 4: Extended OmniPath Network Analysis"""
    fig = plt.figure(figsize=(16, 10))

    conn = sqlite3.connect(EXTENDED_OMNIPATH_DB)

    # A. Network Layer Distribution
    ax1 = plt.subplot(221)

    query = """
    SELECT e.layer, COUNT(*) as count
    FROM edge e
    JOIN node n1 ON e.interactor_a_node_id = n1.id
    JOIN node n2 ON e.interactor_b_node_id = n2.id
    WHERE n1.type = 'protein' AND n2.type = 'protein'
    AND n1.primary_id_type = 'uniprot_id' AND n2.primary_id_type = 'uniprot_id'
    GROUP BY e.layer
    ORDER BY count DESC
    """

    df_layers = pd.read_sql_query(query, conn)

    if not df_layers.empty:
        bars = ax1.bar(range(len(df_layers)), df_layers['count'], color='#FF6B6B')
        ax1.set_xticks(range(len(df_layers)))
        ax1.set_xticklabels(df_layers['layer'], rotation=45, ha='right')
        ax1.set_ylabel('Number of PPI Edges', fontsize=11)
        ax1.set_title('A. Network Layer Distribution', fontsize=12, fontweight='bold')
        ax1.set_yscale('log')

        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}', ha='center', va='bottom', fontsize=8)

    # B. Source Database Contribution
    ax2 = plt.subplot(222)

    query = """
    SELECT e.source_db, COUNT(*) as count
    FROM edge e
    JOIN node n1 ON e.interactor_a_node_id = n1.id
    JOIN node n2 ON e.interactor_b_node_id = n2.id
    WHERE n1.type = 'protein' AND n2.type = 'protein'
    AND n1.primary_id_type = 'uniprot_id' AND n2.primary_id_type = 'uniprot_id'
    GROUP BY e.source_db
    ORDER BY count DESC
    """

    df_sources = pd.read_sql_query(query, conn)

    if not df_sources.empty:
        source_counts = {}
        for idx, row in df_sources.iterrows():
            sources = str(row['source_db']).split('|')
            for source in sources:
                source = source.strip()
                if source and source != '-':
                    source_counts[source] = source_counts.get(source, 0) + row['count']

        if source_counts:
            plt.pie(source_counts.values(), labels=source_counts.keys(),
                   autopct='%1.1f%%', startangle=90)
            ax2.set_title('B. PPI Edge Sources', fontsize=12, fontweight='bold')

    # C. Ferroptosis Proteins in Extended Network
    ax3 = plt.subplot(223)

    query = """
    SELECT role_in_ferroptosis, COUNT(*) as count
    FROM node
    WHERE type = 'protein'
    AND primary_id_type = 'uniprot_id'
    AND role_in_ferroptosis IS NOT NULL
    AND role_in_ferroptosis != ''
    GROUP BY role_in_ferroptosis
    """

    df_roles = pd.read_sql_query(query, conn)

    if not df_roles.empty:
        role_counts = {}
        for idx, row in df_roles.iterrows():
            roles = str(row['role_in_ferroptosis']).split('|')
            for role in roles:
                role = role.strip()
                if role and role != '-':
                    role_counts[role] = role_counts.get(role, 0) + row['count']

        if role_counts:
            roles = list(role_counts.keys())
            counts = list(role_counts.values())

            bars = ax3.bar(range(len(roles)), counts, color='#4ECDC4')
            ax3.set_xticks(range(len(roles)))
            ax3.set_xticklabels(roles, rotation=45, ha='right')
            ax3.set_ylabel('Number of Proteins', fontsize=11)
            ax3.set_title('C. Ferroptosis Proteins', fontsize=12, fontweight='bold')

            for bar in bars:
                height = bar.get_height()
                ax3.text(bar.get_x() + bar.get_width()/2., height,
                        f'{int(height)}', ha='center', va='bottom')

    # D. Network Connectivity Comparison
    ax4 = plt.subplot(224)

    databases = ['Merged', 'Extended + OmniPath']
    db_paths = [MERGED_DB, EXTENDED_OMNIPATH_DB]

    connectivity_data = []
    for db_path in db_paths:
        if db_path.exists():
            conn_temp = sqlite3.connect(db_path)

            query = """
            SELECT COUNT(*) as edge_count
            FROM edge e
            JOIN node n1 ON e.interactor_a_node_id = n1.id
            JOIN node n2 ON e.interactor_b_node_id = n2.id
            WHERE n1.type = 'protein' AND n2.type = 'protein'
            AND n1.primary_id_type = 'uniprot_id' AND n2.primary_id_type = 'uniprot_id'
            """
            edges = conn_temp.execute(query).fetchone()[0]

            query = """
            SELECT COUNT(*) as node_count
            FROM node
            WHERE type = 'protein' AND primary_id_type = 'uniprot_id'
            """
            nodes = conn_temp.execute(query).fetchone()[0]

            avg_degree = (2 * edges) / nodes if nodes > 0 else 0
            connectivity_data.append(avg_degree)
            conn_temp.close()
        else:
            connectivity_data.append(0)

    bars = ax4.bar(databases, connectivity_data, color=['#45B7D1', '#FFD93D'])
    ax4.set_ylabel('Average Degree', fontsize=11)
    ax4.set_title('D. Network Connectivity', fontsize=12, fontweight='bold')

    for bar in bars:
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}', ha='center', va='bottom')

    conn.close()

    plt.suptitle('Figure 4: Extended OmniPath Protein-Protein Interaction Network Analysis',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    return fig

if __name__ == "__main__":
    fig1 = create_figure_1()
    plt.savefig('figure_1_extended_omnipath_integration.png', dpi=300, bbox_inches='tight')
    plt.savefig('figure_1_extended_omnipath_integration.pdf', bbox_inches='tight')
    plt.show()

    fig4 = create_figure_4()
    plt.savefig('figure_4_extended_omnipath_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig('figure_4_extended_omnipath_analysis.pdf', bbox_inches='tight')
    plt.show()

    print("Extended OmniPath figures saved as PNG and PDF files")
