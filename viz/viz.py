import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from matplotlib_venn import venn3

OUTPUTS_DIR = Path('outputs')
KEGG_DB = OUTPUTS_DIR / 'kegg.db'
FERRDB_DB = OUTPUTS_DIR / 'ferrdb_network.db'
FERREG_DB = OUTPUTS_DIR / 'ferreg_network.db'
MERGED_DB = OUTPUTS_DIR / 'merged_ferroptosis_network.db'
EXTENDED_OMNIPATH_DB = OUTPUTS_DIR / 'merged_ferroptosis_w_omnipat.db'

plt.style.use('default')
colors = ['#E74C3C', '#2ECC71', '#3498DB', '#F39C12', '#9B59B6']

def get_db_stats(db_path):
    if not db_path.exists():
        return {'nodes': 0, 'edges': 0}

    conn = sqlite3.connect(db_path)
    try:
        nodes = conn.execute(
            "SELECT COUNT(*) FROM node WHERE type = 'protein' AND primary_id_type = 'uniprot_id'"
        ).fetchone()[0]
        edges = conn.execute("""
            SELECT COUNT(*) FROM edge e
            JOIN node n1 ON e.interactor_a_node_id = n1.id
            JOIN node n2 ON e.interactor_b_node_id = n2.id
            WHERE n1.type = 'protein' AND n2.type = 'protein'
            AND n1.primary_id_type = 'uniprot_id' AND n2.primary_id_type = 'uniprot_id'
        """).fetchone()[0]
    except:
        nodes = edges = 0
    conn.close()
    return {'nodes': nodes, 'edges': edges}

def get_protein_names(db_path):
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

def create_combined_figure():
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # A. Source Database Overlap (Venn Diagram)
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

    venn3(subsets=(only_kegg, only_ferrdb, kegg_ferrdb, only_ferreg, kegg_ferreg, ferrdb_ferreg, all_three),
          set_labels=('KEGG', 'FerrDB', 'FerReg'), set_colors=colors[:3], alpha=0.8, ax=ax1)
    ax1.set_title('A. Source Database Overlap', fontsize=14, fontweight='bold', pad=20)

    # B. Database Statistics Comparison
    kegg_stats = get_db_stats(KEGG_DB)
    ferrdb_stats = get_db_stats(FERRDB_DB)
    ferreg_stats = get_db_stats(FERREG_DB)
    merged_stats = get_db_stats(MERGED_DB)
    extended_stats = get_db_stats(EXTENDED_OMNIPATH_DB)

    databases = ['KEGG', 'FerrDB', 'FerReg', 'Merged', 'Extended\n+OmniPath']
    all_stats = [kegg_stats, ferrdb_stats, ferreg_stats, merged_stats, extended_stats]

    x = np.arange(len(databases))
    width = 0.35

    nodes_counts = [s['nodes'] for s in all_stats]
    edges_counts = [s['edges'] for s in all_stats]

    bars1 = ax2.bar(x - width/2, nodes_counts, width, label='UniProt Proteins',
                    color=colors[0], alpha=0.8, edgecolor='black', linewidth=1)
    bars2 = ax2.bar(x + width/2, edges_counts, width, label='PPI Edges',
                    color=colors[1], alpha=0.8, edgecolor='black', linewidth=1)

    ax2.set_xlabel('Database', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax2.set_title('B. Database Statistics', fontsize=14, fontweight='bold', pad=20)
    ax2.set_xticks(x)
    ax2.set_xticklabels(databases, fontsize=10, fontweight='bold')
    legend = ax2.legend(fontsize=11)
    for text in legend.get_texts():
        text.set_fontweight('bold')
    ax2.set_yscale('log')
    ax2.grid(axis='y', alpha=0.3)

    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax2.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                        f'{int(height):,}', ha='center', va='bottom', fontsize=10,
                        fontweight='bold', rotation=45)

    conn = sqlite3.connect(EXTENDED_OMNIPATH_DB)

    # C. Network Layer Distribution
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
        bars = ax3.bar(range(len(df_layers)), df_layers['count'],
                      color=colors[2], alpha=0.8, edgecolor='black', linewidth=1)
        ax3.set_xticks(range(len(df_layers)))
        ax3.set_xticklabels(df_layers['layer'], fontsize=11, fontweight='bold')
        ax3.set_ylabel('Number of PPI Edges', fontsize=12, fontweight='bold')
        ax3.set_title('C. Network Layer Distribution', fontsize=14, fontweight='bold', pad=20)
        ax3.set_yscale('log')
        ax3.grid(axis='y', alpha=0.3)

        for bar in bars:
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{int(height):,}', ha='center', va='bottom',
                    fontsize=10, fontweight='bold')

    # D. Source Database Contribution
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
            total_edges = sum(source_counts.values())
            filtered_counts = {k: v for k, v in source_counts.items() if (v/total_edges) >= 0.01}

            wedges, texts, autotexts = ax4.pie(filtered_counts.values(),
                                              labels=filtered_counts.keys(),
                                              autopct='%1.1f%%',
                                              startangle=90,
                                              colors=colors[:len(filtered_counts)],
                                              explode=[0.05] * len(filtered_counts))

            for text in texts:
                text.set_fontsize(11)
                text.set_fontweight('bold')
            for autotext in autotexts:
                autotext.set_fontsize(10)
                autotext.set_fontweight('bold')
                autotext.set_color('white')

            ax4.set_title('D. PPI Edge Sources', fontsize=14, fontweight='bold', pad=20)

    conn.close()

    fig.suptitle('Ferroptosis Network Database Integration Analysis',
                fontsize=16, fontweight='bold', y=0.96)

    plt.tight_layout()
    plt.subplots_adjust(top=0.90, wspace=0.25, hspace=0.35)

    return fig


if __name__ == "__main__":
    fig = create_combined_figure()
    plt.savefig(OUTPUTS_DIR / 'figures/combined_network_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

    print("=== Publication-style figure saved ===")
