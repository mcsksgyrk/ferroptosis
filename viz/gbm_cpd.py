from database.external_db import DBconnector
from config import OUTPUTS_DIR

db = DBconnector(OUTPUTS_DIR / 'merged_ferroptosis_w_omnipath.db')

# 1. All GBM-associated edges with full node info
gbm_edges = db.query_to_dataframe("""
    SELECT
        e.interactor_a_node_name as source_id,
        e.interactor_b_node_name as target_id,
        n1.display_name as source_display,
        n2.display_name as target_display,
        n1.type as source_type,
        n2.type as target_type,
        n1.role_in_ferroptosis as source_role,
        n2.role_in_ferroptosis as target_role,
        e.interaction_types,
        e.layer,
        e.source_db
    FROM disease_edge de
    JOIN disease d ON de.disease_id = d.id
    JOIN edge e ON de.edge_id = e.id
    JOIN node n1 ON e.interactor_a_node_name = n1.name
    JOIN node n2 ON e.interactor_b_node_name = n2.name
    WHERE d.disease_id = 'ICD-11: 2A00'
""")

# 2. Unique GBM nodes categorized by type
gbm_nodes = db.query_to_dataframe("""
    SELECT DISTINCT
        n.name, n.display_name, n.type,
        n.role_in_ferroptosis, n.source_db
    FROM disease_edge de
    JOIN disease d ON de.disease_id = d.id
    JOIN edge e ON de.edge_id = e.id
    JOIN node n ON n.name IN (e.interactor_a_node_name, e.interactor_b_node_name)
    WHERE d.disease_id = 'ICD-11: 2A00'
""")

# 3. Which GBM proteins sit in which ferroptosis layer
gbm_protein_layers = db.query_to_dataframe("""
    SELECT DISTINCT
        n.name, n.display_name, n.role_in_ferroptosis, e.layer,
        CASE WHEN n.source_db LIKE '%KEGG%' THEN 1 ELSE 0 END as is_kegg_core
    FROM disease_edge de
    JOIN disease d ON de.disease_id = d.id
    JOIN edge e ON de.edge_id = e.id
    JOIN node n ON n.name IN (e.interactor_a_node_name, e.interactor_b_node_name)
    WHERE d.disease_id = 'ICD-11: 2A00'
    AND n.type = 'protein'
""")

# 4. Compound -> direct target (within GBM subnetwork)
gbm_compound_targets = db.query_to_dataframe("""
    SELECT DISTINCT
        n1.display_name as compound,
        n2.display_name as direct_target,
        n2.name as target_uniprot,
        n2.role_in_ferroptosis as target_role,
        e.interaction_types as compound_effect,
        e.layer
    FROM disease_edge de
    JOIN disease d ON de.disease_id = d.id
    JOIN edge e ON de.edge_id = e.id
    JOIN node n1 ON e.interactor_a_node_name = n1.name
    JOIN node n2 ON e.interactor_b_node_name = n2.name
    WHERE d.disease_id = 'ICD-11: 2A00'
    AND n1.type = 'small_molecule'
""")

# 5. Downstream of compound targets
target_uniprots = gbm_compound_targets.target_uniprot.unique().tolist()
if target_uniprots:
    placeholders = ','.join([f"'{uid}'" for uid in target_uniprots])
    downstream = db.query_to_dataframe(f"""
        SELECT DISTINCT
            n_src.display_name as source_protein,
            n_dst.display_name as downstream_target,
            n_dst.type as downstream_type,
            n_dst.role_in_ferroptosis as downstream_role,
            e.interaction_types,
            e.layer
        FROM edge e
        JOIN node n_src ON e.interactor_a_node_name = n_src.name
        JOIN node n_dst ON e.interactor_b_node_name = n_dst.name
        WHERE e.interactor_a_node_name IN ({placeholders})
    """)

print(f"GBM edges: {len(gbm_edges)}")
print(f"GBM nodes: {len(gbm_nodes)}")
print(f"Compound targets: {len(gbm_compound_targets)}")
print(f"Downstream: {len(downstream)}")
print(f"Node types:\n{gbm_nodes.type.value_counts()}")


import networkx as nx
import matplotlib.pyplot as plt

color_map = {
    'protein': '#4C72B0',
    'small_molecule': '#DD8452',
    'miRNA': '#55A868',
    'lncRNA': '#C44E52',
    'not sure': '#8C8C8C'
}

# --- Figure 1: Full GBM subnetwork ---

G1 = nx.DiGraph()

node_type_lookup = dict(zip(gbm_nodes.name, gbm_nodes.type))
node_display_lookup = dict(zip(gbm_nodes.name, gbm_nodes.display_name))

for _, row in gbm_nodes.iterrows():
    G1.add_node(row.display_name, type=row.type)

for _, row in gbm_edges.iterrows():
    G1.add_edge(row.source_display, row.target_display)

node_colors = [color_map.get(G1.nodes[n].get('type', ''), '#8C8C8C') for n in G1.nodes()]

fig1, ax1 = plt.subplots(figsize=(16, 12))
pos = nx.spring_layout(G1, k=2, seed=42)
nx.draw_networkx_nodes(G1, pos, node_color=node_colors, node_size=300, ax=ax1)
nx.draw_networkx_edges(G1, pos, arrows=True, edge_color='#CCCCCC', ax=ax1)
nx.draw_networkx_labels(G1, pos, font_size=6, ax=ax1)

legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=10, label=t)
                  for t, c in color_map.items()]
ax1.legend(handles=legend_handles, loc='upper left')
ax1.set_title('GBM-associated ferroptosis subnetwork')
plt.tight_layout()
plt.savefig('gbm_full_network.png', dpi=200)
plt.show()

# --- Figure 2: Compound chains (only downstream within GBM subnetwork) ---

gbm_node_names = set(gbm_nodes.display_name)
downstream_filtered = downstream[downstream.downstream_target.isin(gbm_node_names)]

G2 = nx.DiGraph()

for _, row in gbm_compound_targets.iterrows():
    G2.add_node(row.compound, type='small_molecule')
    G2.add_node(row.direct_target, type='protein')
    G2.add_edge(row.compound, row.direct_target)

for _, row in downstream_filtered.iterrows():
    if row.source_protein in G2.nodes():
        G2.add_node(row.downstream_target, type=row.downstream_type or 'protein')
        G2.add_edge(row.source_protein, row.downstream_target)

node_colors_2 = [color_map.get(G2.nodes[n].get('type', ''), '#8C8C8C') for n in G2.nodes()]

fig2, ax2 = plt.subplots(figsize=(16, 12))
pos2 = nx.spring_layout(G2, k=2, seed=42)
nx.draw_networkx_nodes(G2, pos2, node_color=node_colors_2, node_size=300, ax=ax2)
nx.draw_networkx_edges(G2, pos2, arrows=True, edge_color='#CCCCCC', ax=ax2)
nx.draw_networkx_labels(G2, pos2, font_size=6, ax=ax2)

ax2.legend(handles=legend_handles, loc='upper left')
ax2.set_title('GBM compound → target → downstream (within GBM subnetwork)')
plt.tight_layout()
plt.savefig('gbm_compound_chains.png', dpi=200)
plt.show()

print(f"Figure 1: {G1.number_of_nodes()} nodes, {G1.number_of_edges()} edges")
print(f"Figure 2: {G2.number_of_nodes()} nodes, {G2.number_of_edges()} edges")
print(f"Downstream edges kept (within GBM): {len(downstream_filtered)} of {len(downstream)}")


# approach 2
lipid_perox_ids = ['O60488', 'P16050', 'Q6P1A2', 'Q15366']
gpx4_branch_ids = ['P36969', 'Q9UPY5']
all_ids = lipid_perox_ids + gpx4_branch_ids

all_placeholders = ','.join([f"'{uid}'" for uid in all_ids])
branch_proteins = db.query_to_dataframe(f"""
    SELECT name, display_name, type, role_in_ferroptosis, source_db
    FROM node
    WHERE name IN ({all_placeholders})
""")
print(branch_proteins)

lipid_placeholders = ','.join([f"'{uid}'" for uid in lipid_perox_ids])
gpx4_placeholders = ','.join([f"'{uid}'" for uid in gpx4_branch_ids])

gbm_compound_names = gbm_compound_targets.compound.unique().tolist()
compound_name_placeholders = ','.join([f"'{c}'" for c in gbm_compound_names])

direct_to_lipid = db.query_to_dataframe(f"""
    SELECT DISTINCT
        n1.display_name as compound,
        n2.display_name as target,
        e.interaction_types,
        e.layer
    FROM edge e
    JOIN node n1 ON e.interactor_a_node_name = n1.name
    JOIN node n2 ON e.interactor_b_node_name = n2.name
    WHERE n1.display_name IN ({compound_name_placeholders})
    AND e.interactor_b_node_name IN ({lipid_placeholders})
""")

direct_to_gpx4 = db.query_to_dataframe(f"""
    SELECT DISTINCT
        n1.display_name as compound,
        n2.display_name as target,
        e.interaction_types,
        e.layer
    FROM edge e
    JOIN node n1 ON e.interactor_a_node_name = n1.name
    JOIN node n2 ON e.interactor_b_node_name = n2.name
    WHERE n1.display_name IN ({compound_name_placeholders})
    AND e.interactor_b_node_name IN ({gpx4_placeholders})
""")

print(f"\n=== Direct compound -> lipid peroxidation ===")
print(direct_to_lipid)
print(f"\n=== Direct compound -> GPX4 branch ===")
print(direct_to_gpx4)

indirect_to_lipid = db.query_to_dataframe(f"""
    SELECT DISTINCT
        n_cpd.display_name as compound,
        n_mid.display_name as intermediate,
        n_lp.display_name as lipid_perox_target,
        e1.interaction_types as compound_effect,
        e2.interaction_types as intermediate_effect
    FROM edge e1
    JOIN edge e2 ON e1.interactor_b_node_name = e2.interactor_a_node_name
    JOIN node n_cpd ON e1.interactor_a_node_name = n_cpd.name
    JOIN node n_mid ON e1.interactor_b_node_name = n_mid.name
    JOIN node n_lp ON e2.interactor_b_node_name = n_lp.name
    WHERE n_cpd.display_name IN ({compound_name_placeholders})
    AND e2.interactor_b_node_name IN ({lipid_placeholders})
""")

indirect_to_gpx4 = db.query_to_dataframe(f"""
    SELECT DISTINCT
        n_cpd.display_name as compound,
        n_mid.display_name as intermediate,
        n_gpx.display_name as gpx4_branch_target,
        e1.interaction_types as compound_effect,
        e2.interaction_types as intermediate_effect
    FROM edge e1
    JOIN edge e2 ON e1.interactor_b_node_name = e2.interactor_a_node_name
    JOIN node n_cpd ON e1.interactor_a_node_name = n_cpd.name
    JOIN node n_mid ON e1.interactor_b_node_name = n_mid.name
    JOIN node n_gpx ON e2.interactor_b_node_name = n_gpx.name
    WHERE n_cpd.display_name IN ({compound_name_placeholders})
    AND e2.interactor_b_node_name IN ({gpx4_placeholders})
""")

print(f"\n=== Indirect compound -> ? -> lipid peroxidation ===")
print(f"Chains found: {len(indirect_to_lipid)}")
print(indirect_to_lipid)
print(f"\n=== Indirect compound -> ? -> GPX4 branch ===")
print(f"Chains found: {len(indirect_to_gpx4)}")
print(indirect_to_gpx4)

compounds_hitting_lipid = set(direct_to_lipid.compound) | set(indirect_to_lipid.compound)
compounds_hitting_gpx4 = set(direct_to_gpx4.compound) | set(indirect_to_gpx4.compound)
both = compounds_hitting_lipid & compounds_hitting_gpx4
only_lipid = compounds_hitting_lipid - compounds_hitting_gpx4
only_gpx4 = compounds_hitting_gpx4 - compounds_hitting_lipid
neither = set(gbm_compound_names) - compounds_hitting_lipid - compounds_hitting_gpx4

print(f"\n=== Compound branch comparison ===")
print(f"Hit BOTH branches: {both}")
print(f"Only lipid peroxidation: {only_lipid}")
print(f"Only GPX4 branch: {only_gpx4}")
print(f"Neither: {neither}")
direct_to_gpx4.to_csv('uhoh.csv')

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def is_direct_interaction(interaction_types):
    if not interaction_types:
        return False
    return 'is_direct:true' in interaction_types


color_map = {
    'protein': '#4C72B0',
    'small_molecule': '#DD8452',
    'miRNA': '#55A868',
    'lncRNA': '#C44E52',
    'not sure': '#8C8C8C'
}

# --- Figure 1: Full GBM subnetwork ---

G1 = nx.DiGraph()

for _, row in gbm_nodes.iterrows():
    G1.add_node(row.display_name, type=row.type)

direct_edges_1 = []
indirect_edges_1 = []

for _, row in gbm_edges.iterrows():
    G1.add_edge(row.source_display, row.target_display)
    edge = (row.source_display, row.target_display)
    if is_direct_interaction(row.interaction_types):
        direct_edges_1.append(edge)
    else:
        indirect_edges_1.append(edge)

node_colors_1 = [color_map.get(G1.nodes[n].get('type', ''), '#8C8C8C') for n in G1.nodes()]

fig1, ax1 = plt.subplots(figsize=(16, 12))
pos1 = nx.spring_layout(G1, k=2, seed=42)
nx.draw_networkx_nodes(G1, pos1, node_color=node_colors_1, node_size=300, ax=ax1)
nx.draw_networkx_edges(G1, pos1, edgelist=direct_edges_1, arrows=True, edge_color='#333333', width=2, ax=ax1)
nx.draw_networkx_edges(G1, pos1, edgelist=indirect_edges_1, arrows=True, edge_color='#AAAAAA', width=1, style='dashed', ax=ax1)
nx.draw_networkx_labels(G1, pos1, font_size=6, ax=ax1)

legend_1 = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=10, label=t)
            for t, c in color_map.items()]
legend_1.append(plt.Line2D([0], [0], color='#333333', linewidth=2, label='Direct interaction'))
legend_1.append(plt.Line2D([0], [0], color='#AAAAAA', linewidth=1, linestyle='dashed', label='Indirect interaction'))
ax1.legend(handles=legend_1, loc='upper left')
ax1.set_title('GBM-associated ferroptosis subnetwork')
plt.tight_layout()
plt.savefig('gbm_full_network.png', dpi=200)
plt.show()

# --- Figure 2: Compound chains within GBM subnetwork ---

gbm_node_names = set(gbm_nodes.display_name)
downstream_filtered = downstream[downstream.downstream_target.isin(gbm_node_names)]

G2 = nx.DiGraph()
direct_edges_2 = []
indirect_edges_2 = []

for _, row in gbm_compound_targets.iterrows():
    G2.add_node(row.compound, type='small_molecule')
    G2.add_node(row.direct_target, type='protein')
    edge = (row.compound, row.direct_target)
    G2.add_edge(*edge)
    if is_direct_interaction(row.compound_effect):
        direct_edges_2.append(edge)
    else:
        indirect_edges_2.append(edge)

for _, row in downstream_filtered.iterrows():
    if row.source_protein in G2.nodes():
        G2.add_node(row.downstream_target, type=row.downstream_type or 'protein')
        edge = (row.source_protein, row.downstream_target)
        G2.add_edge(*edge)
        if is_direct_interaction(row.interaction_types):
            direct_edges_2.append(edge)
        else:
            indirect_edges_2.append(edge)

node_colors_2 = [color_map.get(G2.nodes[n].get('type', ''), '#8C8C8C') for n in G2.nodes()]

fig2, ax2 = plt.subplots(figsize=(16, 12))
pos2 = nx.spring_layout(G2, k=2, seed=42)
nx.draw_networkx_nodes(G2, pos2, node_color=node_colors_2, node_size=300, ax=ax2)
nx.draw_networkx_edges(G2, pos2, edgelist=direct_edges_2, arrows=True, edge_color='#333333', width=2, ax=ax2)
nx.draw_networkx_edges(G2, pos2, edgelist=indirect_edges_2, arrows=True, edge_color='#AAAAAA', width=1, style='dashed', ax=ax2)
nx.draw_networkx_labels(G2, pos2, font_size=6, ax=ax2)

ax2.legend(handles=legend_1, loc='upper left')
ax2.set_title('GBM compound -> target -> downstream (within GBM subnetwork)')
plt.tight_layout()
plt.savefig('gbm_compound_chains.png', dpi=200)
plt.show()

# --- Figure 3: Branch comparison ---

group_colors = {
    'lipid_perox': '#E74C3C',
    'gpx4_branch': '#4C72B0',
    'both': '#7B2D8E',
    'only_gpx4': '#3498DB',
    'neither': '#8C8C8C'
}

G3 = nx.DiGraph()

for _, row in branch_proteins.iterrows():
    if row['name'] in lipid_perox_ids:
        G3.add_node(row.display_name, group='lipid_perox')
    else:
        G3.add_node(row.display_name, group='gpx4_branch')

all_gbm_cpds = set(gbm_compound_names)
for cpd in all_gbm_cpds:
    if cpd in both:
        G3.add_node(cpd, group='both')
    elif cpd in only_gpx4:
        G3.add_node(cpd, group='only_gpx4')
    else:
        G3.add_node(cpd, group='neither')

direct_edges_3 = []
indirect_edges_3 = []

for _, row in direct_to_lipid.iterrows():
    edge = (row.compound, row.target)
    G3.add_edge(*edge)
    if is_direct_interaction(row.interaction_types):
        direct_edges_3.append(edge)
    else:
        indirect_edges_3.append(edge)

for _, row in direct_to_gpx4.iterrows():
    edge = (row.compound, row.target)
    G3.add_edge(*edge)
    if is_direct_interaction(row.interaction_types):
        direct_edges_3.append(edge)
    else:
        indirect_edges_3.append(edge)

added_indirect = set()
for _, row in indirect_to_lipid.iterrows():
    edge = (row.compound, row.lipid_perox_target)
    if edge not in [(e[0], e[1]) for e in direct_edges_3] and edge not in added_indirect:
        G3.add_edge(*edge)
        indirect_edges_3.append(edge)
        added_indirect.add(edge)

for _, row in indirect_to_gpx4.iterrows():
    edge = (row.compound, row.gpx4_branch_target)
    if edge not in [(e[0], e[1]) for e in direct_edges_3] and edge not in added_indirect:
        G3.add_edge(*edge)
        indirect_edges_3.append(edge)
        added_indirect.add(edge)

node_colors_3 = [group_colors.get(G3.nodes[n].get('group', ''), '#8C8C8C') for n in G3.nodes()]

fig3, ax3 = plt.subplots(figsize=(16, 12))
pos3 = nx.spring_layout(G3, k=2.5, seed=42)
nx.draw_networkx_nodes(G3, pos3, node_color=node_colors_3, node_size=400, ax=ax3)
nx.draw_networkx_edges(G3, pos3, edgelist=direct_edges_3, arrows=True, edge_color='#333333', width=2, ax=ax3)
nx.draw_networkx_edges(G3, pos3, edgelist=indirect_edges_3, arrows=True, edge_color='#AAAAAA', width=1, style='dashed', ax=ax3)
nx.draw_networkx_labels(G3, pos3, font_size=7, ax=ax3)

legend_3 = [
    mpatches.Patch(color='#E74C3C', label='Lipid peroxidation target'),
    mpatches.Patch(color='#4C72B0', label='GPX4 branch target'),
    mpatches.Patch(color='#7B2D8E', label='Compound: hits both branches'),
    mpatches.Patch(color='#3498DB', label='Compound: only GPX4 branch'),
    mpatches.Patch(color='#8C8C8C', label='Compound: neither branch'),
    plt.Line2D([0], [0], color='#333333', linewidth=2, label='Direct interaction'),
    plt.Line2D([0], [0], color='#AAAAAA', linewidth=1, linestyle='dashed', label='Indirect interaction'),
]
ax3.legend(handles=legend_3, loc='upper left')
ax3.set_title('GBM compounds: GPX4 defense branch vs lipid peroxidation branch')
plt.tight_layout()
plt.savefig('gbm_branch_comparison.png', dpi=200)
plt.show()

print(f"Figure 1: {G1.number_of_nodes()} nodes, {G1.number_of_edges()} edges ({len(direct_edges_1)} direct, {len(indirect_edges_1)} indirect)")
print(f"Figure 2: {G2.number_of_nodes()} nodes, {G2.number_of_edges()} edges ({len(direct_edges_2)} direct, {len(indirect_edges_2)} indirect)")
print(f"Figure 3: {G3.number_of_nodes()} nodes, {G3.number_of_edges()} edges ({len(direct_edges_3)} direct, {len(indirect_edges_3)} indirect)")
