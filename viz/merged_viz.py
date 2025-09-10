import matplotlib.pyplot as plt
import numpy as np

# Set style for publication quality
plt.style.use('default')

# Create figure with 2 panels
fig = plt.figure(figsize=(12, 6))

# Panel A: Network Size Comparison
ax1 = plt.subplot(1, 2, 1)
networks = ['Ferroptosis', 'Autophagy']
protein_counts = [1880, 16470]
colors = ['#E74C3C', '#2ECC71']  # Red and green for good contrast

bars = ax1.bar(networks, protein_counts, color=colors, alpha=0.8, edgecolor='black', linewidth=1)
for bar, count in zip(bars, protein_counts):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
             f'{count:,}', ha='center', va='bottom', fontweight='bold', fontsize=12)

ax1.set_ylabel('Number of Proteins', fontsize=12, fontweight='bold')
ax1.set_title('A. Network Size Comparison', fontsize=14, fontweight='bold', pad=20)
ax1.set_ylim(0, max(protein_counts) * 1.15)
ax1.grid(axis='y', alpha=0.3)

# Panel B: Final Database Edge Composition
ax2 = plt.subplot(1, 2, 2)
edge_sources = ['ARN\n(Autophagy)', 'Ferroptosis\nNetworks', 'Cross-Network\nConnections']
edge_counts = [37984, 8990, 3763]
colors = ['#3498DB', '#F39C12', '#9B59B6']  # Blue, orange, purple for good distinction

# Create horizontal bar chart for better label visibility
bars = ax2.barh(edge_sources, edge_counts, color=colors, alpha=0.8, edgecolor='black', linewidth=1)
for bar, count in zip(bars, edge_counts):
    width = bar.get_width()
    ax2.text(width + width*0.01, bar.get_y() + bar.get_height()/2.,
             f'{count:,}', ha='left', va='center', fontweight='bold', fontsize=12)

ax2.set_xlabel('Number of Edges', fontsize=12, fontweight='bold')
ax2.set_title('B. Final Database Edge Composition', fontsize=14, fontweight='bold', pad=20)
ax2.set_xlim(0, max(edge_counts) * 1.2)
ax2.grid(axis='x', alpha=0.3)

# Add overall title and adjust layout
fig.suptitle('Ferroptosis-Autophagy Network Integration Analysis',
             fontsize=16, fontweight='bold', y=0.95)

plt.tight_layout()
plt.subplots_adjust(top=0.85, wspace=0.3)
plt.show()

# Print summary statistics for manuscript
print("=== Summary Statistics for Manuscript ===")
print(f"Total proteins in integrated database: {1880 + 16470:,}")
print(f"Cross-network interactions: 3,763")
print(f"Cross-network density: {3763 / (1880 * 16470) * 100:.4f}% of possible connections")
print(f"Total edges in final database: {37984 + 8990 + 3763:,}")
print(f"Cross-network edges represent {3763 / (37984 + 8990 + 3763) * 100:.1f}% of all edges")
