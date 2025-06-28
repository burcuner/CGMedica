# Filter top disrupted genes by log2FC and p-value
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
import string

# Define realistic gene names as in your actual RNA/proteomics dataset
realistic_gene_symbols = [
    'RHOB', 'ALDH3B1', 'ALOX12B', 'GFAP', 'SOD2', 'IL1B', 'BDNF', 'MAPK1', 'TNF', 'NFE2L2',
    'CXCL10', 'CD44', 'MMP9', 'NOS2', 'VEGFA', 'CASP3', 'TLR4', 'FOS', 'HMOX1', 'STAT3',
    'TGFB1', 'AIF1', 'PTGS2', 'P2RX7', 'APP', 'C1QA', 'CSF1', 'CX3CR1', 'HIF1A', 'RPS6KB1',
    'PLA2G4A', 'NFATC1', 'CREB1', 'NRF1', 'RELA', 'ATP1A1', 'CAMK2A', 'PARK7', 'GSK3B', 'AKT1',
    'TRPV4', 'GRIA1', 'CACNA1C', 'CYBB', 'PINK1', 'OPA1', 'MFN2', 'FIS1', 'DRP1', 'JUN'
]

# Generate test data using actual gene names
np.random.seed(42)
volcano_df = pd.DataFrame({
    'Gene': realistic_gene_symbols[:50],
    'log2FC': np.random.normal(loc=0, scale=2, size=50),
    '-log10(p-value)': np.random.uniform(0, 10, size=50)
})

# Ensure '-log10(p-value)' is numeric and doesn't contain NaNs
volcano_df['-log10(p-value)'] = pd.to_numeric(volcano_df['-log10(p-value)'], errors='coerce')
volcano_df = volcano_df.dropna(subset=['-log10(p-value)'])

# Recalculate p-value from -log10(p-value)
volcano_df['p-value'] = 10 ** (-volcano_df['-log10(p-value)'])

# Filter and sort top disrupted genes
top_disrupted_genes = volcano_df[
    (volcano_df['log2FC'].abs() > 1) & (volcano_df['p-value'] < 0.05)
].copy()

top_disrupted_genes = top_disrupted_genes.sort_values(['p-value', 'log2FC'], ascending=[True, False])
prioritized_genes = top_disrupted_genes.head(20).reset_index(drop=True)

# Display table using print
print("\nTop 20 Prioritized Disrupted Genes Post-TBI:\n")
print(prioritized_genes.to_string(index=False))

# Volcano Plot
plt.figure(figsize=(10, 6))
sns.scatterplot(
    data=volcano_df,
    x='log2FC',
    y='-log10(p-value)',
    hue=(volcano_df['p-value'] < 0.05) & (volcano_df['log2FC'].abs() > 1),
    palette={True: 'red', False: 'gray'}
)
plt.axhline(-np.log10(0.05), linestyle='--', color='blue', linewidth=1)
plt.axvline(1, linestyle='--', color='green', linewidth=1)
plt.axvline(-1, linestyle='--', color='green', linewidth=1)

# Label top significant genes
for _, row in prioritized_genes.iterrows():
    plt.text(row['log2FC'], row['-log10(p-value)'], row['Gene'], fontsize=8, ha='right')

plt.title('Volcano Plot of Simulated Gene Expression Data')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 p-value')
plt.grid(True)
plt.tight_layout()
plt.show()

# Add mock functional pathway annotation results using actual genes
mock_pathways = pd.DataFrame({
    'Pathway': [
        'Oxidative Stress Response', 'Neuroinflammation', 'Mitochondrial Dysfunction',
        'Microglial Activation', 'Regulation of Apoptosis', 'Cytokine Signaling',
        'Blood-Brain Barrier Disruption', 'Cellular Senescence'
    ],
    'Enrichment Source': ['GO:BP', 'GO:BP', 'KEGG', 'GO:BP', 'Reactome', 'Reactome', 'GO:BP', 'Reactome'],
    'Adjusted P-value': [0.0003, 0.0012, 0.0025, 0.0041, 0.0067, 0.0083, 0.0124, 0.0189],
    'Overlap Genes': [
        ['HMOX1', 'SOD2'], ['IL1B', 'TNF'], ['OPA1', 'MFN2'],
        ['CX3CR1', 'CSF1'], ['CASP3'], ['CXCL10'], ['GFAP'], ['CD44']
    ]
})

# Display mock functional annotation table
print("\nMock Functional Pathway Annotation:\n")
print(mock_pathways.to_string(index=False))
