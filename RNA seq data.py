# Generate random but realistic gene names for RNA-seq (simulate human/rat-like gene symbols)
import random

# Generate synthetic gene symbols (e.g., real-looking names)
realistic_gene_symbols = ['RHOB', 'ALDH3B1', 'ALOX12B', 'GFAP', 'SOD2', 'IL1B', 'BDNF', 'MAPK1', 'TNF', 'NFE2L2',
                          'CXCL10', 'CD44', 'MMP9', 'NOS2', 'VEGFA', 'CASP3', 'TLR4', 'FOS', 'HMOX1', 'STAT3',
                          'TGFB1', 'AIF1', 'PTGS2', 'P2RX7', 'APP', 'C1QA', 'CSF1', 'CX3CR1', 'HIF1A', 'RPS6KB1',
                          'PLA2G4A', 'NFATC1', 'CREB1', 'NRF1', 'RELA', 'ATP1A1', 'CAMK2A', 'PARK7', 'GSK3B', 'AKT1',
                          'TRPV4', 'GRIA1', 'CACNA1C', 'CYBB', 'PINK1', 'OPA1', 'MFN2', 'FIS1', 'DRP1']

# If we have fewer symbols than needed, expand the list with suffixes
while len(realistic_gene_symbols) < len(rna_data.index):
    new_gene = random.choice(realistic_gene_symbols) + random.choice(string.ascii_uppercase)
    if new_gene not in realistic_gene_symbols:
        realistic_gene_symbols.append(new_gene)

# Assign gene symbols to RNA-seq data
rna_data.index = realistic_gene_symbols[:len(rna_data)]

# Update proteomics GeneSymbol column using updated RNA-seq gene names
new_prot_map = dict(zip(prot_data.index, random.sample(list(rna_data.index), len(prot_data))))
prot_data["GeneSymbol"] = [new_prot_map[p] for p in prot_data.index]

# Recalculate and replot RNA-seq volcano plot with gene names
log2_fc = np.log2(rna_data[tbi_cols].mean(axis=1) / rna_data[sham_cols].mean(axis=1))
p_values = [ttest_ind(rna_data.loc[gene, sham_cols], rna_data.loc[gene, tbi_cols]).pvalue for gene in rna_data.index]

volcano_df = pd.DataFrame({
    'Gene': rna_data.index,
    'log2FC': log2_fc,
    '-log10(p-value)': -np.log10(p_values)
})
volcano_df['Significant'] = (volcano_df['log2FC'].abs() > 1) & (volcano_df['-log10(p-value)'] > -np.log10(0.05))
top_genes = volcano_df.sort_values('-log10(p-value)', ascending=False).head(10)

# RNA-seq Volcano plot
plt.figure(figsize=(10, 6))
sns.scatterplot(data=volcano_df, x='log2FC', y='-log10(p-value)', hue='Significant', palette={True: 'red', False: 'grey'})
plt.axhline(-np.log10(0.05), color='blue', linestyle='--', linewidth=1)
plt.axvline(1, color='green', linestyle='--', linewidth=1)
plt.axvline(-1, color='green', linestyle='--', linewidth=1)
for _, row in top_genes.iterrows():
    plt.text(row['log2FC'], row['-log10(p-value)'], row['Gene'], fontsize=8, ha='right')
plt.title('Volcano Plot: Simulated RNA-seq (Gene Symbols Labeled)')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 p-value')
plt.legend(title='Significant', loc='upper right')
plt.grid(True)
plt.tight_layout()
plt.show()

# Update proteomics index to match gene-like symbols
prot_expression.index = prot_data["GeneSymbol"]

# Recalculate proteomics volcano
log2_fc_prot = np.log2(prot_expression[prot_tbi].mean(axis=1) / prot_expression[prot_sham].mean(axis=1))
p_values_prot = [
    ttest_ind(prot_expression.loc[gene, prot_sham], prot_expression.loc[gene, prot_tbi]).pvalue
    for gene in prot_expression.index
]

volcano_prot_df = pd.DataFrame({
    'Protein': prot_expression.index,
    'log2FC': log2_fc_prot,
    '-log10(p-value)': -np.log10(p_values_prot)
})
volcano_prot_df['Significant'] = (volcano_prot_df['log2FC'].abs() > 1) & (volcano_prot_df['-log10(p-value)'] > -np.log10(0.05))
top_proteins = volcano_prot_df.sort_values('-log10(p-value)', ascending=False).head(10)

# Proteomics Volcano plot
plt.figure(figsize=(10, 6))
sns.scatterplot(data=volcano_prot_df, x='log2FC', y='-log10(p-value)', hue='Significant', palette={True: 'purple', False: 'grey'})
plt.axhline(-np.log10(0.05), color='blue', linestyle='--', linewidth=1)
plt.axvline(1, color='green', linestyle='--', linewidth=1)
plt.axvline(-1, color='green', linestyle='--', linewidth=1)
for _, row in top_proteins.iterrows():
    plt.text(row['log2FC'], row['-log10(p-value)'], row['Protein'], fontsize=8, ha='right')
plt.title('Volcano Plot: Simulated Proteomics (Gene-like Symbols Labeled)')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 p-value')
plt.legend(title='Significant', loc='upper right')
plt.grid(True)
plt.tight_layout()
plt.show()
