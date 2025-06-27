# Prepare proteomics data with gene-like labels
prot_expression = prot_data.drop(columns="GeneSymbol")
prot_expression.index = prot_data["GeneSymbol"]

# Proteomics volcano values
log2_fc_prot = np.log2(prot_expression[tbi_cols].mean(axis=1) / prot_expression[sham_cols].mean(axis=1))
p_values_prot = [
    ttest_ind(prot_expression.loc[gene, sham_cols], prot_expression.loc[gene, tbi_cols]).pvalue
    for gene in prot_expression.index
]

volcano_prot_df = pd.DataFrame({
    'Protein': prot_expression.index,
    'log2FC': log2_fc_prot,
    '-log10(p-value)': -np.log10(p_values_prot)
})
volcano_prot_df['Significant'] = (volcano_prot_df['log2FC'].abs() > 1) & (volcano_prot_df['-log10(p-value)'] > -np.log10(0.05))
top_proteins = volcano_prot_df.sort_values('-log10(p-value)', ascending=False).head(10)

# Proteomics Volcano plot with labels
plt.figure(figsize=(10, 6))
sns.scatterplot(data=volcano_prot_df, x='log2FC', y='-log10(p-value)', hue='Significant', palette={True: 'purple', False: 'grey'})
plt.axhline(-np.log10(0.05), color='blue', linestyle='--', linewidth=1)
plt.axvline(1, color='green', linestyle='--', linewidth=1)
plt.axvline(-1, color='green', linestyle='--', linewidth=1)
for _, row in top_proteins.iterrows():
    plt.text(row['log2FC'], row['-log10(p-value)'], row['Protein'], fontsize=8, ha='right')
plt.title('Volcano Plot: Simulated Proteomics (Gene Symbols Labeled)')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 p-value')
plt.legend(title='Significant', loc='upper right')
plt.grid(True)
plt.tight_layout()
plt.show()
