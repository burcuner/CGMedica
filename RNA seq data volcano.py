import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

# Separate Sham and TBI samples
sham_cols = [col for col in rna_data.columns if "Sham" in col]
tbi_cols = [col for col in rna_data.columns if "TBI" in col]

# Calculate log2 fold change and p-values
log2_fc = np.log2(rna_data[tbi_cols].mean(axis=1) / rna_data[sham_cols].mean(axis=1))
p_values = [ttest_ind(rna_data.loc[gene, sham_cols], rna_data.loc[gene, tbi_cols]).pvalue for gene in rna_data.index]

# Adjust for plotting
volcano_df = pd.DataFrame({
    'Gene': rna_data.index,
    'log2FC': log2_fc,
    '-log10(p-value)': -np.log10(p_values)
})
volcano_df['Significant'] = (volcano_df['log2FC'].abs() > 1) & (volcano_df['-log10(p-value)'] > -np.log10(0.05))

# Plot
plt.figure(figsize=(10, 6))
sns.scatterplot(data=volcano_df, x='log2FC', y='-log10(p-value)', hue='Significant', palette={True: 'red', False: 'grey'})
plt.axhline(-np.log10(0.05), color='blue', linestyle='--', linewidth=1)
plt.axvline(1, color='green', linestyle='--', linewidth=1)
plt.axvline(-1, color='green', linestyle='--', linewidth=1)
plt.title('Volcano Plot: Simulated RNA-seq TBI vs. Sham')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 p-value')
plt.legend(title='Significant', loc='upper right')
plt.tight_layout()
plt.grid(True)
plt.show()
