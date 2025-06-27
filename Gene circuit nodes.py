Python 3.13.2 (tags/v3.13.2:4f8bb39, Feb  4 2025, 15:23:48) [MSC v.1942 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license()" for more information.
>>> # Filter top disrupted genes by log2FC and p-value
... volcano_df['p-value'] = 10 ** (-volcano_df['-log10(p-value)'])
... top_disrupted_genes = volcano_df[
...     (volcano_df['log2FC'].abs() > 1) & (volcano_df['p-value'] < 0.05)
... ].copy()
... 
... # Sort by significance (lowest p-value) and highest fold change
... top_disrupted_genes = top_disrupted_genes.sort_values(['p-value', 'log2FC'], ascending=[True, False])
... 
... # Select top 20 for prioritization
... prioritized_genes = top_disrupted_genes.head(20).reset_index(drop=True)
... 
... import ace_tools as tools; tools.display_dataframe_to_user(name="Top 20 Prioritized Disrupted Genes Post-TBI", dataframe=prioritized_genes)
... prioritized_genes.head(10)
