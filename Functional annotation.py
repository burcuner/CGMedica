Python 3.13.2 (tags/v3.13.2:4f8bb39, Feb  4 2025, 15:23:48) [MSC v.1942 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license()" for more information.
>>> # Fallback to manual annotation: generate mock pathway results for demonstration
... mock_pathways = pd.DataFrame({
...     'Pathway': [
...         'Oxidative Stress Response', 'Neuroinflammation', 'Mitochondrial Dysfunction',
...         'Microglial Activation', 'Regulation of Apoptosis', 'Cytokine Signaling',
...         'Blood-Brain Barrier Disruption', 'Cellular Senescence'
...     ],
...     'Enrichment Source': ['GO:BP', 'GO:BP', 'KEGG', 'GO:BP', 'Reactome', 'Reactome', 'GO:BP', 'Reactome'],
...     'Adjusted P-value': [0.0003, 0.0012, 0.0025, 0.0041, 0.0067, 0.0083, 0.0124, 0.0189],
...     'Overlap Genes': [
...         ['Gene0278', 'Gene0381'], ['Gene0156', 'Gene0169'], ['Gene0039', 'Gene0365'],
...         ['Gene0073', 'Gene0289'], ['Gene0351'], ['Gene0240'], ['Gene0132'], ['Gene0401']
...     ]
... })
... 
... # Save mock results
... mock_pathway_path = "/mnt/data/mock_pathway_annotation_results.csv"
... mock_pathways.to_csv(mock_pathway_path, index=False)
... 
... import ace_tools as tools; tools.display_dataframe_to_user(name="Mock Functional Pathway Annotation", dataframe=mock_pathways)
... mock_pathway_path
