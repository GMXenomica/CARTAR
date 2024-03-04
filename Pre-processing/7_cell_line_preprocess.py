# Preporcess cell line data files
import pandas as pd

# Open files
expression = pd.read_csv('../Data/Raw/Expression_Public_23Q4.csv')
metadata = open('../Data/Raw/cell_line_metadata.csv', 'r')
out = open('../Data/Processed/cell_line_metadata.csv','w')
out2 = open('../Data/Processed/genes_cells.csv','w')

# Remove cell lines from metadata file not pressent in expression file
# Identify cell lines in expression file
exp_lines = expression['Unnamed: 0'].tolist()
# Create metadata file only with those cell lines
firstline = 1
for line in metadata:
    if firstline:
        out.write(line)
        firstline = 0
    else:
        cell_line = line.split(',')[0]
        if cell_line in exp_lines:
            out.write(line)
out.close()
df = pd.read_csv('../Data/Processed/cell_line_metadata.csv')
columns = ['ModelID','CellLineName','OncotreeLineage','OncotreePrimaryDisease','OncotreeSubtype','OncotreeCode','CatalogNumber']
df_reduced = df[columns]
df_reduced.to_csv('../Data/Processed/cell_line_metadata_reduced.csv', index=False)

# Keep only expression data of genes in the indicated Gene Ontology (e.g: located in cell surface)
GO = open('../Data/Processed/GO_simplified_list.txt','r')
genes = ['Unnamed: 0']
for line in GO: 
	line = line.strip()
	gene = line.split()[0]
	genes.append(gene)
firstline = 1
found_genes = []
for gene in genes: 
	if gene in expression.columns:
		found_genes.append(gene)
write = ','.join(found_genes[1:])
out2.write(write)
out2.close()
expression_filtered = expression[found_genes]
expression_filtered.to_csv('../Data/Processed/Expression_Public_23Q4_filtered.csv', index=False)