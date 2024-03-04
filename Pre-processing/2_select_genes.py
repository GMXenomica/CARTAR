## Keep only genes of interest (e.g: located in cell surface)

# Remove duplicated lines from GO file
GO = open('../Data/Raw/GO_0005886.txt','r') # Modify to GO list of interest
out = open('../Data/Processed/GO_simplified_list.txt','w')
lines = []
for line in GO:
	if line not in lines:
	    lines.append(line)
for line in lines:
	out.write(line)
out.close()

# Select genes located in the GO list
GO = open('../Data/Processed/GO_simplified_list.txt','r')
genes = []
for line in GO: 
	line = line.strip()
	gene = line.split()[0]
	genes.append(gene)

# GTEX
# remove genes that are not in our GO list
gtex = open('../Data/Processed/gtexTpm.csv','r')
no_membrane_genes = open('../Data/Processed/no_membrane_genes.csv','w')
exclusion = ''
out_gtex = open('../Data/Processed/gtexTpm_selected_v1.csv','w')
n = 0
for line in gtex: 
	# write header
	if n == 0:
		out_gtex.write(line)
		n = 1
	# selected genes information
	else: 
		gene = line.split(',')[0]
		if gene in genes:
			out_gtex.write(line)
		else:
			exclusion += gene + ','
no_membrane_genes.write(exclusion.strip(','))
out_gtex.close()
no_membrane_genes.close()

# TCGA
# remove genes that are not in our GO list
tcga = open('../Data/Processed/tcgaTpm.csv','r')
out_tcga = open('../Data/Processed/tcgaTpm_selected_v1.csv','w')
n = 0
for line in tcga: 
	# write header
	if n == 0:
		out_tcga.write(line)
		n = 1
	# selected genes information
	else: 
		gene = line.split(',')[0]
		if gene in genes:
			out_tcga.write(line)
out_tcga.close()