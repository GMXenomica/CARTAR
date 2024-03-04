from math import log2
import pandas as pd

# Transform expression data to log2(TPM+1)
data = open('../Data/Processed/targetable_gene_Tpm_TumorVsControl.csv','r')
out = open('../Data/Processed/combined_log2_expression.csv','w')
n = 0
for line in data:
    if n == 0:
        out.write(line)
        n = 1
    else: 
        line = line.strip()
        fields = line.split(',')
        write = fields[0]
        for field in fields[1:]:
            field = float(field)
            expression = log2(field +1)
            write += ',' + str(expression)
        out.write(write + '\n')
out.close()

# Remove columns of Metastasic samples and tumors without control data
data = pd.read_csv('../Data/Processed/combined_log2_expression.csv', index_col= False)
out = pd.DataFrame()
for column in data.columns:
    if column != 'MESO' and column != 'UVM':
        if 'Metastatic' not in column:
            out[column] = data[column]
out.to_csv('../Data/Processed/combined_log2_expression_v2.csv', index=False) 

# Calculate log2(FC) as log2(TPM+1) tumor expression - log2(TPM+1) control expression 
data = open('../Data/Processed/combined_log2_expression_v2.csv','r')
out = open('../Data/Processed/log2FC_expression.csv','w')
n = 0
for line in data: 
    line = line.strip()
    fields = line.split(',')
    if n == 0:
        write = fields[0]
        for field in fields[1::2]:
            write += ',' + field 
        n = 1
        out.write(write + '\n')
    else: 
        write = fields[0]
        k = 2
        for field in fields[1::2]:
            logFC = float(field) - float(fields[k])
            k += 2
            write += ',' + str(logFC)
        out.write(write + '\n')
out.close()