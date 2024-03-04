## Combine all TCGA samples belonging to the same tumor group
## Combine all GTEX samples belonging to the same tissue
## Add GTEX samples to the corresponding control group of TCGA tumors

import pandas as pd

#TCGA
# Create a dictionary containing the list of samples belonging to the same tumor group
tcga = pd.read_csv('../Data/Processed/tcgaTpm_selected_v3.csv')
tumors = {}
n = 0 
for column in tcga.columns:
    if n == 0: 
        n = 1
    else: 
        stype = column.split('.')[0]
        fields = stype.split('_')
        tumor = fields[0]
        group = fields[1]
        if 'Primary' in group or 'Recurrent' in group:
            tag = ''
        elif group == 'Metastatic': 
            tag = '_Metastatic'
        elif 'Normal' in group:
            tag = '_Normal'
        sample_group = tumor + tag
        if sample_group not in tumors.keys():
            tumors[sample_group] = [column]
        else:
            tumors[sample_group].append(column)

#GTEX
# Create a dictionary containing the list of samples belonging to the same tissue
gtex = pd.read_csv('../Data/Processed/gtexTpm_selected_v3.csv')
tissues = {}
n = 0 
for column in gtex.columns:
    if n == 0: 
        n = 1
    else: 
        if '.' not in column:
            tissues[column] = [column]   
        else: 
           tissue = column.split('.')[0]
           tissues[tissue].append(column) 
# Calculate median expression for the gtex data
result = pd.DataFrame()
result['gene'] = gtex['gene']
out = open('../Data/Processed/replicates.tsv','w') # Save in a txt file the number of replicates of each tissue/tumor
for tissue in tissues.keys():
    replicate = len(tissues[tissue])
    out.write(tissue + '\t' + str(replicate) + '\n')
    median_per_tissue = gtex[tissues[tissue]].median(axis=1)
    result[str(tissue)] = median_per_tissue
result.to_csv('../Data/Processed/gtex_targetable_gene_Tpm_by_tissue.csv', index=False)
for tumor in tumors.keys():
    replicate = len(tumors[tumor])
    out.write(tumor + '\t' + str(replicate) + '\n')
out.close()

# Create a dictionary with all the gtex tissues added to control samples
gtex_tcga = {'ACC':'Adrenal Gland','BLCA':'Bladder','BRCA':'Breast','CESC':'Cervix Uteri', 'COAD':'Colon', 'DLBC':'Blood',
                'ESCA':'Esophagus','GBM':'Brain','KICH':'Kidney','KIRC':'Kidney','KIRP':'Kidney','LAML':'Bone Marrow','LGG':'Brain',
                'LIHC':'Liver','LUAD':'Lung','LUSC':'Lung','OV':'Ovary','PAAD':'Pancreas','PRAD':'Prostate','READ':'Colon','SKCM':'Skin',
                'STAD':'Stomach','TGCT':'Testis','THCA':'Thyroid','THYM':'Blood','UCEC':'Uterus','UCS':'Uterus'}

#TCGA-GTEX
# Combine the columns of tcga and gtex data
combine = pd.merge(tcga, gtex, on='gene')
combine.to_csv('../Data/Processed/targetable_genes_gtex_tcga.csv', index=False) 

# Add the GTEX samples that act as control samples of any tumor
# Add GTEX samples to the dicitinary in the tumor for which they can be used as control
for tumor in gtex_tcga.keys():
    control = gtex_tcga[tumor]
    for tissue in tissues.keys():
        if tissue == control and tissue not in tcga.keys():
            for sample in tissues[tissue]:    
                tcga[sample] = gtex[sample]
tcga.to_csv('../Data/Processed/tcga_gtex_combined_data.csv', index=False)
# Add the GTEX control samples to the corresponding group in the dictionary
for control in gtex_tcga.keys():
    tissue = gtex_tcga[control]
    group = control + '_Normal'
    if group in tumors.keys():
        tumors[group].extend(tissues[tissue])
    else:
        tumors[group] = tissues[tissue]
# Calculate median expression of each group 
combined = pd.read_csv('../Data/Processed/tcga_gtex_combined_data.csv')
out = open('../Data/Processed/group_replicates.tsv','w') # Save in a txt file the number of replicates of each group
result = pd.DataFrame()
for tumor in tumors.keys():
    replicate = len(tumors[tumor])
    out.write(tumor + '\t' + str(replicate) + '\n')
    median_per_group = combined[tumors[tumor]].median(axis=1)
    result[str(tumor)] = median_per_group
out.close()
result_final = pd.DataFrame()
sorted_columns = sorted(result.columns)
result_final = result[sorted_columns]
result_final.insert(0, 'gene', combined['gene'])
result_final.to_csv('../Data/Processed/targetable_gene_Tpm_TumorVsControl.csv', index=False)
# Change negative value to expression data to minimum expression (1E-08)
data = open('../Data/Processed/targetable_gene_Tpm_TumorVsControl.csv','r')
out = open('../Data/Processed/targetable_gene_Tpm_TumorVsControl_final.csv','w')
n = 0
for line in data:
    if n == 0:
        out.write(line)
        n = 1
    else:
        line = line.strip()
        fields = line.split(',')
        write = fields[0] + ','
        for field in fields[1:]:
           if '-' in field:
               write += '1.0e-08,'
           else:
               write += field + ','
        write = write.strip(',')
        out.write(write + '\n')
out.close()      