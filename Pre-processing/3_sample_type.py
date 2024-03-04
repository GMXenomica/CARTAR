## Identify to which tumor/tissue belongs each sample and remove those without this information

import pandas as pd

# Add cancer type abbreviation to TCGA_phenotype file
data = pd.read_csv('../Data/Raw/TCGA_phenotype_denseDataOnlyDownload.tsv',sep='\t')
disease_match = {'kidney chromophobe':'KICH', 'colon adenocarcinoma':'COAD', 
                 'lung squamous cell carcinoma':'LUSC', 'bladder urothelial carcinoma':'BLCA', 
                 'diffuse large B-cell lymphoma':'DLBC', 'adrenocortical cancer':'ACC',
                 'cholangiocarcinoma':'CHOL', 'kidney papillary cell carcinoma':'KIRP', 
                 'esophageal carcinoma':'ESCA', 'liver hepatocellular carcinoma':'LIHC', 
                 'uterine carcinosarcoma':'UCS', 'glioblastoma multiforme':'GBM', 
                 'stomach adenocarcinoma':'STAD', 'uveal melanoma':'UVM', 
                 'pancreatic adenocarcinoma':'PAAD', 'thyroid carcinoma':'THCA', 
                 'sarcoma':'SARC', 'uterine corpus endometrioid carcinoma':'UCEC', 
                 'thymoma':'THYM', 'kidney clear cell carcinoma':'KIRC', 
                 'acute myeloid leukemia':'LAML', 'skin cutaneous melanoma':'SKCM', 
                 'lung adenocarcinoma':'LUAD', 'mesothelioma':'MESO', 
                 'brain lower grade glioma':'LGG', 'breast invasive carcinoma':'BRCA', 
                 'testicular germ cell tumor':'TGCT', 'rectum adenocarcinoma':'READ', 
                 'head & neck squamous cell carcinoma':'HNSC', 
                 'cervical & endocervical cancer':'CESC', 
                 'pheochromocytoma & paraganglioma':'PCPG', 
                 'prostate adenocarcinoma':'PRAD', 
                 'ovarian serous cystadenocarcinoma':'OV'}
primary_disease = []
for each in data['_primary_disease'].tolist():
    primary_disease.append(disease_match[each])
data['primary_disease'] = primary_disease
data.to_csv('../Data/Processed/TCGA_sample_cancertype.tsv', sep='\t', index=False)

# TCGA
# Create a dictionay assigning each sample to the tumor type and remove samples annotated as 'Addiotional'
phenotype = open('../Data/Processed/TCGA_sample_cancertype.tsv','r')
samples_cancer = {}
n = 0 
for line in phenotype: 
    line = line.strip()
    if n == 0:
        n = 1
    else:
        fields = line.split('\t')
        stype = fields [2]
        if stype and 'Additional' not in stype:
            sample = fields[0]
            cancer = fields[-1]
            sample_type = cancer + '_' + stype
            samples_cancer[sample] = sample_type
# Remove samples with no tumor type information
tcga_matrix = pd.read_csv('../Data/Processed/tcgaTpm_selected_v1.csv') 
remove = []
n = 0
for column in tcga_matrix.columns:
    if n == 0: 
        n = 1
    else: 
        if column not in samples_cancer.keys():
            remove.append(column)
for sample in remove: 
    del tcga_matrix[sample]
tcga_matrix.to_csv('../Data/Processed/tcgaTpm_selected_v2.csv', index=False)
# Replace the sample name by the tumor type
tcga = open('../Data/Processed/tcgaTpm_selected_v2.csv','r')
out = open('../Data/Processed/tcgaTpm_selected_v3.csv','w')
n = 0
for line in tcga: 
    if n == 0: 
        line = line.strip()
        samples = line.split(',')
        n = 1
        header = ''
        x = 0
        for sample in samples: 
            if x == 0:
                header += sample + ','
                x = 1
            else: 
                header += samples_cancer[sample] + ','
        header = header.strip(',')
        out.write(header + '\n')
    else: 
        out.write(line)
out.close()
       
# GTEX
# Create a dictionay assigning each sample to the tissue type
phenotype = open('../Data/Raw/GTEX_phenotype.tsv','r')
samples_tissue = {}
n = 0 
for line in phenotype: 
    line = line.strip()
    if n == 0:
        n = 1
    else:
        fields = line.split('\t')
        tissue = fields[2]
        if tissue != '<not provided>':
            sample = fields[0]
            samples_tissue[sample] = tissue
# Remove samples with no tissue type information
gtex_matrix = pd.read_csv('../Data/Processed/gtexTpm_selected_v1.csv')
remove = []
n = 0
for column in gtex_matrix.columns:
    if n == 0: 
        n = 1
    else: 
        if column not in samples_tissue.keys():
            remove.append(column)
for sample in remove: 
    del gtex_matrix[sample]
gtex_matrix.to_csv('../Data/Processed/gtexTpm_selected_v2.csv', index=False)
# Replace the sample name by the tissue type
gtex = open('../Data/Processed/gtexTpm_selected_v2.csv','r')
out = open('../Data/Processed/gtexTpm_selected_v3.csv','w')
n = 0
for line in gtex: 
    if n == 0: 
        line = line.strip()
        samples = line.split(',')
        n = 1
        header = ''
        x = 0
        for sample in samples: 
            if x == 0:
                header += sample + ','
                x = 1
            else: 
                header += samples_tissue[sample] + ','
        header = header.strip(',')
        out.write(header + '\n')
    else: 
        out.write(line)
out.close()