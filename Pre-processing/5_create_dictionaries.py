import pickle
import statistics
from scipy.stats import mannwhitneyu

# Create dicitionaries with expression values of all genes for all tumoral groups
# List with all tumor abbreviations
tumors = ['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS']

# Read expression data for the samples
data = open(r'C:\Users\34656\Pipeline for TCGA_GTEX\Data\Processed\tcga_gtex_combined_data.csv','r')
tumor_replicates = {} # Dictionary with the following strcutre {tumors_group:[position_indexes]}
result1 = {} # Dictionary with the following strcutre {genes:tumors:{groups:[TPM_expression_values]}} A-C
result2 = {} #D-J
result3 = {} #K-N
result4 = {} #O-R
result5 = {} #S-T
result6 = {} #U-Z
firstline = 1 
n = 0
# Identifiy the position of each tumor_group in the lines
for line in data:
    fields = line.strip().split(',')
    if firstline:
        for field in fields[1:]:
            n += 1
            field = field.split('.')[0]
            abr = field.split('_')[0]
            if abr in tumors:
                if 'Primary' in field:
                    name = abr + '_Tumor'
                    if name not in tumor_replicates.keys():
                        tumor_replicates[name] = [n]
                    else:
                        tumor_replicates[name].append(n)
                elif 'Normal' in field:
                    name = abr + '_Normal'
                    if name not in tumor_replicates.keys():
                        tumor_replicates[name] = [n]
                    else:
                        tumor_replicates[name].append(n)
        firstline = 0
    # Idnetify the expression values of the introduced gene belonging to each tumor_group
    else: 
        gene = fields[0]
        if 'A' <= gene[0] <= 'C':
            result1[gene] = {}
            for sample in tumor_replicates.keys():
                name = sample.split('_')
                abr = name[0]
                group = name[1]
                if abr not in result1[gene].keys():
                    result1[gene][abr] = {}
                result1[gene][abr][group] = []
                for index in tumor_replicates[sample]:
                    result1[gene][abr][group].append(float(fields[index]))
        elif 'D' <= gene[0] <= 'J':
            result2[gene] = {}
            for sample in tumor_replicates.keys():
                name = sample.split('_')
                abr = name[0]
                group = name[1]
                if abr not in result2[gene].keys():
                    result2[gene][abr] = {}
                result2[gene][abr][group] = []
                for index in tumor_replicates[sample]:
                    result2[gene][abr][group].append(float(fields[index]))
        elif 'K' <= gene[0] <= 'N':
            result3[gene] = {}
            for sample in tumor_replicates.keys():
                name = sample.split('_')
                abr = name[0]
                group = name[1]
                if abr not in result3[gene].keys():
                    result3[gene][abr] = {}
                result3[gene][abr][group] = []
                for index in tumor_replicates[sample]:
                    result3[gene][abr][group].append(float(fields[index]))
        elif 'O' <= gene[0] <= 'R':
            result4[gene] = {}
            for sample in tumor_replicates.keys():
                name = sample.split('_')
                abr = name[0]
                group = name[1]
                if abr not in result4[gene].keys():
                    result4[gene][abr] = {}
                result4[gene][abr][group] = []
                for index in tumor_replicates[sample]:
                    result4[gene][abr][group].append(float(fields[index]))
        elif 'S' <= gene[0] <= 'T':
            result5[gene] = {}
            for sample in tumor_replicates.keys():
                name = sample.split('_')
                abr = name[0]
                group = name[1]
                if abr not in result5[gene].keys():
                    result5[gene][abr] = {}
                result5[gene][abr][group] = []
                for index in tumor_replicates[sample]:
                    result5[gene][abr][group].append(float(fields[index]))
        elif 'U' <= gene[0] <= 'Z':
            result6[gene] = {}
            for sample in tumor_replicates.keys():
                name = sample.split('_')
                abr = name[0]
                group = name[1]
                if abr not in result6[gene].keys():
                    result6[gene][abr] = {}
                result6[gene][abr][group] = []
                for index in tumor_replicates[sample]:
                    result6[gene][abr][group].append(float(fields[index]))
# Order tumors alphabetically
for gene, tumors in result1.items():
    result1[gene] = {tumor: tumors[tumor] for tumor in sorted(tumors)}
for gene, tumors in result2.items():
    result2[gene] = {tumor: tumors[tumor] for tumor in sorted(tumors)}
for gene, tumors in result3.items():
    result3[gene] = {tumor: tumors[tumor] for tumor in sorted(tumors)}
for gene, tumors in result4.items():
    result4[gene] = {tumor: tumors[tumor] for tumor in sorted(tumors)}
for gene, tumors in result5.items():
    result5[gene] = {tumor: tumors[tumor] for tumor in sorted(tumors)}
for gene, tumors in result6.items():
    result6[gene] = {tumor: tumors[tumor] for tumor in sorted(tumors)}
#Save dictionary in file
with open('../Data/Processed/tcga_AC.pkl', 'wb') as archivo:
    pickle.dump(result1, archivo)
with open('../Data/Processed/tcga_DJ.pkl', 'wb') as archivo:
    pickle.dump(result2, archivo)
with open('../Data/Processed/tcga_KN.pkl', 'wb') as archivo:
    pickle.dump(result3, archivo)
with open('../Data/Processed/tcga_OR.pkl', 'wb') as archivo:
    pickle.dump(result4, archivo)
with open('../Data/Processed/tcga_ST.pkl', 'wb') as archivo:
    pickle.dump(result5, archivo)
with open('../Data/Processed/tcga_UZ.pkl', 'wb') as archivo:
    pickle.dump(result6, archivo)

# Create dicitionaries with expression values of all genes for all GTEX tissues  
# List with all gtex tissues
tissues = ['Blood','Blood Vessel','Brain','Thyroid','Pancreas','Muscle','Lung','Skin','Colon','Nerve','Adipose Tissue','Ovary','Heart','Breast','Pituitary','Testis','Vagina','Esophagus','Small Intestine','Spleen','Adrenal Gland','Stomach','Uterus','Liver','Bone Marrow','Salivary Gland','Prostate','Kidney','Bladder','Fallopian Tube','Cervix Uteri']
# Function for statistical analysis of significance of tumoral group vs gtex data and customize plot function
# Read expression data for the samples
data = open(r'C:\Users\34656\Pipeline for TCGA_GTEX\Data\Processed\targetable_genes_gtex_tcga.csv','r')
tissue_replicates = {} # Dictionary with the following strcutre {tissue:[position_indexes]}
result1 = {} # Dictionary with the following strcutre {genes:{tissue:[TPM_expression_values]}} A-C
result2 = {} #D-J
result3 = {} #K-N
result4 = {} #O-R
result5 = {} #S-T
result6 = {} #U-Z
firstline = 1
n = 0
for line in data:
    fields = line.strip().split(',')
    if firstline:
        for field in fields[1:]:
            n += 1
            tissue = field.split('.')[0]
            # Identify the tissue present in each field
            if tissue in tissues:
                if tissue not in tissue_replicates.keys():
                    tissue_replicates[tissue] = [n]
                else:
                    tissue_replicates[tissue].append(n)
        firstline = 0
    else: 
        gene = fields[0]
        if 'A' <= gene[0] <= 'C':
            result1[gene] = {}
            for tissue in tissue_replicates.keys():
                if tissue not in result1[gene].keys():
                    result1[gene][tissue] = []
                for index in tissue_replicates[tissue]:
                    result1[gene][tissue].append(float(fields[index]))
        elif 'D' <= gene[0] <= 'J':
            result2[gene] = {}
            for tissue in tissue_replicates.keys():
                if tissue not in result2[gene].keys():
                    result2[gene][tissue] = []
                for index in tissue_replicates[tissue]:
                    result2[gene][tissue].append(float(fields[index]))
        elif 'K' <= gene[0] <= 'N':
            result3[gene] = {}
            for tissue in tissue_replicates.keys():
                if tissue not in result3[gene].keys():
                    result3[gene][tissue] = []
                for index in tissue_replicates[tissue]:
                    result3[gene][tissue].append(float(fields[index]))
        elif 'O' <= gene[0] <= 'R':
            result4[gene] = {}
            for tissue in tissue_replicates.keys():
                if tissue not in result4[gene].keys():
                    result4[gene][tissue] = []
                for index in tissue_replicates[tissue]:
                    result4[gene][tissue].append(float(fields[index]))
        elif 'S' <= gene[0] <= 'T':
            result5[gene] = {}
            for tissue in tissue_replicates.keys():
                if tissue not in result5[gene].keys():
                    result5[gene][tissue] = []
                for index in tissue_replicates[tissue]:
                    result5[gene][tissue].append(float(fields[index]))
        elif 'U' <= gene[0] <= 'Z':
            result6[gene] = {}
            for tissue in tissue_replicates.keys():
                if tissue not in result6[gene].keys():
                    result6[gene][tissue] = []
                for index in tissue_replicates[tissue]:
                    result6[gene][tissue].append(float(fields[index]))
# Order tissues alphabetically
for gene, tissues in result1.items():
    result1[gene] = {tissue: tissues[tissue] for tissue in sorted(tissues)}
for gene, tissues in result2.items():
    result2[gene] = {tissue: tissues[tissue] for tissue in sorted(tissues)}
for gene, tissues in result3.items():
    result3[gene] = {tissue: tissues[tissue] for tissue in sorted(tissues)}
for gene, tissues in result4.items():
    result4[gene] = {tissue: tissues[tissue] for tissue in sorted(tissues)}
for gene, tissues in result5.items():
    result5[gene] = {tissue: tissues[tissue] for tissue in sorted(tissues)}
for gene, tissues in result6.items():
    result6[gene] = {tissue: tissues[tissue] for tissue in sorted(tissues)}
# Save dictionary in file
with open('../Data/Processed/gtex_AC.pkl', 'wb') as archivo:
    pickle.dump(result1, archivo)
with open('../Data/Processed/gtex_DJ.pkl', 'wb') as archivo:
    pickle.dump(result2, archivo)
with open('../Data/Processed/gtex_KN.pkl', 'wb') as archivo:
    pickle.dump(result3, archivo)
with open('../Data/Processed/gtex_OR.pkl', 'wb') as archivo:
    pickle.dump(result4, archivo)
with open('../Data/Processed/gtex_ST.pkl', 'wb') as archivo:
    pickle.dump(result5, archivo)
with open('../Data/Processed/gtex_UZ.pkl', 'wb') as archivo:
    pickle.dump(result6, archivo)

# Create SKCM metastatic dicitionary
data = open(r'c:\Users\34656\Pipeline for TCGA_GTEX\Data\Processed\tcga_gtex_combined_data.csv','r')
gtex_tcga = {'SKCM':'Skin'}
tissue_replicates = {} # Dictionary identifying the field in which is group is found ({group:[indexes]})
result = {} # Dictionary identifying the expression values of each group (gene{tumor:{group:[TPM_expression_values]}})
firstline = 1
n = 0
for line in data:
    fields = line.strip().split(',')
    # Select samples belonging to SKCM and control group
    if firstline:
        for field in fields[1:]:
            n += 1
            field = field.split('.')[0]
            abr = field.split('_')[0]
            # Identify the group
            if abr == 'SKCM':
                if 'Primary' in field:
                    name = 'Primary'
                    if name not in tissue_replicates.keys():
                        tissue_replicates[name] = [n]
                    else:
                        tissue_replicates[name].append(n)
                elif 'Normal' in field:
                    name = 'Normal'
                    if name not in tissue_replicates.keys():
                        tissue_replicates[name] = [n]
                    else:
                        tissue_replicates[name].append(n)
                elif 'Metastatic' in field:
                    name = 'Metastatic'
                    if name not in tissue_replicates.keys():
                        tissue_replicates[name] = [n]
                    else:
                        tissue_replicates[name].append(n)
            # Identify control samples from GTEX
            elif field in gtex_tcga.values():   
                for abr, tissue in gtex_tcga.items():
                    if field == tissue and abr == 'SKCM':
                        name = 'Normal'
                        if name not in tissue_replicates.keys():
                            tissue_replicates[name] = [n]
                        else:
                            tissue_replicates[name].append(n)
        firstline = 0
    else: 
        # Select expresion values of SKCM according to the group
        gene = fields[0]
        result[gene] = {}
        result[gene]['SKCM'] = {}
        for group in tissue_replicates.keys():
            result[gene]['SKCM'][group] = []
            for index in tissue_replicates[group]:
                result[gene]['SKCM'][group].append(float(fields[index]))
with open('../Data/Processed/SKCM.pkl', 'wb') as archivo:
    pickle.dump(result, archivo)

# Create a dictionary with median TPM expression values for all genes and groups 
tumors = ['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS']
gtex_tcga = {'ACC':'Adrenal Gland','BLCA':'Bladder','BRCA':'Breast','CESC':'Cervix Uteri', 'COAD':'Colon', 'DLBC':'Blood',
                'ESCA':'Esophagus','GBM':'Brain','KICH':'Kidney','KIRC':'Kidney','KIRP':'Kidney','LAML':'Bone Marrow','LGG':'Brain',
                'LIHC':'Liver','LUAD':'Lung','LUSC':'Lung','OV':'Ovary','PAAD':'Pancreas','PRAD':'Prostate','READ':'Colon','SKCM':'Skin',
                'STAD':'Stomach','TGCT':'Testis','THCA':'Thyroid','THYM':'Blood','UCEC':'Uterus','UCS':'Uterus'}

# Read expression data for the samples
data = open(r'C:\Users\34656\Pipeline for TCGA_GTEX\Data\Processed\tcga_gtex_combined_data.csv','r')
tumor_replicates = {} # Dictionary with the following strcutre {tumors_group:[position_indexes]}
result = {} # Ditionary with following structure {gene:{tumor:{group:[median,sample_size]}}}
firstline = 1
n = 0
# Identifiy the position of each tumor_group in the lines
for line in data:
    fields = line.strip().split(',')
    if firstline:
        for field in fields[1:]:
            n += 1
            field = field.split('.')[0]
            abr = field.split('_')[0]
            if abr in tumors:
                if 'Primary' in field:
                    name = abr + '_Tumor'
                    if name not in tumor_replicates.keys():
                        tumor_replicates[name] = [n]
                    else:
                        tumor_replicates[name].append(n)
                elif 'Normal' in field:
                    name = abr + '_Normal'
                    if name not in tumor_replicates.keys():
                        tumor_replicates[name] = [n]
                    else:
                        tumor_replicates[name].append(n)
            else:                     
                for tumor, tissue in gtex_tcga.items():
                    if abr == tissue:
                        name = tumor + '_Normal'
                        if name not in tumor_replicates.keys():
                            tumor_replicates[name] = [n]
                        else:
                            tumor_replicates[name].append(n)
        firstline = 0
    # Idnetify the expression values of the introduced gene belonging to each tumor_group
    else: 
        gene = fields[0]
        result[gene] = {}
        for sample in tumor_replicates.keys():
            values = []
            name = sample.split('_')
            abr = name[0]
            group = name[1]
            if abr not in result[gene].keys():
                result[gene][abr] = {}
            for index in tumor_replicates[sample]:
                values.append(float(fields[index]))
            median = statistics.median(values)
            result[gene][abr][group] = []
            result[gene][abr][group].append(median)
            result[gene][abr][group].append(len(values))
with open('../Data/Processed/median.pkl', 'wb') as archivo:
    pickle.dump(result, archivo)   

# Create dictionary with tumor median and p-values comparing 'Primary tumor' and 'Normal' samples

# Identifiy the position of each tumor_group in the lines
for line in data:
    fields = line.strip().split(',')
    if firstline:
        firstline = 0
    # Idnetify the expression values of the introduced gene belonging to each tumor_group
    else: 
        gene = fields[0]
        result[gene] = {}
        for tumor in tumors:
            T = tumor + '_Tumor'
            N = tumor + '_Normal'
            N_values = []
            T_values = []
            for index in tumor_replicates[T]:
                T_values.append(float(fields[index]))
            for index in tumor_replicates[N]:
                N_values.append(float(fields[index]))
            _, p_value = mannwhitneyu(T_values, N_values)
            result[gene][tumor] = p_value
with open('../Data/Processed/p_value.pkl', 'wb') as archivo:
    pickle.dump(result, archivo)  