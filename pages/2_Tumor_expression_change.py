from unittest import skip
import streamlit as st
import pandas as pd
from math import log2
import requests
import base64

st.set_page_config(page_title='CARTAR', page_icon='logo.png',layout='wide')
mystyle = '''
    <style>
        p {
            text-align: justify;
        }
    </style>
    '''

st.logo('logo_v2.png', icon_image='logo.png')

st.markdown(
    """
    <style>
    [data-testid="stElementToolbar"] {
        display: none;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# Create footer
def create_footer():
    footer_container = st.container()
    with footer_container:
        st.markdown("<br>" * 1, unsafe_allow_html=True)  # Añade espacio en blanco
        st.markdown("""
        <div style="background-color: #f0f2f6; padding: 10px; text-align: center; font-size: 10.5px;">
            How to cite: Miguel Hernandez-Gamarra, Alba Salgado-Roo, Eduardo Dominguez, Elena María Goiricelaya Seco, Sara Veiga-Rúa, Lucía F Pedrera-Garbayo, Ángel Carracedo, Catarina Allegue, CARTAR: a comprehensive web tool for identifying potential targets in chimeric antigen receptor therapies using TCGA and GTEx data, Briefings in Bioinformatics, Volume 25, Issue 4, July 2024, bbae326, <a href="https://doi.org/10.1093/bib/bbae326">https://doi.org/10.1093/bib/bbae326</a>.
        </div>
        """, unsafe_allow_html=True)

st.markdown(mystyle, unsafe_allow_html=True)
st.title('Gene expression change across tumors')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.write('This tool can be used to explore fold change expression values ("Primary tumor" compared to "Control" samples) for a gene or gene set of interest in the desired tumors. This will provide a table with the fold cahnge value for each indicated gene in each selected tumor. This will allow to get preliminary information for candidate target genes to check if they are overexpressed in certain tumor and if it can be used to treat more than one tumor.')
st.set_option('deprecation.showPyplotGlobalUse', False)

tumor_options  = ['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS']
scale_options = ['FC','log2(FC)']
# Abbreviation dictionary
abbreviations = {'ACC':'Adrenocortical carcinoma','BLCA':'Bladder Urothelial Carcinoma','BRCA':'Breast invasive carcinoma',
                 'CESC':'Cervical squamous cell carcinoma and endocervical adenocarcinoma','CHOL':'Cholangio carcinoma',
                 'COAD':'Colon adenocarcinoma','DLBC':'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma','ESCA':'Esophageal carcinoma',
                 'GBM':'Glioblastoma multiforme','HNSC':'Head and Neck squamous cell carcinoma','KICH':'Kidney Chromophobe',
                 'KIRC':'Kidney renal clear cell carcinoma','KIRP':'Kidney renal papillary cell carcinoma',
                 'LAML':'Acute Myeloid Leukemia','LGG':'Brain Lower Grade Glioma','LIHC':'Liver hepatocellular carcinoma',
                 'LUAD':'Lung adenocarcinoma','LUSC':'Lung squamous cell carcinoma','OV':'Ovarian serous cystadenocarcinoma',
                 'PAAD':'Pancreatic adenocarcinoma','PCPG':'Pheochromocytoma and Paraganglioma','PRAD':'Prostate adenocarcinoma',
                 'READ':'Rectum adenocarcinoma','SARC':'Sarcoma','SKCM':'Skin Cutaneous Melanoma','STAD':'Stomach adenocarcinoma',
                 'TGCT':'Testicular Germ Cell Tumors','THCA':'Thyroid carcinoma','THYM':'Thymoma',
                 'UCEC':'Uterine Corpus Endometrial Carcinoma','UCS':'Uterine Carcinosarcoma'}

genes = st.text_input('Enter gene symbols of interest (separated by commas or spaces):').upper().strip(' ')
experimental_pm_file = open('Data/HPA_evidence_pm.csv','r')
for line in experimental_pm_file:
    experimental_pm_genes = line.split(',')
correct_genes = []
# Identify if indicated gene is present in the data
exclude = open('Data/no_membrane_genes.csv','r')
for line in exclude:
    no_membrane = line.split(',')
data = pd.read_csv('Data/log2FC_expression.csv')
if ' ' in genes and ',' not in genes:
    genes = genes.replace(' ',',')
elif ', ' in genes or ' ,' in genes:
    genes = genes.replace(' ','')
if  genes == '':
    st.error('Introduce at least one gene symbol. You can try CEACAM6,DPEP1')
else:
    genes = genes.split(',')
    for gene in genes:
        if gene == '':
            pass
        elif gene not in data['gene'].values:
            if gene in no_membrane:
                st.error(f'The protein encoded by {gene} is not located at the membrane')
            else:
                st.error(f'{gene} gene symbol not found')
        else:
            correct_genes.append(gene)
HPA_membrane = ""
for gene in correct_genes:
    if gene in experimental_pm_genes:
        HPA_membrane += f'{gene} , '
if HPA_membrane[-2:] == ', ':
    HPA_membrane = HPA_membrane[:-2]
parts = HPA_membrane.rsplit(',', 1)
HPA_membrane = ' and'.join(parts)
tumors = st.multiselect('Select tumors (optional)', tumor_options)
# Expander to show abbreviation meaning
with st.expander('Extension of tumor abbreviations\' meaning'):
    for abbreviation, meaning in abbreviations.items():
        st.write(f"**{abbreviation}:** {meaning}")
if not tumors :
    tumors = tumor_options
    st.info('If no tumors are selected, all tumors will be displayed')
selection = st.radio('Select scale', scale_options)
if selection == 'FC':
    scale = 'FC'
elif selection == 'log2(FC)':
    scale = 'log2(FC)'
st.info('FC = Fold Change')

if st.button('Show Fold Change'):
    # If there is at least one valid gene
    if correct_genes:
        data = open('Data/log2FC_expression.csv','r')
        t_data = {'Tumor':[]}
        for gene in correct_genes:
            t_data[gene] = []
        firstline = 1
        indexes = [] # List indicating the index corresponding to the introdueced tumors
        for line in data:
            line = line.strip()
            fields = line.split(',')
            # Idnetify the name of the intoduced tumors and their indexes
            if firstline: 
                for field in fields[1:]:
                    if field in tumors:
                        indexes.append(f'{field}_{fields.index(field)}')
                firstline = 0
            # Identify the corresponding fold change for the introduced tumors in the desired scale
            else: 
                for gene in correct_genes:
                    if fields[0] == gene:
                        if scale == 'FC':
                            for index in indexes:
                                tumor = index.split('_')[0]
                                n = int(index.split('_')[1])
                                t_data[gene].append(2**float((fields[n])))
                        else:
                            for index in indexes:
                                tumor = index.split('_')[0]
                                n = int(index.split('_')[1])
                                t_data[gene].append(float(fields[n]))
        for index in indexes:
            tumor = index.split('_')[0]
            t_data['Tumor'].append(tumor)
        # Show table with the results
        table_data = pd.DataFrame(t_data)
        st.header('Data table', divider='rainbow')
        st.write(
            f'The {scale} expression for the specified genes between the "Primary tumor" and "Control" samples is displayed in the table below. To determine whether the expression difference is statistically significant across these conditions, refer to the [**Tumor Gene Expression Tool**](https://cartar-car-targets.streamlit.app/Tumor_gene_expression). For insights into whether any gene is expressed in healthy GTEx tissues, visit the [**Tissue Gene Expression Tool**](https://cartar-car-targets.streamlit.app/Tissue_gene_expression) to assess its specificity. Click on the column names to sort the tumors based on the respective column in ascending or descending order.'
        )
        if 'and' in HPA_membrane:
            st.write(
                f'**{HPA_membrane} have been experimetally reported to be located in the plasma membrane by the Human Protein Atlas.**'
            )
        else:
            st.write(
                f'**{HPA_membrane} has been experimetally reported to be located in the plasma membrane by the Human Protein Atlas.**'
            )   
        st.dataframe(table_data, hide_index=True)
        table = table_data.to_csv(encoding='utf-8', index=False)
        b64 = base64.b64encode(table.encode()).decode()
        href = f'<a href="data:file/csv;base64,{b64}" download="table.csv">Download CSV File</a>'
        st.markdown(href, unsafe_allow_html=True)
    else:
        if genes == '':
            st.error('No gene symbol was introduced')
        else:
            st.error('Please check if the intorduce gene symbol is found or if gene symbols are introduced separated by commas without spaces or separated by one space as in the following example: **EGFR,FGFR1,CD19** or **EGFR FGFR1 CD19**')
create_footer()
