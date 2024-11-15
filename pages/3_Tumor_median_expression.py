import streamlit as st
import numpy as np
import math
from math import log2
import os
import pandas as pd
import plotly.express as px
from io import StringIO
import requests
import csv
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
        st.markdown("<br>" * 1, unsafe_allow_html=True)
        st.markdown("""
                <style>
                .footer-content {
                    background-color: #f0f2f6;
                    color: #262730;
                    padding: 10px;
                    text-align: center;
                    font-size: 10.5px;
                }
                @media (prefers-color-scheme: dark) {
                    .footer-content {
                        background-color: #262730;
                        color: #ffffff;
                    }
                    .footer-content a {
                        color: #4da6ff;
                    }
                }
                </style>
                <div class="footer-content">
                    How to cite: Miguel Hernandez-Gamarra, Alba Salgado-Roo, Eduardo Dominguez, Elena María Goiricelaya Seco, Sara Veiga-Rúa, Lucía F Pedrera-Garbayo, Ángel Carracedo, Catarina Allegue, CARTAR: a comprehensive web tool for identifying potential targets in chimeric antigen receptor therapies using TCGA and GTEx data, Briefings in Bioinformatics, Volume 25, Issue 4, July 2024, bbae326, <a href="https://doi.org/10.1093/bib/bbae326">https://doi.org/10.1093/bib/bbae326</a>.
                </div>
                """, unsafe_allow_html=True)

st.markdown(mystyle, unsafe_allow_html=True)
st.title('Gene median expression across tumors')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.write('This tool can be used to easily visualize a barplot showcasing the median expression value of a selected gene across primary tumor, metastatic (if available), and control samples across tumors of interest. Besides, median expression values and sample size of each group are reported in table format. This will allow to get information about a candidate target gene identified with the **Tumor-associated antigens tool** to check its expression in all sample groups and whether it can be used to treat more than one tumor.')
st.set_option('deprecation.showPyplotGlobalUse', False)

tumor_options  = ['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS']
scale_options = ['TPM','log2(TPM+1)']
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

gene = st.text_input('Enter gene symbol').upper().strip(' ')
experimental_pm_file = open('Data/HPA_evidence_pm.csv','r')
for line in experimental_pm_file:
    experimental_pm_genes = line.split(',')
# Identify if indicated gene is present in the data
data = pd.read_csv('Data/log2FC_expression.csv')
exclude = open('Data/no_membrane_genes.csv','r')
for line in exclude:
    no_membrane = line.split(',')
if gene == '':
    st.error('Introduce gene symbol. You can try CEACAM6')
elif gene != '' and gene not in data['gene'].values:
    if gene in no_membrane:
        st.error(f' The protein encoded by {gene} is not located at the membrane')
    else:
        st.error(f'{gene} gene symbol not found')
tumors = st.multiselect('Select tumors (optional)', tumor_options)
# Expander to show abbreviation meaning
with st.expander('Extension of tumor abbreviations\' meaning'):
    for abbreviation, meaning in abbreviations.items():
        st.write(f"**{abbreviation}:** {meaning}")
if not tumors :
    tumors = tumor_options
    st.info('If no tumors are selected, all tumors will be displayed')
selection = st.radio('Select scale', scale_options)
if selection == 'TPM':
    scale = 'TPM'
elif selection == 'log2(TPM+1)':
    scale = 'log2(TPM+1)'
st.info('TPM = Transcript Per Million')

if st.button('Create barplot'):
    if gene != '' and gene in data['gene'].values: 
        # Get the required data
        t_data = {'Tumor':[],'Metastatic median':[], 'Primary tumor median':[], 'Control median':[]}         
        categories = [] # List conatining the name of introduced valid tumors
        for tumor in tumors:
            categories.append(tumor)
            t_data['Tumor'].append(tumor)
        categories = sorted(categories) # Sort tumor named alphabetically
        metastatic = [] # List with median values of metastatic sample groups (NA if not available)
        primary = [] # List with median values of primary tumor sample groups
        normal = [] # List with median values of control sample groups
        data = open('Data/targetable_gene_Tpm_TumorVsControl_final.csv','r') # Read expression data for the samples
        m_indexes = [] # List with indexes of metastatic sample groups
        p_indexes = [] # List with indexes values of primary tumor sample groups
        n_indexes = [] # List with indexes values of control sample groups
        firstline = 1
        for line in data:
            fields = line.split(',')
            # Identify position of fields belonging to each sample group of the selected tumors in the firstline of the file
            if firstline: 
                n = 0
                for field in fields[1:]:
                    n += 1
                    if field in tumors:
                        control = fields[n+1]
                        p_indexes.append(fields.index(field))
                        if 'Metastatic' in control:
                            m_indexes.append(n+1)
                            n_indexes.append(n+2)
                        else:
                            m_indexes.append(-1)
                            n_indexes.append(n+1)
                firstline = 0
            else: 
                # Identify median value for each sample group in the selected tumors and selected scale, and add it to corresponding list
                if gene == fields[0]:
                    if scale == 'log2(TPM+1)':
                        for index in m_indexes:
                            if index == -1:
                                metastatic.append(np.nan) # Add NA if metastatic data not available
                                t_data['Metastatic median'].append(np.nan)
                            else:
                                metastatic.append(log2(float(fields[index])+1))
                                t_data['Metastatic median'].append(log2(float(fields[index])+1))
                        for index in p_indexes:
                            primary.append(log2(float(fields[index])+1))
                            t_data['Primary tumor median'].append(log2(float(fields[index])+1))
                        for index in n_indexes: 
                            normal.append(log2(float(fields[index])+1))
                            t_data['Control median'].append(log2(float(fields[index])+1))
                    else:
                        for index in m_indexes:
                            if index == -1:
                                metastatic.append(np.nan)
                                t_data['Metastatic median'].append(np.nan)
                            else:
                                metastatic.append(float(fields[index]))
                                t_data['Metastatic median'].append(float(fields[index]))
                        for index in p_indexes:
                            primary.append(float(fields[index]))
                            t_data['Primary tumor median'].append(float(fields[index]))
                        for index in n_indexes: 
                            normal.append(float(fields[index]))    
                            t_data['Control median'].append(float(fields[index]))
        classes = ['Metastatic', 'Primary', 'Control']
        values = [None]*len(categories)
        for i in range(len(categories)):
            for k in range(len(classes)):
                if k == 0:
                    values[i] = [metastatic[i]]
                elif k == 1:
                    values[i].append(primary[i])
                elif k == 2:
                    values[i].append(normal[i])
        # Create table with the results
        table_data = pd.DataFrame(t_data)
        # Creation of the dictionary containing all needed information for the pltypeot
        data = {'Tumor type': [], 'Value': [], 'Sample': [], 'Label': []}
        for j, category in enumerate(categories):
            for i, condition in enumerate(classes):
                data['Tumor type'].append(category)
                data['Sample'].append(condition)
                data['Value'].append(values[j][i])
                data['Label'].append(f'{category}: {round(values[j][i],2)}')
        # Creation of the interactive bar plot  
        df = pd.DataFrame(data)
        colors = ['#8C8C8C', '#1FA698', '#D9B991']
        fig = px.bar(df, x='Tumor type', y='Value', color='Sample', color_discrete_sequence=colors,custom_data=['Label'])
        fig.update_traces(hovertemplate='%{customdata}')
        # Customize graph titile and axis names
        if scale == 'log2(TPM+1)':
            fig.update_layout(title=f'{gene} expression comparison between tumoral conditions', title_x=0.165, xaxis_title='Tumor type', yaxis_title= f'{gene} expression in log2(TPM+1)', barmode='group')
        else:
            fig.update_layout(title=f'{gene} expression comparison between tumoral conditions', title_x=0.165, xaxis_title='Tumor type', yaxis_title= f'{gene} expression in TPM', barmode='group')
        # Show the graph and table
        st.header('Interactive barplot', divider='rainbow')
        st.plotly_chart(fig,use_container_width=True)
        st.write(
            f'The interactive barplot above illustrates the expression of {gene} in {scale} across "Metastatic," "Primary tumor," and "Control" samples from selected tumors. The visualization allows for a comparison of expression between these sample groups. Note that no statistical significance is indicated here. For information on the statistical significance between "Primary tumor" and "Control" samples, please refer to [**Tumor gene expression tool**](https://cartar-car-targets.streamlit.app/Tumor_gene_expression).'
        )
        st.header('Data table', divider='rainbow')
        st.write(
            f'The median values for each sample group in the selected tumors are displayed in the table below. The median values for each sample group in the selected tumors are displayed in the table below. If you want to explore the log2(Fold Change) between \"Primary tumor\" and \"Control\" samples, along with the statistical significance, please visit [**Tumor Gene Expression Tool**](https://cartar-car-targets.streamlit.app/Tumor_gene_expression). You can click on the column names to sort the tumors based on that column in ascending or descending order. Please note that **p-values under 0.001 are rounded to 0**; for the complete decimal value, click on the respective cell.'
        )
        if gene in experimental_pm_genes:
            st.write(
                f'**{gene} has been experimetally reported to be located in the plasma membrane by the Human Protein Atlas.**'
            )
        st.dataframe(table_data, hide_index=True)
        table = table_data.to_csv(encoding='utf-8', index=False)
        b64 = base64.b64encode(table.encode()).decode()
        href = f'<a href="data:file/csv;base64,{b64}" download="table.csv">Download CSV File</a>'
        st.markdown(href, unsafe_allow_html=True) 
    elif gene == '':
        st.error('No gene symbol was introduced')
    elif gene not in data['gene'].values:
        st.error(f'{gene} gene symbol not found')  
create_footer()
