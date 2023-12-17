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

st.set_page_config(page_title='CARTAR', page_icon='logo.png',layout='wide')
mystyle = '''
    <style>
        p {
            text-align: justify;
        }
    </style>
    '''
@st.cache_data
def get_base64_of_bin_file(png_file):
    with open(png_file, "rb") as f:
        data = f.read()
    return base64.b64encode(data).decode()

def build_markup_for_logo(
    png_file,
    background_position="50% 10%",
    image_width="90%",
    image_height="",
):
    binary_string = get_base64_of_bin_file(png_file)
    return """
            <style>
                [data-testid="stSidebarNav"] {
                    background-image: url("data:image/png;base64,%s");
                    background-repeat: no-repeat;
                    background-position: %s;
                    background-size: %s %s;
                }
                [data-testid="stSidebarNav"]::before {
                content: "";
                margin-left: 20px;
                margin-top: 20px;
                font-size: 15px;
                position: relative;
                top: 100px;
                }
            </style>
            """ % (
        binary_string,
        background_position,
        image_width,
        image_height,
    )

def add_logo(png_file):
    logo_markup = build_markup_for_logo(png_file)
    st.markdown(
        logo_markup,
        unsafe_allow_html=True,
    )

add_logo('logo_v2.png')

st.markdown(mystyle, unsafe_allow_html=True)
st.title('Gene median expression across tumors')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.write(
    'This tool will generate an interactive barplot comparing median expression values of the chosen gene in the selected scale for the \"Metastatic\", \"Primary tumor\" and \"Control\" samples of the indicated tumors and a table with all relevant data. Data is obtained from TCGA project and GTEX expression data has also been included in the correponding control sample group. If no tumor is indicated it will report the expression data of the indicated gene in all tumors.'
)
st.set_option('deprecation.showPyplotGlobalUse', False)

tumor_options  = ['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS']
scale_options = ['TPM','log2(TPM+1)']

gene = st.text_input('Select gene').upper()
# Identify if indicated gene is present in the data
data = pd.read_csv('Data/log2FC_expression.csv')
if gene == '':
    st.error('Indicate a gene')
elif gene != '' and gene not in data['gene'].values:
    st.error('Indicated gene not found')
tumors = st.multiselect('Choose tumors', tumor_options)
if not tumors :
    tumors = tumor_options
    st.success('If not tumors are seleceted all tumors will be shown')
selection = st.radio('Choose scale', scale_options)
if selection == 'TPM':
    scale = 'TPM'
elif selection == 'log2(TPM+1)':
    scale = 'log2(TPM+1)'

if st.button('Create barplot'):
    if gene != '': 
        # Get the required data
        t_data = {'Tumor':[],'Metastatic median':[], 'Primary tumor median':[], 'Control median':[]}         
        categories = [] # List conatining the name of introduced valid tumors
        for tumor in tumors:
            categories.append(tumor)
            t_data['Tumor'].append(tumor)
        categories = sorted(categories) # Sort tumor named alphabetically
        metastatic = [] # List with median values of metastatic sample groups (NA if not available)
        primary = [] # List with median values of primary tumor sample groups
        normal = [] # List with median values of normal sample groups
        data = open('Data/targetable_gene_Tpm_TumorVsControl_final.csv','r') # Read expression data for the samples
        m_indexes = [] # List with indexes of metastatic sample groups
        p_indexes = [] # List with indexes values of primary tumor sample groups
        n_indexes = [] # List with indexes values of normal sample groups
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
        classes = ['Metastatic', 'Primary', 'Normal']
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
            f'The above figure shows an interactive barplot for {gene} expression in {scale} across \"Metastatic\", \"Primary tumor\" and \"Control\" samples of selected tumors to compare the expression between these sample groups. No statistical significance is shown. If you are interested in knowing the statistical significance beteween \"Primary tumor\" and \"Control\" samples you must go to [**Tumor gene expression tool**](http://localhost:8501/Tumor_gene_expression).'
        )
        st.header('Data table', divider='rainbow')
        st.write(
            f'Median of each sample group for the selected tumors are shown in the table below. If you are interested in knowing the log2(Fold Change) between \"Primary tumor\" and \"Control\" samples, as well as the statistical significance, you must go to [**Tumor gene expression tool**](http://localhost:8501/Tumor_gene_expression). You can click in the column names to order the tumors according to that column from higher to lower or viceversa. By clicking in a cell you can see the value with all the decimals.'
        )
        st.dataframe(table_data, hide_index=True)
        st.write(
            'This table can be download in CSV format to be open in Excel or oher softwares compatible with this format.'
        )  
