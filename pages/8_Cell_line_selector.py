import streamlit as st
import pandas as pd
from math import log2
import requests
import csv
import plotly.express as px

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
st.title('Cancer cell line selector')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.write(
    'This tool facilitates the identification of cancer cell lines from the Cancer Cell Line Encyclopedia (CCLE) expressing the indicated gene above or below the indicated threshold. You can specify the cancer cell line lineage of interest, but if none is specified it will provide information anout all lineages. The expression data for cancer cell lines is sourced from [DepMap Public 23Q2 Primary Files](https://depmap.org/portal/download/all/). The tool generates an interactive barplot illustrating the gene\'s expression in those cell lines by lineage, meeting the established threshold, accompanied by a table containing pertinent data for all relevant cell lines.'
)
st.set_option('deprecation.showPyplotGlobalUse', False)

tissue_options = ['Ovary/Fallopian Tube', 'Myeloid', 'Bowel', 'Skin', 'Bladder/Urinary Tract', 'Lung', 'Kidney', 'Breast', 'Lymphoid', 'Pancreas', 'CNS/Brain', 'Soft Tissue', 'Bone', 'Fibroblast', 'Esophagus/Stomach', 'Thyroid', 'Peripheral Nervous System', 'Pleura', 'Prostate', 'Biliary Tract', 'Head and Neck', 'Uterus', 'Ampulla of Vater', 'Liver', 'Cervix', 'Eye', 'Vulva/Vagina', 'Adrenal Gland', 'Testis', 'Other', 'Normal']
expression_options = ['Overexpression', 'Underexpression']
scale_options = ['TPM','log2(TPM+1)']
data = open('Data/genes_cells.csv','r')
#data = open(r'C:\Users\34656\Pipeline for TCGA_GTEX\Data\Processed\genes_cells.csv', 'r') # Open cell line expression data
data_list = []
for line in data:
    data_list.extend(line.split(','))
gene = st.text_input('Select gene').upper()  # Introduce gene
if gene == '':
    st.error('Introduce a gene')
if 'MORF' not in gene:
    if 'ORF' in gene:
        gene = gene.replace('ORF','orf')
if gene != '' and gene not in data_list:
    st.error(f'{gene} gene not found')
tumors = st.multiselect('Select cancer cell line lineage (optional)', tissue_options) # Introduce tumor tissues or none if you want all tumors
if not tumors:
    #tumors = tissue_options
    st.success('If no cancer cell lineage is selected, all lineages will be displayed')
expression =  st.radio('Select underexpression or overexpression', expression_options)# Introduce under to look for cell with expression values below the threshold or over for expression values aboce the threshold
scale =  st.radio('Select scale', scale_options) # Introduce scale 
value =  st.text_input('Select threshold expression value:') # Introduce threshold value
if value:
    if ',' in value:
        value = value.replace(',','.')
    value = float(value)
    if scale == 'TPM':
        threshold = log2(value+1)
    else:
        threshold = value

if st.button('Find cell lines'):
    url_data = 'https://gitlab.com/gmx2/CARTAR/-/raw/main/Expression_Public_23Q2_filtered.csv'
    data = pd.read_csv(url_data, index_col=0)
    #data = pd.read_csv(r'C:\Users\34656\Pipeline for TCGA_GTEX\Data\Processed\Expression_Public_23Q2_filtered.csv', index_col=0) # Open cell line expression data
    url_metadata = 'https://gitlab.com/gmx2/CARTAR/-/raw/main/cell_line_metadata.csv'
    metadata = pd.read_csv(url_metadata)
    #metadata = pd.read_csv(r'C:\Users\34656\Pipeline for TCGA_GTEX\Data\Processed\cell_line_metadata.csv')
    # Identify cell lines belonging to specified tumor
    cell_lines = []
    if tumors:
        indexes = metadata[metadata['OncotreeLineage'].isin(tumors)].index.tolist() # Create a list with the row number corresponding to specified tumors
        cell_lines.extend(metadata.loc[indexes, 'ModelID'].tolist())
    else:
        cell_lines.extend(metadata['ModelID'].tolist())
    # If inidcated gene in dataset
    if gene in data.columns:
        # Identify cell lines with expression related to threshold
        data = pd.DataFrame({gene: data[gene]})
        data2 = data.T
        desired_lines = [] # List to store the cell lines with desired expression
        values = [] # List to store the cell line expression
        for column in data2.columns:
            if column in cell_lines:
                exp = float(data2[column].tolist()[0])
                if expression == 'Underexpression':
                    if exp <= threshold:
                        desired_lines.append(column)
                        values.append(exp)
                else:
                    if exp >= threshold:
                        desired_lines.append(column)   
                        values.append(exp)       
    else:
        st.error(f'{gene} gene data not available')
    # Select column corresponding to the indicated gene and rows for desired cell lines
    data3 = data.loc[cell_lines, [gene]]
    # Identify maximum and minimum expression values for the gene
    max_value = data3[gene].max()
    min_value = data3[gene].min()
    if scale == 'TPM':
        max_value = 2**(max_value)-1
        min_value = 2**(min_value)-1
    # Obtain all relevant information about the selected cell lines
    table = {'Cell line':[],f'{gene} expression':[],'Catalog Number':[],'Lineage':[],'Primary Disease':[],'Subtype':[],'Code':[]}
    desired_rows = metadata[metadata['ModelID'].isin(desired_lines)].index.tolist()
    table['Cell line'].extend(metadata.loc[desired_rows, 'CellLineName'].tolist())
    if scale == 'TPM':
        log_values = []
        for value in values:
            log_values.append(2**value-1)
        table[f'{gene} expression'] = log_values
    else:
        table[f'{gene} expression'] = values
    table['Catalog Number'].extend(metadata.loc[desired_rows, 'CatalogNumber'].tolist())
    table['Lineage'].extend(metadata.loc[desired_rows, 'OncotreeLineage'].tolist())
    table['Primary Disease'].extend(metadata.loc[desired_rows, 'OncotreePrimaryDisease'].tolist())
    table['Subtype'].extend(metadata.loc[desired_rows, 'OncotreeSubtype'].tolist())
    table['Code'].extend(metadata.loc[desired_rows, 'OncotreeCode'].tolist())
    table_data = pd.DataFrame(table)
    if table_data.size > 0:
        # Generate interactive barplot
        fig_data = {'Cell line': table['Cell line'], 'Value': table[f'{gene} expression'], 'Lineage':table['Lineage'], 'Label':[f'{cell_line}: {round(expression, 2)} {scale}' for cell_line, expression in zip(table['Cell line'], table[f'{gene} expression'])]}
        df = pd.DataFrame(fig_data)
        df = df.sort_values(by='Value', ascending=False)
        fig = px.bar(df, x='Cell line', y='Value', color='Lineage', color_discrete_sequence=['#CC6262', '#72CC62', '#6A76FC', '#FED4C4', '#FE00CE', '#0DF9FF', '#F6F926', '#FF9616', '#479B55', '#EEA6FB', '#DC587D', '#D626FF', '#6E899C', '#00B5F7', '#B68E00', '#C9FBE5', '#FF0092', '#22FFA7', '#E3EE9E', '#86CE00', '#BC7196', '#7E7DCD', '#FC6955', '#E48F72'], custom_data=['Label'])
        fig.update_traces(hovertemplate='%{customdata}')
        # Customize graph titile and axis names
        if scale == 'log2(TPM+1)':
            fig.update_layout(title=f'{gene} expression comparison between desired cell lines', title_x=0.25, xaxis_title='Cell lines', yaxis_title= f'{gene} expression in log2(TPM+1)')
        else:
            fig.update_layout(title=f'{gene} expression comparison between desired cell lines', title_x=0.25, xaxis_title='Cell lines', yaxis_title= f'{gene} expression in TPM')
        # Show interactive plot
        st.header('Interactive barplot', divider='rainbow')
        st.plotly_chart(fig,use_container_width=True)
        st.write(
            f'The above figure shows an interactive barplot for {gene} expression in {scale} across all cell lines belonging to the indicated lineage and meeting the stablished threshold. This can be useful to identify candidate target and negative control cell lines.'
        )
        # Show table
        st.header('Data table', divider='rainbow')
        st.write(
            f'All relevant data for the selected cell lines is presented in the table below, encompassing lineage, primary diseases, and disease subtypes. You can click on the column names to arrange the genes based on that column either in ascending or descending order. Further details for each cell line can be viewed by clicking on the respective cell, revealing the values with all the decimals.'
        )
        st.dataframe(table_data, hide_index=True)
        st.write(
            'This table can be downloaded in CSV format.'
        ) 
    else:
        if expression == 'Underexpression':
            st.error(f'No cell lines with {gene} expression under {value} {scale} have been found. Minimum value expression for {gene} in selected cell lines is {min_value} {scale}')
        else:
            st.error(f'No cell lines with {gene} expression over {value} {scale} have been found. Maximum value expression for {gene} in selected cell lines is {max_value} {scale}')