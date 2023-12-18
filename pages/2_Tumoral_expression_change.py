import streamlit as st
import pandas as pd
from math import log2
import requests

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
st.title('Gene expression change across tumors')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.write(
    'This tool will generate a table indicating the fold change (FC) expression between \"Primary tumor\" and \"Control\" samples for the specified genes, considering the selected scale and tumor types. If no tumor is specified, it will display the expression changes across all tumors.'
)
st.set_option('deprecation.showPyplotGlobalUse', False)

tumor_options  = ['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS']
scale_options = ['FC','log2(FC)']
genes = st.text_input('Enter genes of interest (separated by commas without spaces):').upper()
correct_genes = []
# Identify if indicated gene is present in the data
data = pd.read_csv('Data/log2FC_expression.csv')
if genes == '':
    st.error('Introduce a gene')
else:
    if ' ' in genes:
        st.error('Remove spaces')
    else:
        genes = genes.split(',')
        for gene in genes:
            if gene not in data['gene'].values:
                st.error(f'{gene} gene not found')
            else:
                correct_genes.append(gene)
tumors = st.multiselect('Select tumors', tumor_options)
if not tumors :
    tumors = tumor_options
    st.info('If no tumors are selected, all tumors will be displayed')
selection = st.radio('Select scale', scale_options)
if selection == 'FC':
    scale = 'FC'
elif selection == 'log2(FC)':
    scale = 'log2(FC)'

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
            # Identify the correspondinf fold change for the introduced tumors in the desired scale
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
            f'The {scale} expression for the specified genes between the "Primary tumor" and "Control" samples is displayed in the table below. To determine whether the expression difference is statistically significant across these conditions, refer to the [**Tumor Gene Expression Tool**](http://localhost:8501/Tumor_gene_expression). For insights into whether any gene is expressed in healthy GTEX tissues, visit the [**Tissue Gene Expression Tool**](http://localhost:8501/Tissue_gene_expression) to assess its specificity. Click on the column names to sort the tumors based on the respective column in ascending or descending order. Please note that **p-values under 0.001 are rounded to 0**; for the complete decimal value, click on the respective cell.'
        )
        st.dataframe(table_data, hide_index=True)
        st.write(
            'This table can be downloaded in CSV format for opening in Excel or other software compatible with this file format.'
        ) 