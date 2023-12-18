import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
from math import log2
import pickle
import numpy as np
import seaborn as sns
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
st.title('Gene correlation for logic gating CAR')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.write(
    'This tool facilitates exploring the expression correlation between the two specified genes across all samples within the chosen tumor (data sourced from the TCGA project) and its corresponding controls (data obtained from both the TCGA project and GTEX tissues).This can be helpful for the design of CAR logic gating strategies.')
st.set_option('deprecation.showPyplotGlobalUse', False)

tumor_options  = ['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS']
scale_options = ['TPM','log2(TPM+1)']
gene1 = st.text_input('Select first gene').upper()

# Identify if indicated gene is present in the data
data = pd.read_csv('Data/log2FC_expression.csv')
if gene1 == '':
    st.error('Introduce a gene')
elif gene1 != '' and gene1 not in data['gene'].values:
    st.error(f'{gene1} gene not found')
gene2 = st.text_input('Select second gene').upper()
if gene2 == '':
    st.error('Introduce a gene')
elif gene2 != '' and gene2 not in data['gene'].values:
    st.error(f'{gene2} gene not found')
tumor = st.selectbox('Choose tumor', tumor_options)
selection = st.radio('Choose scale', scale_options)
if selection == 'TPM':
    scale = 'TPM'
elif selection == 'log2(TPM+1)':
    scale = 'log2(TPM+1)'

# Dictionary with the tumor:GTEX healthy tissue correspondance
gtex_tcga = {'ACC':'Adrenal Gland','BLCA':'Bladder','BRCA':'Breast','CESC':'Cervix Uteri', 'COAD':'Colon', 'DLBC':'Blood',
                    'ESCA':'Esophagus','GBM':'Brain','KICH':'Kidney','KIRC':'Kidney','KIRP':'Kidney','LAML':'Bone Marrow','LGG':'Brain',
                    'LIHC':'Liver','LUAD':'Lung','LUSC':'Lung','OV':'Ovary','PAAD':'Pancreas','PRAD':'Prostate','READ':'Colon','SKCM':'Skin',
                    'STAD':'Stomach','TGCT':'Testis','THCA':'Thyroid','THYM':'Blood','UCEC':'Uterus','UCS':'Uterus'}

if st.button(f'Show correlation'):
    if gene1 != '' and gene2 != '':
        # If gene in data
        if gene1 in data['gene'].values and gene2 in data['gene'].values:  
            if 'A' <= gene1[0] <= 'C':
                gtex1 = 'Data/gtex_AC.pkl'
                tcga1 = 'Data/tcga_AC.pkl'
            elif 'D' <= gene1[0] <= 'J':
                gtex1 = 'Data/gtex_DJ.pkl'
                tcga1 = 'Data/tcga_DJ.pkl'
            elif 'K' <= gene1[0] <= 'N':
                gtex1 = 'Data/gtex_KN.pkl'
                tcga1 = 'Data/tcga_KN.pkl'
            elif 'O' <= gene1[0] <= 'R':   
                gtex1 = 'Data/gtex_OR.pkl'
                tcga1 = 'Data/tcga_OR.pkl'
            elif 'S' <= gene1[0] <= 'T':   
                gtex1 = 'Data/gtex_ST.pkl'
                tcga1 = 'Data/tcga_ST.pkl'
            elif 'U' <= gene1[0] <= 'Z':   
                gtex1 = 'Data/gtex_UZ.pkl'
                tcga1 = 'Data/tcga_UZ.pkl'
            if 'A' <= gene2[0] <= 'C':
                gtex2 = 'Data/gtex_AC.pkl'
                tcga2 = 'Data/tcga_AC.pkl'
            elif 'D' <= gene2[0] <= 'J':
                gtex2 = 'Data/gtex_DJ.pkl'
                tcga2 = 'Data/tcga_DJ.pkl'
            elif 'K' <= gene2[0] <= 'N':
                gtex2 = 'Data/gtex_KN.pkl'
                tcga2 = 'Data/tcga_KN.pkl'
            elif 'O' <= gene2[0] <= 'R':   
                gtex2 = 'Data/gtex_OR.pkl'
                tcga2 = 'Data/tcga_OR.pkl'
            elif 'S' <= gene2[0] <= 'T':   
                gtex2 = 'Data/gtex_ST.pkl'
                tcga2 = 'Data/tcga_ST.pkl'
            elif 'U' <= gene2[0] <= 'Z':   
                gtex2 = 'Data/gtex_UZ.pkl'
                tcga2 = 'Data/tcga_UZ.pkl'
            # Create the dicitionary with all the data
            groups = [] # Groups of tumor (Primary or Normal)
            values1 = [] # Gene1 expression values 
            values2 = [] # Gene2 expression values 
            #If both genes in same file get the desired data
            if gtex1 == gtex2:
                with open(tcga1, 'rb') as archivo:
                    tcga = pickle.load(archivo)
                for group in tcga[gene1][tumor].keys():
                    for value in tcga[gene1][tumor][group]:
                        groups.append(group)
                        if scale == 'log2(TPM+1)':
                            value = log2(value+1)
                        values1.append(value)
                for group in tcga[gene2][tumor].keys():
                    for value in tcga[gene2][tumor][group]:
                        if scale == 'log2(TPM+1)':
                            value = log2(value+1)
                        values2.append(value)
                if tumor in gtex_tcga.keys():            
                    with open(gtex1, 'rb') as archivo:
                        gtex = pickle.load(archivo)
                    for value in gtex[gene1][gtex_tcga[tumor]]:
                        groups.append('Normal')
                        if scale == 'log2(TPM+1)':
                            value = log2(value+1)
                        values1.append(value)   
                    for value in gtex[gene2][gtex_tcga[tumor]]:
                        if scale == 'log2(TPM+1)':
                            value = log2(value+1)
                        values2.append(value)   
            #If genes in different files get data from respective files
            else:
                with open(tcga1, 'rb') as archivo:
                    tcga1 = pickle.load(archivo)
                for group in tcga1[gene1][tumor].keys():
                    for value in tcga1[gene1][tumor][group]:
                        groups.append(group)
                        if scale == 'log2(TPM+1)':
                            value = log2(value+1)
                        values1.append(value)
                with open(tcga2, 'rb') as archivo:
                    tcga2 = pickle.load(archivo)
                for group in tcga2[gene2][tumor].keys():
                    for value in tcga2[gene2][tumor][group]:
                        if scale == 'log2(TPM+1)':
                            value = log2(value+1)
                        values2.append(value)
                if tumor in gtex_tcga.keys():            
                    with open(gtex1, 'rb') as archivo:
                        gtex1 = pickle.load(archivo)
                    for value in gtex1[gene1][gtex_tcga[tumor]]:
                        groups.append('Normal')
                        if scale == 'log2(TPM+1)':
                            value = log2(value+1)
                        values1.append(value)  
                    with open(gtex2, 'rb') as archivo:
                        gtex2 = pickle.load(archivo)    
                    for value in gtex2[gene2][gtex_tcga[tumor]]:
                        if scale == 'log2(TPM+1)':
                            value = log2(value+1)
                        values2.append(value)  
            print(len(groups))
            print(len(values1))
            print(len(values2))
            data = {'Sample':groups, f'{gene1} expression':values1, f'{gene2} expression':values2}
            # Create correlation plot
            df = pd.DataFrame(data)          
            fig = px.scatter(data, x=f'{gene2} expression', y=f'{gene1} expression', color='Sample', color_discrete_sequence=['#20b2aa','#d2b48c'])
            # Customize graph titile and axis names
            if scale == 'log2(TPM+1)':
                fig.update_layout(title=f'{gene1} and {gene2} expression correlation in {tumor} samples', title_x=0.2, xaxis_title=f'{gene2} expression in log2(TPM+1)', yaxis_title= f'{gene1} expression in log2(TPM+1)')
            else:
                fig.update_layout(title=f'{gene1} and {gene2} expression correlation in {tumor} samples', title_x=0.2, xaxis_title=f'{gene2} expression in TPM', yaxis_title= f'{gene1} expression in TPM')
            # Create table with the results
            table_data = pd.DataFrame(data)
            if scale == 'TPM':
                plt.ylabel(f'{gene1} expression in TPM')
                plt.xlabel(f'{gene2} expression in TPM')
            else: 
                plt.ylabel(f'{gene1} expression in log2(TPM+1)')
                plt.xlabel(f'{gene2} expression in log2(TPM+1)')
            plt.xticks(rotation=45, ha='right')
            # Show the graph and table
            st.header('Correlation plot', divider='rainbow')
            st.plotly_chart(fig,use_container_width=True)
            st.write(
                f'The presented figure illustrates the expression correlation between {gene1} and {gene2} in {scale} within the {tumor}, showcasing the comparison of expression levels between "Primary Tumor" and "Normal" samples.'
            )
            st.header('Data table', divider='rainbow')
            st.write(
                f'The table below displays the expression values of each gene for all samples. You can enhance your exploration by clicking on the column names to sort the tumors based on that column either from highest to lowest or vice versa. Please note that **p-values under 0.001 are rounded to 0**; for the complete decimal value, click on the respective cell.''
            )
            st.dataframe(table_data, hide_index=True)
            st.write(
                'This table can be downloaded in CSV format for opening in Excel or other software compatible with this file format.'
            )   