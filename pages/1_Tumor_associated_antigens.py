import streamlit as st
import pandas as pd
from math import log2
import math
import numpy as np
import plotly.express as px
from scipy import stats
from scipy.stats import mannwhitneyu
import pickle

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


def get_tcga_info(list):
    for gene in list:
        for group in tcga[gene][tumor].keys():
            for value in tcga[gene][tumor][group]:
                genes.append(gene)
                groups.append(group)
                values.append(value)

def get_gtex_info(list):
    for gene in list:
        tissue = gtex_tcga[tumor]
        group = 'Normal'
        for value in gtex[gene][tissue]:
            genes.append(gene)
            groups.append(group)
            values.append(value)

st.markdown(mystyle, unsafe_allow_html=True)
st.title('Tumor-associated antigens identification')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.write(
    'This tool allows the identification of tumor-associated antigens of the indicated tumor. You can introduce the desired Fold Change (FC) between \"Primary tumor\" and \"Control\" samples to be used as a threshold. \"Primary tumor\" data is obtained form TCGA, while \"Control\" data is obtained form the combinaton of TCGA and GTEX. This tool will return a table with all genes with a FC higher than the one specified and the associated p-values.'
)
st.set_option('deprecation.showPyplotGlobalUse', False)

tumor_options = ['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS']
tumor = st.selectbox('Choose tumor', tumor_options)
if not tumor:
    st.error('Select tumor of interest')
FC = st.text_input('Introduce Fold Change value used as threshold to identify tumor-associated genes:')
if FC:
    if ',' in FC:
        FC = FC.replace(',','.')
    FC = float(FC)
    if FC < 1.8:
        st.info('Please take into account that this could take some time')
    if FC < 1:
        st.error('Fold changes under 1 are not recommended')


# Dictionary with the tumor:GTEX healthy tissue correspondance
gtex_tcga = {'ACC':'Adrenal Gland','BLCA':'Bladder','BRCA':'Breast','CESC':'Cervix Uteri', 'COAD':'Colon', 'DLBC':'Blood',
                'ESCA':'Esophagus','GBM':'Brain','KICH':'Kidney','KIRC':'Kidney','KIRP':'Kidney','LAML':'Bone Marrow','LGG':'Brain',
                'LIHC':'Liver','LUAD':'Lung','LUSC':'Lung','OV':'Ovary','PAAD':'Pancreas','PRAD':'Prostate','READ':'Colon','SKCM':'Skin',
                'STAD':'Stomach','TGCT':'Testis','THCA':'Thyroid','THYM':'Blood','UCEC':'Uterus','UCS':'Uterus'}

if st.button('Show tumor-associated genes'):
    log2FC = log2(FC)
    # Show genes overexpressed by the indicated FC
    data = pd.read_csv('Data/log2FC_expression.csv')
    final = {'Gene':[],'log2(FC)':[],'FC':[],f'{tumor} median':[], f'{tumor} sample size':[], 'Control median':[],'Control sample size':[],'Significance':[]} 
    figure = {'Gene':[],'log2(FC)':[],'pvalue':[],'Legend':[]}
    for column in data.columns:
        if column == tumor:
            condition = data[column] >= log2FC
            indexes = data.loc[condition].index
            max_FC = 2**data[column].max()
            if indexes.size > 0:
                for index in indexes: 
                    gene = data['gene'][index]
                    log2_FC = data[column][index]
                    FC = 2 ** log2_FC
                    final['Gene'].append(gene)
                    figure['Gene'].append(gene)
                    final['log2(FC)'].append(log2_FC)
                    figure['log2(FC)'].append(log2_FC)
                    final['FC'].append(FC) 
                # Get requested information
                genes = [] # List with genes
                groups = [] # Gruops of tumor (Primary or Normal)
                values = [] # Expression values
                # Determine p-vlaue of differential expression
                AC_genes = []
                DJ_genes = []
                KN_genes = []
                OR_genes = []
                ST_genes = []
                UZ_genes = []
                for gene in final['Gene']:
                    if 'A' <= gene[0] <= 'C':
                        AC_genes.append(gene)
                        data = 'Data/tcga_AC.pkl'
                    elif 'D' <= gene[0] <= 'J':   
                        DJ_genes.append(gene)
                        data = 'Data/tcga_DJ.pkl'
                    elif 'K' <= gene[0] <= 'N':
                        KN_genes.append(gene)
                        data = 'Data/tcga_KN.pkl'
                    elif 'O' <= gene[0] <= 'R':
                        OR_genes.append(gene)
                        data = 'Data/tcga_OR.pkl'
                    elif 'S' <= gene[0] <= 'T':
                        ST_genes.append(gene)
                        data = 'Data/tcga_ST.pkl'
                    elif 'U' <= gene[0] <= 'Z':
                        UZ_genes.append(gene)
                        data = 'Data/tcga_UZ.pkl'
                if AC_genes != []:
                    with open('Data/tcga_AC.pkl', 'rb') as archivo:
                        tcga = pickle.load(archivo)
                    get_tcga_info(AC_genes)
                if DJ_genes != []:
                    with open('Data/tcga_DJ.pkl', 'rb') as archivo:
                        tcga = pickle.load(archivo)
                    get_tcga_info(DJ_genes)
                if KN_genes != []:
                    with open('Data/tcga_KN.pkl', 'rb') as archivo:
                        tcga = pickle.load(archivo)
                    get_tcga_info(KN_genes)
                if OR_genes != []:
                    with open('Data/tcga_OR.pkl', 'rb') as archivo:
                        tcga = pickle.load(archivo)
                    get_tcga_info(OR_genes)
                if ST_genes != []:
                    with open('Data/tcga_ST.pkl', 'rb') as archivo:
                        tcga = pickle.load(archivo)
                    get_tcga_info(ST_genes)
                if UZ_genes != []:
                    with open('Data/tcga_UZ.pkl', 'rb') as archivo:
                        tcga = pickle.load(archivo)
                    get_tcga_info(UZ_genes)
     
                if tumor in gtex_tcga.keys():
                    if AC_genes != []:
                        with open('Data/gtex_AC.pkl', 'rb') as archivo:
                            gtex = pickle.load(archivo)
                        get_gtex_info(AC_genes)
                    if DJ_genes != []:
                        with open('Data/gtex_DJ.pkl', 'rb') as archivo:
                            gtex = pickle.load(archivo)
                        get_gtex_info(DJ_genes)
                    if KN_genes != []:
                        with open('Data/gtex_KN.pkl', 'rb') as archivo:
                            gtex = pickle.load(archivo)
                        get_gtex_info(KN_genes)
                    if OR_genes != []:
                        with open('Data/gtex_OR.pkl', 'rb') as archivo:
                            gtex = pickle.load(archivo)
                        get_gtex_info(OR_genes)
                    if ST_genes != []:
                        with open('Data/gtex_ST.pkl', 'rb') as archivo:
                            gtex = pickle.load(archivo)
                        get_gtex_info(ST_genes)
                    if UZ_genes != []:
                        with open('Data/gtex_UZ.pkl', 'rb') as archivo:
                            gtex = pickle.load(archivo)
                        get_gtex_info(UZ_genes)
                data = {'Gene': genes, 'Type':groups, 'Values':values}
                df = pd.DataFrame(data)
                genes = df['Gene'].unique()
                p_values = [] 
                for k in range(len(genes)):
                    gene = genes[k]
                    tumor_data = df[df['Gene'] == gene]
                    types = tumor_data['Type'].unique()
                    for i in range(len(types)):
                        for j in range(i + 1, len(types)):
                            type1 = types[i]
                            type2 = types[j]
                            values1 = tumor_data[tumor_data['Type'] == type1]['Values']
                            values2 = tumor_data[tumor_data['Type'] == type2]['Values']
                            _, p_value = mannwhitneyu(values1, values2)
                            p_values.append(p_value) 
                            final[f'{tumor} median'].append(np.median(values1))
                            final[f'{tumor} sample size'].append(len(values1))
                            final['Control median'].append(np.median(values2))
                            final['Control sample size'].append(len(values2))
                            figure['pvalue'].append(-(math.log10(p_value)))
                            if p_value < 0.001:
                                figure['Legend'].append('Significant')
                                final['Significance'].append('<0.001')
                            elif p_value < 0.01:
                                figure['Legend'].append('Significant')
                                final['Significance'].append('<0.01')
                            elif p_value < 0.05:
                                figure['Legend'].append('Significant')
                                final['Significance'].append('<0.05')
                            else: 
                                figure['Legend'].append('No significant')
                                final['Significance'].append('No significant')
                p_adjusted = stats.false_discovery_control(p_values)
                final['p_adjusted'] = p_adjusted
            else:
                st.error(f'Used FC out of range. Maximum value in this tumor is {max_FC}') 
    table = pd.DataFrame(final)
    if table.size > 0:
        t_data = table.sort_values(by='FC', ascending=False)
        # Show plot and table with the results
        table_data = pd.DataFrame(t_data)
        figure = pd.DataFrame(figure)
        colors = ['rgba(238, 46, 46, 0.83)','rgba(156, 165, 196, 0.95)']
        fig = px.scatter(figure, x='log2(FC)', y='pvalue', color='Legend', color_discrete_sequence=colors,
        custom_data=['Gene'])
        fig.update_traces(hovertemplate='%{customdata}')
        fig.update_layout(title=f'Tumor-associated antigens in {tumor}', title_x=0.35, xaxis_title='log2(FC)', yaxis_title= '-log10(p-value)')
        st.header('Volcano plot', divider='rainbow')
        st.plotly_chart(fig,use_container_width=True)
        st.write(
            f'Interactive volcano plot showing the log2(FC) between \"Primary {tumor} tumor\" and \"Control\" samples and the p-value for all genes with a FC higher than the specified threshold. The name of the gene is displayed by hovering the cursor over the plot. Significant genes (p-value < 0.05) are shown in :red[**red**] while not singiciant genes are shown in :grey[**grey**].'
        )
        st.header('Data table', divider='rainbow')
        st.write(
            f'All relevant data is shown in the table below, including the log2(FC) between \"Primary {tumor} tumor\" and \"Control\" samples and p-value for each gene exceeding the inidcated threshold. You can click in the column names to order the genes according to that column from higher to lower or viceversa. **p-values of 0 are rounded**, you can see the value with all decimals by clicking in the cell.')
        st.dataframe(table_data, hide_index=True)
        st.write(
            'This table can be download in CSV format to be open in Excel or oher softwares compatible with this format.'
        )    