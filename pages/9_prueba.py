from itertools import count
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
from math import log2
import numpy as np
import seaborn as sns
import requests
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

st.markdown(mystyle, unsafe_allow_html=True)
st.title('Gene expression across tumors')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.write(
    'This tool will generate the indicated plot comparing expression values of the chosen gene in the selected scale comparing \"Primary tumor\" to \"Normal\" samples and a table with all relevant data. Data is obtained from TCGA project and GTEX expression data has also been included in the correponding control group. If no tumor is indicated it will report the expression data of the indicated gene in all tumors.'
)
st.set_option('deprecation.showPyplotGlobalUse', False)

tumor_options  = ['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS']
scale_options = ['TPM','log2(TPM+1)']
plot_options = ['Boxplot','Violinplot','Dotplot']

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
selection = st.radio('Choose plot', plot_options)
if selection == 'Boxplot':
    plot = 'Boxplot'
elif selection == 'Violinplot':
    plot = 'Violinplot'
elif selection == 'Dotplot':
    plot = 'Dotplot'
selection2 = st.radio('Choose scale', scale_options)
st.info('log2(TPM+1) scale is recommended for better visualization due to outliers')
if selection2 == 'TPM':
    scale = 'TPM'
elif selection2 == 'log2(TPM+1)':
    scale = 'log2(TPM+1)'

# Generate the data for the plot function
def plot_data(data):
    df = pd.DataFrame(data)
    group_order = ['Tumor', 'Normal']
    df['Type'] = pd.Categorical(df['Type'], categories=group_order, ordered=True)
    df = df.sort_values(by=['Tumor', 'Type'])
    return(df)

# Calculate statistical significance and customize plot function
def plot_significance(tumors,y):
    # Creat table data
    data = {'Tumor':[],'Tumor median':[], 'Tumor sample size':[], 'Control median':[], 'Control sample size':[],'log2(Fold Change)':[],'Significance':[],'p-value':[]}
    # Create list to store percentile90 of each group
    percentile90 = []
    # Compare statistical difference between groups of a tumor type
    tumors = df['Tumor'].unique()
    for k in range(len(tumors)):
        tumor = tumors[k]
        tumor_data = df[df['Tumor'] == tumor]
        types = tumor_data['Type'].unique()
        positions = np.arange(len(categories))
        for i in range(len(types)):
            for j in range(i + 1, len(types)):
                type1 = types[i]
                type2 = types[j]
                values1 = tumor_data[tumor_data['Type'] == type1]['Values']
                percentile_90_1 = np.percentile(values1, 90)
                data['Tumor sample size'].append(len(values1))
                values2 = tumor_data[tumor_data['Type'] == type2]['Values']
                percentile_90_2 = np.percentile(values2, 90)
                if percentile_90_2 > percentile_90_1:
                    percentile_90_1 = percentile_90_2
                percentile90.append(percentile_90_1)
                data['Control sample size'].append(len(values2))
                _, p_value = mannwhitneyu(values1, values2)
                data['Tumor'].append(tumor)
                data['Tumor median'].append(np.median(values1))
                data['Control median'].append(np.median(values2))
                data['p-value'].append(p_value)
                if p_value < 0.001:
                    data['Significance'].append('<0.001')
                elif p_value < 0.01:
                    data['Significance'].append('<0.01')
                elif p_value < 0.05:
                    data['Significance'].append('<0.05')
                else: 
                    data['Significance'].append('No significant')
                if scale == 'TPM':
                    logvalues1 = [log2(value1 + 1) for value1 in values1]
                    logvalues2 = [log2(value2 + 1) for value2 in values2]
                    data['log2(Fold Change)'].append(np.median(logvalues1)-np.median(logvalues2))
                else:
                    data['log2(Fold Change)'].append(np.median(values1)-np.median(values2))
                # Identify significant differences
                if p_value < 0.05:
                    color = 'red' if np.median(values1) > np.median(values2) else 'green'
                # Add * representation of significance
                if p_value < 0.001:
                    if y == 1:
                        plt.text(positions[k], percentile_90_1, '***', color=color, ha='center', fontsize=8)
                    else:
                        plt.text(positions[k], max(df['Values']), '***', color=color, ha='center', fontsize=8)
                elif p_value < 0.01:
                    if y == 1:
                        plt.text(positions[k], percentile_90_1, '**', color=color, ha='center', fontsize=8)
                    else:
                        plt.text(positions[k], max(df['Values']), '**', color=color, ha='center', fontsize=8)
                elif p_value < 0.05:
                    if y == 1:
                        plt.text(positions[k], percentile_90_1, '*', color=color, ha='center', fontsize=8)
                    else:
                        plt.text(positions[k], max(df['Values']), '*', color=color, ha='center', fontsize=8)
    # Create table with the results
    table_data = pd.DataFrame(data)
    # Customize the plot
    plt.title(f'{gene} expression difference between tumoral conditions', y=1.03)
    if y == 1:
        if scale == 'TPM':
            plt.ylim(0,max(percentile90)+20)
        else:
            plt.ylim(0,max(percentile90)+0.5)
    else:
        plt.ylim(0, max(df['Values']))
    if scale == 'TPM':
        plt.ylabel(f'{gene} expression in TPM')
    else: 
        plt.ylabel(f'{gene} expression in log2(TPM+1)')
    plt.xlabel('Tumor')
    plt.xticks(rotation=45, ha='right')
    plt.subplots_adjust(left=0.067, bottom=0.135, right=0.968, top=0.91)
    # Show the graph and table
    st.header(plot, divider='rainbow')
    st.pyplot()
    st.write(
        f'The above figure shows the {plot} for {gene} expression in {scale} across selected tumors and compares the expression between \"Primary tumor\" and \"Normal\" samples. Statistical significance is shown between each indicated tumor and its corresponding control (***: p_value < 0.001, **: p_value < 0.01, *: p_value < 0.05). It will be shown in :red[**red**] when the gene is **overexpressed** in \"Primary tumor\" samples compared to \"Normal\" samples, and in :green[**green**] when **underexpressed**.'
    )
    st.header('Data table', divider='rainbow')
    st.write(
        f'All relevant data is shown in the table below, including the log2(Fold Change) for each comparison calculated as log2(TPM+1) median of primary tumor expression - log2(TPM+1) median of control expression, and the p-value of each comparison. You can click in the column names to order the tumors according to that column from higher to lower or viceversa. By clicking in a cell you can see the value with all the decimals.'
    )
    st.dataframe(table_data, hide_index=True)
    st.write(
        'This table can be download in CSV format to be open in Excel or oher softwares compatible with this format.'
    )   
    
# Dictionary with the tumor:GTEX healthy tissue correspondance
gtex_tcga = {'ACC':'Adrenal Gland','BLCA':'Bladder','BRCA':'Breast','CESC':'Cervix Uteri', 'COAD':'Colon', 'DLBC':'Blood',
                    'ESCA':'Esophagus','GBM':'Brain','KICH':'Kidney','KIRC':'Kidney','KIRP':'Kidney','LAML':'Bone Marrow','LGG':'Brain',
                    'LIHC':'Liver','LUAD':'Lung','LUSC':'Lung','OV':'Ovary','PAAD':'Pancreas','PRAD':'Prostate','READ':'Colon','SKCM':'Skin',
                    'STAD':'Stomach','TGCT':'Testis','THCA':'Thyroid','THYM':'Blood','UCEC':'Uterus','UCS':'Uterus'}

if st.button(f'Create {plot}'):
    if gene != '':
        # If gene and tumor abreviation in data
        if gene in data['gene'].values:  
            # Load dictionary with tcga data  
            if 'A' <= gene[0] <= 'C'
                data = 'Data/tcga_AC.pkl'
            elif 'D' <= gene[0] <= 'J':   
                data = 'Data/tcga_DJ.pkl'
            elif 'K' <= gene[0] <= 'N':
                data = 'Data/tcga_KN.pkl'
            elif 'O' <= gene[0] <= 'R':
                data = 'Data/tcga_OR.pkl'
            elif 'S' <= gene[0] <= 'T':
                data = 'Data/tcga_ST.pkl'
            elif 'U' <= gene[0] <= 'Z':
                data = 'Data/tcga_UZ.pkl'
            with open(data, 'rb') as archivo:
                tcga = pickle.load(archivo)
            # Get requested information
            categories = [] # List with tumor types
            groups = [] # Gruops of tumor (Primary or Normal)
            values = [] # Expression values
            for tumor in tumors:
                for group in tcga[gene][tumor].keys():
                    for value in tcga[gene][tumor][group]:
                        categories.append(tumor)
                        groups.append(group)
                        values.append(value)
            # Load dictionary with gtex data
            x = 1
            for tumor in tumors:
                if tumor in gtex_tcga.keys():
                    if x:
                        with open('data/tissue_combined.pkl', 'rb') as archivo:
                            gtex = pickle.load(archivo)
                        x = 0
                    group = 'Normal'
                    for value in gtex[gene][gtex_tcga[tumor]]:
                        categories.append(gtex_tcga[tumor])
                        groups.append(group)
                        values.append(value)                    
            
            data = {'Tumor': categories, 'Type':groups, 'Values':values}
            # Create the boxplot
            if plot == 'Boxplot':
                df = plot_data(data)
                plt.figure()
                sns.set(style='whitegrid')
                sns.boxplot(x='Tumor', y='Values', hue='Type', data=df, palette={'Tumor': 'lightseagreen', 'Normal': 'tan'}, flierprops={'marker': '.'},whis=(10, 90),showfliers=False)
                xmin, xmax, ymin, ymax = plt.axis()
                # Calculate statistical significance and customize the plot
                plot_significance(tumors,1)
            # Create the violin plot
            if plot == 'Violinplot':
                df = plot_data(data)
                plt.figure()
                sns.violinplot(x='Tumor', y='Values', hue='Type', data=df, split=True, inner='quartile', density_norm='width',palette={'Tumor': 'lightseagreen', 'Normal': 'tan', 'Metastatic':'grey'})
                xmin, xmax, ymin, ymax = plt.axis()
                # Calculate statistical significance and customize plot
                plot_significance(tumors,0)         
            # Create the dotplot
            if plot == 'Dotplot':
                data = {'Tumor': categories, 'Type':groups, 'Values':values}
                df = plot_data(data)
                plt.figure()
                sns.stripplot(x='Tumor', y='Values', jitter=True, hue='Type', data=data, size=4, palette={'Tumor': 'lightseagreen', 'Normal': 'tan', 'Metastatic':'grey'})
                xmin, xmax, ymin, ymax = plt.axis()
                # Calculate statistical significance and customize the plot
                plot_significance(tumors,0)   
