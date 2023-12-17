import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
from math import log2
import numpy as np
import seaborn as sns
import pickle
from scipy.stats import mannwhitneyu

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
st.title('Gene expression across GTEX tissues')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.write(
    'This tool will generate the indicated plot comparing expression values of the chosen gene in the selected scale of the indicated tumor (data obtained from TCGA project) to GTEX tissue samples and a table with all relevant data. If no tumor is indicated it will report the expression data of the indicated gene in all GTEX tissues.'
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
tumor = st.selectbox('Choose tumor (optional)', options=[None] + tumor_options)
if not tumor:
    tumor = ''
    st.error('If no tumor is indicated this will only show gene expression across GTEX tissues')
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

# Create a list with all gtex tissues
tissues = ['Blood','Blood Vessel','Brain','Thyroid','Pancreas','Muscle','Lung','Skin','Colon','Nerve','Adipose Tissue','Ovary','Heart','Breast','Pituitary','Testis','Vagina','Esophagus','Small Intestine','Spleen','Adrenal Gland','Stomach','Uterus','Liver','Bone Marrow','Salivary Gland','Prostate','Kidney','Bladder','Fallopian Tube','Cervix Uteri']
# Function for statistical analysis of significance of tumoral group vs gtex data and customize plot function
def statistics(y):
    if tumor != '':
        if scale == 'TPM':
            n = 15
        else:
            n = 0.5
        plt.axvline(x=0.5, color='red', linestyle='--') # Add a line to separte the tumor group fromt the GTEX data
        # Creat table data
        data = {'GTEX tissue':[],f'{tumor} median':[],f'{tumor} sample size':[], 'Tissue median':[],'Tissue sample size':[], 'log2(Fold Change)':[],'Significance':[],'p-value':[]}
        # Statistical analysis of significance of tumoral group vs gtex data
        for i in range(len(groups)):
            if groups[i] != tumor + ' Tumor':
                percentile_90 = np.percentile(groups_data[i], 90)
                data[f'{tumor} sample size'].append(len(groups_data[0]))
                data['Tissue sample size'].append(len(groups_data[i]))
                _, p_value = mannwhitneyu(groups_data[i], groups_data[0])
                data['GTEX tissue'].append(groups[i])
                data[f'{tumor} median'].append(np.median(groups_data[0]))
                data['Tissue median'].append(np.median(groups_data[i]))
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
                    logvalues1 = [log2(value1 + 1) for value1 in groups_data[0]]
                    logvalues2 = [log2(value2 + 1) for value2 in groups_data[i]]
                    data['log2(Fold Change)'].append(np.median(logvalues1)-np.median(logvalues2))
                else:
                    data['log2(Fold Change)'].append(np.median(groups_data[0])-np.median(groups_data[i]))
                # Identify significant differences
                if p_value < 0.05:
                    color = 'red' if np.median(groups_data[i]< np.median(groups_data[0])) else 'green'
                # Add * representation of significance
                if p_value < 0.001:
                    if y == 1:
                        plt.text(positions[i], percentile_90, '***', color=color,ha='center', fontsize=8)
                    else:
                        plt.text(positions[i], max(groups_data[i]) + n, '***', color=color,ha='center',fontsize=8)
                elif p_value < 0.01:
                    if y == 1:
                        plt.text(positions[i], percentile_90, '**', color=color,ha='center', fontsize=8)
                    else:
                        plt.text(positions[i], max(groups_data[i]) + n, '***', color=color,ha='center',fontsize=8)
                elif p_value < 0.05:
                    if y == 1:
                        plt.text(positions[i], percentile_90, '*', color=color,ha='center', fontsize=8)
                    else:
                        plt.text(positions[i], max(groups_data[i]) + n, '***', color=color,ha='center',fontsize=8)
    else:
        data = {'GTEX tissue':[],'Tissue median':[],'Tissue sample size':[]}
        for i in range(len(groups)):
            data['GTEX tissue'].append(groups[i])
            data['Tissue median'].append(np.median(groups_data[i]))
            data['Tissue sample size'].append(len(groups_data[i]))
    # Create table with the results
    table_data = pd.DataFrame(data)  
    # Customize plot 
    if y == 0:
        plt.ylim(0, max(df['Values']) + 2)
    plt.xlabel('Tissues')
    if tumor != '':
        plt.title(f'{gene} expression comparison between {tumor} and GTEX data', y=1.03)
    else: 
        plt.title(f'{gene} expression in GTEX data', y=1.03)
    if scale == 'TPM':
        plt.ylabel(f'{gene} expression in TPM')
    else: 
        plt.ylabel(f'{gene} expression in log2(TPM+1)')
    plt.xticks(rotation=45,ha='right')
    plt.subplots_adjust(left=0.067, bottom=0.2, right=0.968, top=0.91)
    # Show the graph and table
    st.header(plot, divider='rainbow')
    st.pyplot()
    if tumor:
        if plot != 'Dotplot':
            st.write(
                f'The above figure shows the {plot} for {gene} expression in {scale} across {tumor} and all GTEX tissues. Statistical significance is shown between the {tumor} and each GTEX tissues (***: p_value < 0.001, **: p_value < 0.01, *: p_value < 0.05). It will be shown in :red[**red**] when the gene is **overexpressed** in {tumor} compared to the GTEX tissue, and in :green[**green**] when **underexpressed**.'
            )
        else:
            st.write(
                f'The above figure shows the {plot} for {gene} expression in {scale} across {tumor} and all GTEX tissues. **Median** of each group is indicated with the black bar. Statistical significance is shown between the {tumor} and each GTEX tissues (***: p_value < 0.001, **: p_value < 0.01, *: p_value < 0.05). It will be shown in :red[**red**] when the gene is **overexpressed** in {tumor} compared to the GTEX tissue, and in :green[**green**] when **underexpressed**.'
            )
    else:
        if plot != 'Dotplot':
            st.write(
                f'The above figure shows the {plot} for {gene} expression in {scale} across all GTEX tissues. No statistical significance is shown as no tumor was indicated.'
            )
        else:
            st.write(
                f'The above figure shows the {plot} for {gene} expression in {scale} across all GTEX tissues. **Median** of each group is indicated with the black bar. Statistical significance is shown between the {tumor} and each GTEX tissues (***: p_value < 0.001, **: p_value < 0.01, *: p_value < 0.05). It will be shown in :red[**red**] when the gene is **overexpressed** in {tumor} compared to the GTEX tissue, and in :green[**green**] when **underexpressed**.'
            )
    st.header('Data table', divider='rainbow')
    if not tumor:
        st.write(
            f'All relevant data is shown in the table below, as no tumor was chosen there is no information comparing the expresion to any tumor. In case you are interested in knowing the log2(FoldChange) of a given tumor expression compared to GTEX tissues you must indicated the tumor of interest. You can click in the column names to order the tissues according to that column from higher to lower or viceversa. By clicking in a cell you can see the value with all the decimals.'
        )
    else:
        st.write(
            f'All relevant data is shown in the table below, including the log2(Fold Change) for each comparison calculated as log2(TPM+1) median of {tumor} expression - log2(TPM+1) median of GTEX tissue expression, and the p-value of each comparison. You can click in the column names to order the tissues according to that column from higher to lower or viceversa. By clicking in a cell you can see the value with all the decimals.'
        )
    st.dataframe(table_data, hide_index=True)
    st.write(
        'This table can be download in CSV format to be open in Excel or oher softwares compatible with this format.'
    )    

if st.button(f'Create {plot}'):
    if gene != '':
        groups = [] # List with the names of the tissue
        values = [] # List with the exoression values
        if tumor != None:
            # Load dictionary with tcga data  
            if 'A' <= gene[0] <= 'C':
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
            for value in tcga[gene][tumor]['Tumor']:
                groups.append(tumor + ' Tumor')
                if scale == 'log2(TPM+1)':
                    value = log2(value+1)
                values.append(value)
        # Load dictionary with gtex data  
        if 'A' <= gene[0] <= 'C':
            data = 'Data/gtex_AC.pkl'
        elif 'D' <= gene[0] <= 'J':   
            ata = 'Data/gtex_DJ.pkl'
        elif 'K' <= gene[0] <= 'N':
            data = 'Data/gtex_KN.pkl'
        elif 'O' <= gene[0] <= 'R':
            data = 'Data/gtex_OR.pkl'
         elif 'S' <= gene[0] <= 'T':
            data = 'Data/gtex_ST.pkl'
         elif 'U' <= gene[0] <= 'Z':
            data = 'Data/gtex_UZ.pkl'
        with open(data, 'rb') as archivo:
            gtex = pickle.load(archivo)
        for tissue in gtex[gene].keys():
            for value in gtex[gene][tissue]:
                groups.append(tissue)
                if scale == 'log2(TPM+1)':
                    value = log2(value+1)
                values.append(value)
        # Create the dicitionary with all the data
        data = {'Groups':groups, 'Values':values} 
        # Select tissues and order them alphabetically 
        indexes = np.unique(data['Groups'], return_index=True)[1]
        groups = [data['Groups'][index] for index in sorted(indexes)]
        df = pd.DataFrame(data)
        # Create a list with the values for each group
        groups_data = []
        for group in groups:
            group_data = np.array(data['Values'])[np.array(data['Groups']) == group]
            groups_data.append(group_data)
        # Calculate the position of each group in the graph
        positions = np.arange(len(groups))
        # Create the boxplot
        if plot == 'Boxplot':
            sns.boxplot(data=df, x='Groups', y='Values', hue='Groups', palette='Spectral', whis=(10, 90) ,legend=False, showfliers=False)
            xmin, xmax, ymin, ymax = plt.axis()
            # Statistical analysis of significance of tumoral group vs gtex data and customize plot function
            statistics(1)
        # Create a violin plot
        if plot == 'Violinplot':
            plt.figure(figsize=(12, 6))
            sns.violinplot(data=df, x='Groups', y='Values', hue='Groups', palette='Spectral', legend=False)
            xmin, xmax, ymin, ymax = plt.axis()
            # Statistical analysis of significance of tumoral group vs gtex data and customize plot function
            statistics(0)
        # Create the dotplot
        if plot == 'Dotplot':
            plt.figure(figsize=(12, 6))
            sns.stripplot(data=df, x='Groups', y='Values', jitter=True, hue='Groups', palette='Spectral', legend=False, size=4)
            xmin, xmax, ymin, ymax = plt.axis()
            # Calculate the medians for each group
            medians = df.groupby('Groups', observed=False, sort=False)['Values'].median()
            # Add a horizontal line for each median within the corresponding group
            n = 0.25
            for tissue, median in medians.items():
                x_start = xmin + n
                x_end = x_start + 0.5
                n += 1
                plt.plot(
                    [x_start, x_end],
                    [median, median], lw=1, c='k', zorder=10000
                )
            # Statistical analysis of significance of tumoral group vs gtex data
            statistics(0)