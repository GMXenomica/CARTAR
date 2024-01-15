import streamlit as st
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
st.title('Metastatic gene expression in SKCM')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.write(
    'This tool will generate the specified plot, comparing expression values of the selected gene in the chosen scale among the "Metastatic," "Primary tumor," and "Normal" samples of Skin Cutaneous Melanoma (SKCM). A table containing all relevant data accompanies this plot. This option is exclusively available for SKCM, owing to the limited number of "Metastatic" samples in other tumors. Data is sourced from the TCGA project, and GTEx expression data has been incorporated into the corresponding control group.'
)
st.set_option('deprecation.showPyplotGlobalUse', False)

scale_options = ['TPM','log2(TPM+1)']
plot_options = ['Boxplot','Violin plot','Dot plot']
gene = st.text_input('Enter gene symbol').upper()
# Identify if indicated gene is present in the data
data = pd.read_csv('Data/log2FC_expression.csv')
if gene == '':
    st.error('Introduce gene symbol')
elif gene != '' and gene not in data['gene'].values:
    st.error(f'{gene} gene symbol not found')
selection = st.radio('Select plot', plot_options)
if selection == 'Boxplot':
    plot = 'Boxplot'
elif selection == 'Violin plot':
    plot = 'Violin plot'
elif selection == 'Dot plot':
    plot = 'Dot plot'
selection2 = st.radio('Select scale', scale_options)
if selection2 == 'TPM':
    scale = 'TPM'
elif selection2 == 'log2(TPM+1)':
    scale = 'log2(TPM+1)'

# Calculate statistical significance and customize plot function
def plot_significance(tumor,y,bottom,top):
    significant_combinations = []
    data = {'Groups compared':[],'Median Group 1':[],'Group 1 sample size':[], 'Median Group 2':[],'Group 2 sample size':[], 'log2(Fold Change)':[],'Significance':[],'p-value':[]}
    # Get the y-axis limits
    y_range = top - bottom
    # Identify groups with statistical difference
    tumor_types = df['Tumor'].unique()
    for i in reversed(range(len(tumor_types))):
        for j in range(i+1, len(tumor_types)):
            tumor1 = tumor_types[i]
            tumor2 = tumor_types[j]
            tumor1_data = df[df['Tumor'] == tumor1]['Values']
            tumor2_data = df[df['Tumor'] == tumor2]['Values']
            _, p_value = mannwhitneyu(tumor1_data, tumor2_data)
            data['Group 1 sample size'].append(len(tumor1_data))
            data['Group 2 sample size'].append(len(tumor2_data))
            data['Groups compared'].append(f'{tumor1} vs {tumor2}')
            data['Median Group 1'].append(np.median(tumor1_data))
            data['Median Group 2'].append(np.median(tumor2_data))
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
                logvalues1 = [log2(value1 + 1) for value1 in tumor1_data]
                logvalues2 = [log2(value2 + 1) for value2 in tumor2_data]
                data['log2(Fold Change)'].append(np.median(logvalues1)-np.median(logvalues2))
            else:
                data['log2(Fold Change)'].append(np.median(tumor1_data)-np.median(tumor2_data))
            if p_value < 0.05:
                if np.median(tumor1_data) > np.median(tumor2_data):
                    k = 0
                else:
                    k = 1
                significant_combinations.append([(i,j),p_value,k])
    # Create table with the results
    table_data = pd.DataFrame(data)
    # Add to the graph the statistical significance bars
    significant_combinations = significant_combinations[::-1]
    for i, significant_combination in enumerate(significant_combinations):
        # Columns corresponding to the datasets of interest
        x1 = significant_combination[0][0]
        x2 = significant_combination[0][1]
        # What level is this bar among the bars above the plot?
        if i == 2:
            level = len(significant_combinations) - i + 1
        else:
            level = len(significant_combinations) - i 
        # Plot the bar
        bar_height = (y_range * 0.07 * level) + top
        bar_tips = bar_height - (y_range * 0.02)
        plt.plot(
            [x1, x1, x2, x2],
            [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k'
        )
        # Significance level
        p = significant_combination[1]
        if p < 0.001:
            sig_symbol = '***'
        elif p < 0.01:
            sig_symbol = '**'
        elif p < 0.05:
            sig_symbol = '*'
        text_height = bar_height + (y_range * 0.0001)
        plt.text((x1 + x2) * 0.5, text_height, sig_symbol, ha='center', va='bottom', color='black')
    # Customize the plot
    plt.title(f'{gene} expression comparison between SKCM conditions', y=1.03)
    if scale == 'TPM':
        plt.ylabel(f'{gene} expression in TPM')
    else: 
        plt.ylabel(f'{gene} expression in log2(TPM+1)')
    plt.xlabel('SKCM group')
    plt.xticks(rotation=45, ha='right')
    plt.subplots_adjust(left=0.067, bottom=0.155, right=0.968, top=0.91)
    # Show the graph and table
    st.header(plot, divider='rainbow')
    st.pyplot()
    st.write(
        f'The above figure illustrates the {plot} for {gene} expression in {scale} across "Metastatic," "Primary tumor," and "Control" SKCM samples, comparing the expression between these groups. Statistical significance is denoted for each SKCM group (***: p_value < 0.001, **: p_value < 0.01, *: p_value < 0.05).'
    )
    st.header('Data table', divider='rainbow')
    st.write(
        f'All relevant information is presented in the table below, encompassing critical aspects such as the log2(Fold Change) for each comparison. Computed as the log2(TPM+1) median of SKCM Group 1 expression minus the log2(TPM+1) median of SKCM Group 2 expression. For ease of exploration, you can click on the column names to arrange the rows based on the selected column, either in ascending or descending order. Please note that **p-values under 0.001 are rounded to 0**; for the complete decimal value, click on the respective cell.'
    )
    st.dataframe(table_data, hide_index=True)
    st.write(
        'This table can be downloaded in CSV format.'
    ) 

if st.button(f'Create {plot}'):
    if gene != '' and gene in data['gene'].values: 
        # Identify GTEx tissue sample corresponding to control group of tumor
        gtex_tcga = {'SKCM':'Skin'}
        with open('Data/SKCM.pkl', 'rb') as archivo:
            SKCM = pickle.load(archivo)
        # Get requested information
        groups = [] # Gruops of tumor (Metastatic, Primary or Normal)
        values = [] # Expression values
        for group in SKCM[gene]['SKCM'].keys():
            for value in SKCM[gene]['SKCM'][group]:
                groups.append(group)
                if scale == 'log2(TPM+1)':
                    value = log2(value+1)
                values.append(value)
        data = {'Tumor': groups, 'Values':values}
        df = pd.DataFrame(data)
        group_order = ['Metastatic', 'Primary', 'Normal']
        df['Tumor'] = pd.Categorical(df['Tumor'], categories=group_order, ordered=True)
        df = df.sort_values(by=['Tumor'])
        # Create the boxplot
        if plot == 'Boxplot':
            plt.figure()
            sns.boxplot(data=df, x='Tumor', y='Values', hue='Tumor', legend=False, showfliers=False, palette={'Primary': 'lightseagreen', 'Normal': 'tan', 'Metastatic':'grey'})
            xmin, xmax, ymin, ymax = plt.axis()
            # Statistical significant differences and customize the plot
            plot_significance('SKCM',0,ymin,ymax)
        # Create the violin plot
        if plot == 'Violin plot':
            sns.violinplot(x='Tumor', y='Values', hue='Tumor', data=df, inner='quartile', density_norm='width',palette={'Primary': 'lightseagreen', 'Normal': 'tan', 'Metastatic':'grey'}, legend=False)
            xmin, xmax, ymin, ymax = plt.axis()
            # Statistical significant differences and customize the plot
            plot_significance('SKCM',1,ymin,ymax)        
        # Create the dotplot
        if plot == 'Dot plot':
            sns.stripplot(x='Tumor', y='Values', jitter=True, data=data, hue='Tumor', size=4, palette={'Primary': 'lightseagreen', 'Normal': 'tan', 'Metastatic':'grey'})
            plt.xlim(-0.5, 2.5)
            xmin, xmax, ymin, ymax = plt.axis()
            # Calculate the medians for each group
            medians = df.groupby('Tumor', observed=False)['Values'].median()
            # Add a horizontal line for each median within the corresponding group
            n = 0.25
            for tumor, median in medians.items():
                x_start = xmin + n
                x_end = x_start + 0.5
                n += 1
                plt.plot(
                    [x_start, x_end],
                    [median, median], lw=1, c='k', zorder=10000
                )
            # Statistical significant differences and customize the plot
            plot_significance('SKCM',1,ymin,ymax)        
    elif gene not in data['gene'].values:
        st.error(f'{gene} gene symbol not found')
    else:
        st.error('No gene symbol was introduced')        