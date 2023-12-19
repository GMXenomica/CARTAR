import math
import pickle
import pandas as pd
from scipy import stats
import streamlit as st
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
st.title('Tumor-associated antigens identification')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.write(
    'This tool facilitates the identification of tumor-associated antigens specific to the indicated tumor. Simply input your desired Fold Change (FC) threshold between "Primary tumor" and "Control" samples. The "Primary tumor" data is sourced from TCGA, while the "Control" data is a combination of TCGA and GTEX. The tool produces a volcano plot highlighting statistical significant genes above or below the specified FC threshold, accompanied by a detailed table providing expression data for these genes and their corresponding p-values.'
)
st.set_option('deprecation.showPyplotGlobalUse', False)


# List with all tumoral options
tumor_options = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
                 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
                 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS']
# List with expression options
limit_options = ['Above', 'Below']
# Select tumor
tumor = st.selectbox('Select tumor', tumor_options)
if not tumor:
    st.error('Introduce tumor of interest')
# Introduce fold change
FC = st.text_input('Enter the Fold Change (FC) value to be used as the threshold for identifying tumor-associated genes:')
if FC:
    FC = float(FC.replace(',', '.')) if ',' in FC else float(FC)
limit =st.radio('Select whether you are interested in genes above or below the specified threshold.', limit_options)

if st.button('Show tumor-associated genes'):
    # Identify genes meeting the stablished threshold
    log2FC = math.log2(FC)
    data = pd.read_csv('Data/log2FC_expression.csv')
    final = {'Gene': [], 'log2(FC)': [], 'FC': [], f'{tumor} median': [], f'{tumor} sample size': [],
             'Control median': [], 'Control sample size': [], 'Significance': [], 'p_adjusted': []}
    figure = {'Gene': [], 'log2(FC)': [], 'pvalue': [], 'Legend': []}
    selected_genes = []
    for index, row in data.iterrows():
        gene = row['gene']
        log2_FC = row[tumor]
        FC = 2 ** log2_FC
        figure['Gene'].append(gene)
        figure['log2(FC)'].append(log2_FC)
        if limit == 'Above':
            if row[tumor] >= log2FC:
                selected_genes.append(gene)
                final['Gene'].append(gene)
                final['log2(FC)'].append(log2_FC)
                final['FC'].append(FC)
        else:
            if row[tumor] <= log2FC:
                selected_genes.append(gene)
                final['Gene'].append(gene)
                final['log2(FC)'].append(log2_FC)
                final['FC'].append(FC)
    if not selected_genes:
        max_FC = 2 ** data[tumor].max()
        min_FC = 2 ** data[tumor].min()
        st.error(f'Fold Change (FC) used is beyond the valid range. The maximum allowable value in this tumor is {max_FC} and the minimum is {min_FC}.')
    else:
        # Determine median and p-value of differential expression
        with open('Data/median.pkl', 'rb') as archivo:
            medians = pickle.load(archivo)
        with open('Data/p_value.pkl', 'rb') as archivo:
            p_values = pickle.load(archivo)
        stored_p_values = []
        for gene in medians.keys():
            p_value = p_values[gene][tumor]
            if gene in selected_genes:
                for group in medians[gene][tumor].keys():
                    if group == 'Tumor':
                        final[f'{tumor} median'].append(medians[gene][tumor][group][0])
                        final[f'{tumor} sample size'].append(medians[gene][tumor][group][1])
                    elif group == 'Normal':
                        final['Control median'].append(medians[gene][tumor][group][0])
                        final['Control sample size'].append(medians[gene][tumor][group][1])
                stored_p_values.append(p_value)
                p_adjusted = stats.false_discovery_control(stored_p_values)
                final['p_adjusted'] = p_adjusted
                if p_value < 0.001:
                    final['Significance'].append('<0.001')
                elif p_value < 0.01:
                    final['Significance'].append('<0.01')
                elif p_value < 0.05:
                    final['Significance'].append('<0.05')
                else: 
                    final['Significance'].append('No significant')
            figure['pvalue'].append(-(math.log10(p_value)))
            if p_value < 0.05 and gene in selected_genes:
                figure['Legend'].append('Above threshold')
            else:
                figure['Legend'].append('Below threshold')
    table = pd.DataFrame(final)
    if table.size > 0:
        t_data = table.sort_values(by='FC', ascending=False)
        # Show plot and table with the results
        table_data = pd.DataFrame(t_data)
        figure = pd.DataFrame(figure)
        color1 = 'rgba(156, 165, 196, 0.95)'
        color2 = 'rgba(179, 48, 41, 0.89)'
        fig = px.scatter(figure, x='log2(FC)', y='pvalue', color='Legend', color_discrete_map={'Above threshold': color2, 'Below threshold':color1,},custom_data=['Gene'])
        fig.update_layout(legend=dict(traceorder='reversed'))
        fig.update_traces(hovertemplate='%{customdata}')
        fig.update_layout(title=f'Tumor-associated antigens in {tumor}', title_x=0.35, xaxis_title='log2(FC)', yaxis_title= '-log10(p-value)')
        st.header('Volcano plot', divider='rainbow')
        st.plotly_chart(fig,use_container_width=True)
        if limit == 'Above':
            text = f'The interactive volcano plot above visualizes the log2(FC) between "Primary {tumor} tumor" and "Control" samples, along with the corresponding p-values for all genes. Significant genes (p-value < 0.05) exhibiting a fold change higher than the specified threshold are highlighted in <span style="color:{color2}">red</span>, while other genes appear in <span style="color:{color1}">grey</span>. Hovering over the plot reveals the name of each gene.'
        else:
            text = f'The interactive volcano plot above visualizes the log2(FC) between "Primary {tumor} tumor" and "Control" samples, along with the corresponding p-values for all genes. Significant genes (p-value < 0.05) exhibiting a fold change lower than the specified threshold are highlighted in <span style="color:{color2}">red</span>, while other genes appear in <span style="color:{color1}">grey</span>. Hovering over the plot reveals the name of each gene.'            
        st.markdown(text, unsafe_allow_html=True)
        st.header('Data table', divider='rainbow')
        if limit == 'Above':
            st.write(
                f'The table below presents all relevant data, encompassing the log2(FC) between "Primary {tumor} tumor" and "Control" samples, along with the p-value for each gene above the specified threshold. You can enhance exploration by clicking on the column names to arrange genes based on that column, either from higher to lower or vice versa. Please note that **the adjusted p-values under 0.001 are rounded to 0**; for the complete decimal value, click on the respective cell.'
            )
        else:
            st.write(
                f'The table below presents all relevant data, encompassing the log2(FC) between "Primary {tumor} tumor" and "Control" samples, along with the p-value for each gene below the specified threshold. You can enhance exploration by clicking on the column names to arrange genes based on that column, either from higher to lower or vice versa. Please note that **p-values under 0.001 are rounded to 0**; for the complete decimal value, click on the respective cell.'
            )           
        st.dataframe(table_data, hide_index=True)
        st.write(
            'This table can be downloaded in CSV format.'
        )    