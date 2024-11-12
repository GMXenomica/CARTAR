import math
import pickle
import pandas as pd
import statsmodels.stats.multitest as smm
import streamlit as st
import plotly.express as px
import base64

st.set_page_config(page_title='CARTAR', page_icon='logo.png',layout='wide')
mystyle = '''
    <style>
        p {
            text-align: justify;
        }
    </style>
    '''

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

st.logo('logo_v2.png', icon_image='logo.png')

st.markdown(mystyle, unsafe_allow_html=True)
st.title('Tumor-associated antigens identification')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.write('This tool can be used to identify overexpressed or underexpressed genes in your chosen tumor by inputting the desired fold change ("Primary tumor" compared to "Control" samples). An interactive volcano plot will be shown where the y-axis represent the -log10(p-value), while the x-axis is the log2(fold change). Genes meeting the specified criteria and with an adjusted p-value < 0.05 are highlighted and shown in a table. Tumoral overexpressed genes (with a fold change higher than one) can be used to identify antigens associated with the indicated tumor. The greater the fold change the more associated to the tumor the shown genes will be. You can also identify genes that are underexpressed in tumor samples compared to control ones in case you are looking for candidate genes to design a NOT-gate CAR.')
         
# List with all tumoral options
tumor_options = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
                 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
                 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS']
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
experimental_pm_file = open('Data/HPA_evidence_pm.csv','r')
for line in experimental_pm_file:
    experimental_pm_genes = line.split(',')
# List with expression options
limit_options = ['Above', 'Below']
# Select tumor
tumor = st.selectbox('Select tumor', tumor_options)
if not tumor:
    st.error('Introduce tumor of interest')
# Expander to show abbreviation meaning
with st.expander('Extension of tumor abbreviations\' meaning'):
    for abbreviation, meaning in abbreviations.items():
        st.write(f"**{abbreviation}:** {meaning}")
# Introduce fold change
FC = st.text_input('Enter the Fold Change (FC) value to be used as the threshold for identifying tumor-associated genes:')
FC = FC.replace(',', '.')
if FC and FC.replace('.', '').isdigit():
    if ',' in FC:
        FC_f = float(FC) 
    else:
        FC_f = float(FC)
limit =st.radio('Select whether you are interested in genes above or below the specified threshold.', limit_options)
    
if st.button('Show tumor-associated genes'):
    if FC.replace('.', '').isdigit() and FC != '':
        # Identify genes meeting the stablished threshold
        log2FC = math.log2(FC_f)
        data = pd.read_csv('Data/log2FC_expression.csv')
        final = {'Gene': [], 'log2(FC)': [], 'FC': [], f'{tumor} median': [], f'{tumor} sample size': [],
                'Control median': [], 'Control sample size': [], 'Significance': [], 'p_value': [], 'p_adjusted': [], 'HPA membrane location': []}
        figure = {'Gene': [], 'log2(FC)': [], 'adjusted_p_value': [], 'Legend': []}
        selected_genes = []
        for index, row in data.iterrows():
            gene = row['gene']
            log2_FC = row[tumor]
            figure['Gene'].append(gene)
            figure['log2(FC)'].append(log2_FC)
            if limit == 'Above':
                if row[tumor] >= log2FC:
                    selected_genes.append(gene)
                    final['Gene'].append(gene)
                    final['log2(FC)'].append(log2_FC)
                    final['FC'].append(2**log2_FC)
                    if gene in experimental_pm_genes:
                        final['HPA membrane location'].append('Yes')
                    else:
                        final['HPA membrane location'].append('No')
            else:
                if row[tumor] <= log2FC:
                    selected_genes.append(gene)
                    final['Gene'].append(gene)
                    final['log2(FC)'].append(log2_FC)
                    final['FC'].append(2**log2_FC)
                    if gene in experimental_pm_genes:
                        final['HPA membrane location'].append('Yes')
                    else:
                        final['HPA membrane location'].append('No')
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
            list_p_values = [p_values[gene][tumor] for gene in p_values]
            adjusted_p_values = smm.multipletests(list_p_values, method='fdr_bh')[1]
            n = 0
            for gene in medians.keys():
                p_value = p_values[gene][tumor]
                adjusted_p_value = adjusted_p_values[n]
                if gene in selected_genes:
                    for group in medians[gene][tumor].keys():
                        if group == 'Tumor':
                            final[f'{tumor} median'].append(medians[gene][tumor][group][0])
                            final[f'{tumor} sample size'].append(medians[gene][tumor][group][1])
                        elif group == 'Normal':
                            final['Control median'].append(medians[gene][tumor][group][0])
                            final['Control sample size'].append(medians[gene][tumor][group][1])
                    final['p_value'].append(p_value)
                    final['p_adjusted'].append(adjusted_p_values[n])
                    if adjusted_p_values[n] < 0.001:
                        final['Significance'].append('<0.001')
                    elif adjusted_p_values[n] < 0.01:
                        final['Significance'].append('<0.01')
                    elif adjusted_p_values[n] < 0.05:
                        final['Significance'].append('<0.05')
                    else: 
                        final['Significance'].append('No significant')
                figure['adjusted_p_value'].append(-(math.log10(adjusted_p_values[n])))
                if adjusted_p_values[n] < 0.05 and gene in selected_genes:
                    figure['Legend'].append('Above threshold')
                else:
                    figure['Legend'].append('Below threshold')
                n += 1
        table = pd.DataFrame(final)
        if table.size > 0:
            t_data = table.sort_values(by='FC', ascending=False)
            # Show plot and table with the results
            table_data = pd.DataFrame(t_data)
            figure = pd.DataFrame(figure)
            color1 = 'rgba(156, 165, 196, 0.95)'
            color2 = 'rgba(179, 48, 41, 0.89)'
            fig = px.scatter(figure, x='log2(FC)', y='adjusted_p_value', color='Legend', color_discrete_map={'Above threshold': color2, 'Below threshold':color1,},custom_data=['Gene'])
            fig.update_layout(legend=dict(traceorder='reversed'))
            fig.update_traces(hovertemplate='%{customdata}')
            fig.update_layout(title=f'Tumor-associated antigens in {tumor}', title_x=0.35, xaxis_title='log2(FC)', yaxis_title= '-log10(adjusted p-value)')
            st.header('Volcano plot', divider='rainbow')
            st.plotly_chart(fig,use_container_width=True)
            if limit == 'Above':
                text = f'The interactive volcano plot above visualizes the log2(FC) between "Primary {tumor} tumor" and "Control" samples, along with the corresponding p-values for all genes. Significant genes (adjusted p-value < 0.05) exhibiting a fold change higher than the specified threshold are highlighted in <span style="color:{color2}">red</span>, while other genes appear in <span style="color:{color1}">grey</span>. Hovering over the plot reveals the name of each gene.'
            else:
                text = f'The interactive volcano plot above visualizes the log2(FC) between "Primary {tumor} tumor" and "Control" samples, along with the corresponding p-values for all genes. Significant genes (adjusted p-value < 0.05) exhibiting a fold change lower than the specified threshold are highlighted in <span style="color:{color2}">red</span>, while other genes appear in <span style="color:{color1}">grey</span>. Hovering over the plot reveals the name of each gene.'            
            st.markdown(text, unsafe_allow_html=True)
            st.header('Data table', divider='rainbow')
            if limit == 'Above':
                st.write(
                    f'The table below presents all relevant data, encompassing the log2(FC) between "Primary {tumor} tumor" and "Control" samples, along with the p-value and adjusted p-value for each gene above the specified threshold. You can enhance exploration by clicking on the column names to arrange genes based on that column, either from higher to lower or vice versa. Please note that **the adjusted p-values under 0.001 are rounded to 0**; for the complete decimal value, click on the respective cell.'
                )
            else:
                st.write(
                    f'The table below presents all relevant data, encompassing the log2(FC) between "Primary {tumor} tumor" and "Control" samples, along with the p-value and adjusted p-value for each gene below the specified threshold. You can enhance exploration by clicking on the column names to arrange genes based on that column, either from higher to lower or vice versa. Please note that **p-values under 0.001 are rounded to 0**; for the complete decimal value, click on the respective cell.'
                )  
            st.write(
                'The **HPA (Human Protein Atlas) membrane location column** will be "Yes" if the protein has been experimentally reported to be located in the plasma membrane and "No" when they are only located in the membrane arocding to the Gene Ontology, without havin been experimentally located in the plasma membrane by Human Protien Atlas'
            )
            st.dataframe(table_data, hide_index=True)
            table = table_data.to_csv(encoding='utf-8', index=False)
            b64 = base64.b64encode(table.encode()).decode()
            href = f'<a href="data:file/csv;base64,{b64}" download="table.csv">Download CSV File</a>'
            st.markdown(href, unsafe_allow_html=True)  
    elif FC == '':
        st.error('No FC was introduced') 
    elif FC[0] == '-':
        st.error('FC must be a positive number') 
    else: 
        st.error('Introduce a numerical FC') 
create_footer()
