import math
import pickle
import pandas as pd
from scipy import stats
import streamlit as st
import plotly.express as px

# Lista de opciones de tumor
tumor_options = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
                 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
                 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS']

# Interfaz para seleccionar el tumor
tumor = st.selectbox('Choose tumor', tumor_options)
if not tumor:
    st.error('Select tumor of interest')

# Interfaz para introducir Fold Change
FC = st.text_input('Introduce Fold Change value used as threshold to identify tumor-associated genes:')
if FC:
    FC = float(FC.replace(',', '.')) if ',' in FC else float(FC)
    if 1 > FC or FC < 1.8:
        st.info('Please take into account that this could take some time')
        if FC < 1:
            st.error('Fold changes under 1 are not recommended')

# Botón para mostrar genes asociados al tumor
if st.button('Show tumor-associated genes'):
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
        if row[tumor] >= log2FC:
            selected_genes.append(gene)
            final['Gene'].append(gene)
            final['log2(FC)'].append(log2_FC)
            final['FC'].append(FC)
    if not selected_genes:
        max_FC = 2 ** data[tumor].max()
        st.error(f'Used FC out of range. Maximum value in this tumor is {max_FC}')
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
        fig = px.scatter(figure, x='log2(FC)', y='pvalue', color='Legend', color_discrete_map={'Below threshold':'rgba(156, 165, 196, 0.95)','Above threshold': 'rgba(179, 48, 41, 0.89)'},
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