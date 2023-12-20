import streamlit as st
import base64

st.set_page_config(
    page_title='CARTAR',
    page_icon='logo.png',
    layout='wide'
)

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

mystyle = '''
    <style>
        p {
            text-align: justify;
        }
    </style>
    '''

st.markdown(mystyle, unsafe_allow_html=True)

st.write('# Welcome to CARTAR')
st.markdown('### Introduction')
st.markdown('The CARTAR web page offers a suite of tools designed to facilitate the identification of potential targets for Chimeric Antigen Receptor (CAR) therapies. Leveraging expression data from The Cancer Genome Atlas (TCGA) project ([tcga_RSEM_gene_tpm](https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz)) and the Genotype-Tissue Expression (GTEx) project ([gtex_RSEM_gene_tpm](https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_RSEM_gene_tpm.gz)), our platform focuses on pinpointing tumor-associated antigens, ensuring target selectivity, and assessing specificity.')
st.markdown('CARTAR is based on RNA sequencing expression data of 10,522 samples from the TCGA project and 7,858 samples from the GTEx project, using a standard pre-processing pipeline.') 
st.markdown('This website is free and open to all users and there is no login requirement. The code used for raw files preprocessing is available at [GitHub](https://github.com/mighg/CARTAR).')
st.markdown('### Tumor-associated antigens')
st.markdown('This tool can be used to identify overexpressed or underexpressed genes in your chosen tumor by inputting the desired fold change ("Primary tumor" compared to "Control" samples). An interactive volcano plot will be shown where the y-axis represent the -log10(p-value), while the x-axis is the log2(fold change). Genes meeting the specified criteria and with an adjusted p-value < 0.05 are highlighted and shown in a table. Tumoral overexpressed genes (with a fold change higher than one) can be used to identify antigens associated with the indicated tumor. The greater the fold change the more associated to the tumor the shown genes will be. You can also identify genes that are underexpressed in tumor samples compared to control ones in case you are looking for candidate genes to design a NOT-gate CAR.')
st.markdown('### Tumor expression change')
st.markdown('This tool can be used to explore fold change expression values ("Primary tumor" compared to "Control" samples) for a gene or gene set of interest in the desired tumors. This will provide a table with the fold cahnge value for each indicated gene in each selected tumor. This will allow to get preliminary information for candidate target genes to check if they are overexpressed in certain tumor and if it can be used to treat more than one tumor.')
st.markdown('### Tumor median expression')
st.markdown('This tool can be used to easily visualize a barplot showcasing the median expression value of a selected gene across primary tumor, metastatic (if available), and control samples across tumors of interest. Besides, median expression values and sample sizesof each group are reported in table format. This will allow to get information about a candidate target gene identify with the**Tumor-associated antigens tool** to check its expression in "Metastatic" samples if available and whether it can be used to treat more than one tumor.')
st.markdown('### Tumor gene expression')
st.markdown('This tool can be used to generate boxplots, violin plots, or dot plots for the expression values of a gene of interest in primary tumor and control samples of desired tumor groups. Besides, median expression values, sample sizes of each group, and statistical significance of differential expression (obtained by Mann–Whitney U test) between "Primary tumor" and "Control" samples are reported in table format. This provide more detailed information than the **Tumor median expression tool** while no information about "Metastatic" samples is included due to its reduced sample size. In case of Skin Cutaneous Melanoma (SKCM) the "Metastatic" sample size (N=366) you canthe comparisson between the three samples group in the **Metastatic gene expression tool**.')
st.markdown('### Tissue gene expression')
st.markdown('This tool can be used to create boxplots, violin plots, or dot plots for the expression values of a gene of interest across all GTEx tissues. When a tumor group is specified, the expression values of "Primary tumor" samples are included in the plot and compared to the expression values in each GTEx tissue. Besides, median expression values, sample sizes of each group, and statistical significance of differential expression when a tumor is introduce are reported in table format. This is critical to assess the specificity of candidate target genes. The CAR therapy will reacognise the target antigen in all expressing cells and it is important to ensure that no vital tissue cells are destroyed.')
st.markdown('### Metastatic gene expression')
st.markdown('This tool can be used to generate boxplots, violin plots, or dot plots for the expression values of a gene of interest for the "Primary tumor", "Metastatic", and "Control" samples of **Skin Cutaneous Melanoma (SKCM)**. SKCM is the only TCGA tumor group with sufficient "Metastatic" sample size (N=366) to get statistical significance of differential expression. Besides, median expression values, sample sizes for each group, and statistical significance of differential expression between sample groups are reported in table format.')
st.markdown('### Logic-gated CAR')
st.markdown('This tool can be used to explore the correlation between two genes expression values in "Primary tumor" and "Control" samples of a specified tumor to assess the potential of their combination for a logic-gated CAR therapy (OR-gate, AND-gate, NOT-gate, or IF-BETTER-gate).')
st.markdown('- **OR-gate:** effector cell activation upon binding to either one of two antigens. This strategy can be used when either one of the antigens is highly expressed in "Primary tumor" samples, while their expression is low in "Control" samples.')
st.markdown('- **AND-gate:** effector cell activation upon binding to both antigens simultaneously. This strategy can be used when both antigens are expressed in "Primary tumor" samples, while none or only one of the antigens is expressed in "Control" samples.')
st.markdown('- **NOT-gate:** effector cell activation upon binding to one antigen while the second antigen is not recognised. This strategy can be used when the antigen triggering the cytotoxic activity is present is certain healthy tissues and you wnat to increase selectivity by including an antigen thar respress this activity and is only expressed in the healthy tissues. To identify this antigens you must go to **Tissue gene expression tool**.')
st.markdown('- **IF-BETTER-gate:** effector cell activation upon binding to high CAR target expression. Activation is not triggered upon binding to low CAR targer expression, unless the second antigen is present.')
st.markdown('This tool can be used to Identify cancer cell lines for testing selected candidate targets. An interactive barplot displays cell lines meeting specified target expression criteria, with median expression and additional information provided in table format. This tool utilizes expression values from the Cancer Cell Line Encyclopedia (CCLE) ([OmicsExpressionProteinCodingGenesTPMLogp1.csv](https://depmap.org/portal/download/all/?release=DepMap+Public+23Q4&file=OmicsExpressionProteinCodingGenesTPMLogp1.csv#:~:text=DepMap%20Public%2023Q4-,OmicsExpressionProteinCodingGenesTPMLogp1.csv,-Gene%20expression%20TPM)).')
st.markdown('Cancer cell lines to test the cytotoxic activity of a CAR therapy targeted against a candidate gene can be selected looking for thecell lines of the tumoral lineage of interest with high expression of the target gene. On the other hand control cell lines can be selected by looking for those with low expression of the target gene.')   
st.markdown('### Contact us')
st.markdown('If any question or request about CARTAR, please contact us: miguel.gamarra@usc.es')
