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
st.title('Welcome to CARTAR')
st.markdown('<style>div.block-container{padding-top:1rem;}</style>',unsafe_allow_html=True)
st.markdown('### Introduction')
st.markdown('The CARTAR web page offers a suite of tools designed to facilitate the identification of potential targets for Chimeric Antigen Receptor (CAR) therapies. Leveraging expression data from The Cancer Genome Atlas (TCGA) project ([tcga_RSEM_gene_tpm](https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz)) and the Genotype-Tissue Expression (GTEx) project ([gtex_RSEM_gene_tpm](https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_RSEM_gene_tpm.gz)), our platform focuses on pinpointing tumor-associated antigens, ensuring target selectivity, and assessing specificity.')
st.markdown('CARTAR is based on RNA sequencing expression data of 10,522 samples from the TCGA project and 7,858 samples from the GTEx project, using a standard pre-processing pipeline.') 
st.markdown('Raw data was meticulously processed to isolate genes located in the plasma membrane (GO:0005886). These membrane-bound proteins serve as prime candidates for CAR therapy targeting. Additionally, the dataset was annotated with tumor types and conditions (primary tumor, metastatic, or control), with GTEx data serving as control samples where applicable.')
st.markdown('This website is free and open to all users and there is no login requirement. The code used for raw files preprocessing is available at [GitHub](https://github.com/mighg/CARTAR).')
st.markdown('### Cookie policy')
st.markdown('Please be aware that this web tool utilizes cookies established by Streamlit for enhanced user experience and functionality. By continuing to use the tool, you agree to the usage of these cookies.')
st.markdown('### CARTAR tools')
st.markdown('- Tumor-associated antigens: identify tumor-associated antigens')
st.markdown('- Tumor expression change: explore fold change expression values of a gene or set of genes accross specified tumors')
st.markdown('- Tumor median expression: visualize median expression values of a candidate gene across primary tumors, metastatic (if available), and control samples')
st.markdown('- Tumor gene expression: analyze expression values of a specific gene across primary tumor and control samples of specified tumors')
st.markdown('- Tissue gene expression: study off-tumor gene expression to assess specificity of candidate target genes')
st.markdown('- Metastatic gene expression: analyze expression values of a specific gene across primary tumor, metastatic and control samples on Skin Cutaneous Melanoma (SKCM)')
st.markdown('- Logic-gated CAR: explore the correlation between the expression values of two genes in primary tumor and control samples of a specified tumor for the design of dual- targeting CAR therapy')
st.markdown('- Cell line selector: identify cancer cell lines with desired expression values of target gene')
st.markdown('### Contact us')
st.markdown('For any questions or requests about CARTAR, please contact us: miguel.gamarra@usc.es')
