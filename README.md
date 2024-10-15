# RDoC_RNAseq
Bulk-RNAseq analyses on post-mortem brain tissue samples from BA32 and BA16 using a dimensional regression design using RDoC scores.

Abstract
=

INTRODUCTION: Neuropsychiatric symptoms are common in people with Alzheimer’s disease (AD) across all severity stages. Their heterogeneous presentation and variable temporal association with cognitive decline suggest shared and distinct biological mechanisms. We hypothesized that specific patterns of gene expression associate with distinct NIMH Research Domain Criteria (RDoC) domains in AD.
METHODS: Post-mortem bulk RNAseq on the insula and anterior cingulate cortex from 60 brain donors representing the spectrum of canonical AD neuropathology combined with natural language processing approaches based on the RDoC Clinical Domains. 
RESULTS: Distinct sets of >100 genes (pFDR<0.05) were specifically associated with at least one clinical domain (Cognitive, Social, Negative, Positive, Arousal). In addition, dysregulation of immune response pathways was shared across domains and brain regions.
DISCUSSION: Our findings provide evidence for distinct transcriptional profiles associated with RDoC domains suggesting that each dimension is characterized by specific sets of genes providing insight into the underlying mechanisms.


Data and analysis
=

# Raw data: 
RNAseq fastq files can be accessed through GEO accession GSE261050
# Processed data: 
tximport-counts3.csv: raw gene count matrix of all samples  
AD_RNAseq_metadata3.xlsx: metadata of all samples  
RDoC.csv: RDoC scores of five domains for each individual, for detailed information about each domain, please refer to NIMH website: https://www.nimh.nih.gov/research/research-funded-by-nimh/rdoc/definitions-of-the-rdoc-domains-and-constructs.  
# Analysis: 
RNAseq_RDoC.R: main analysis including data preprocessing, PCA analysis, SVA analysis, differential gene analysis based on a case-ctrl group comparison design and dimensional regression analysis over each domain of RDoC scores using R pachage ImpulseDE2.  
plot_PCA_deseq2_modified.R: additional customized function required in PCA analysis part

