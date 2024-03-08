# RDoC_RNAseq
Bulk-RNAseq analyses on post-mortem brain tissue samples from BA32 and BA16 using a classical cases-control group comparison design and a dimensional regression design using RDoC scores.

Abstract
=

Alzheimer Disease (AD) is primarily defined by cognitive impairment. However, neuropsychiatric symptoms (NPS) in AD can profoundly shape clinical trajectories. For example, aggressive behavior, psychosis, apathy, and depression strain healthcare and support provider resources, and lead to accelerated cognitive decline and earlier death. It is thus critical to better understand NPS in AD and their underlying biological mechanisms. Human postmortem brain tissue research is pivotal for a deeper knowledge on the molecular mechanism of AD and NPS, however, postmortem studies typically focus on categorical phenotypes such as the presence of a disease. To overcome limitations of categorical study designs and the lack of extended molecular studies on NPS in AD, we build on our prior work and implemented a proven natural language processing (NLP)-based approach to profile NIMH Research Domain Criteria (RDoC) dimensions in n=101 brain donors across a broad spectrum of AD neuropathology. We then associated the derived dimensional phenotypes with bulk RNAseq data from the insula (Brodmann Area 16) and dorsal anterior cingulate cortex (Brodmann Area 32), two brain regions important for motivation, decision making, and affective processing among others. Our findings show distinct transcriptional profiles associated with RDoC criteria dimensions and provide insight into genes and gene pathways implicated in NPS in AD. The transcriptional profiles associated with RDoC criteria dimensions augment the mechanistic insight into AD derived from prior studies and provide novel approaches for much needed interventions.

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

