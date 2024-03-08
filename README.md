# RDoC_RNAseq
=
Bulk-RNAseq analyses on post-mortem brain tissue samples from BA32 and BA16 using a classical cases-control design and a dimensional design using RDoC scores.

Abstract
=

Alzheimer Disease (AD) is primarily defined by cognitive impairment. However, neuropsychiatric symptoms (NPS) in AD can profoundly shape clinical trajectories. For example, aggressive behavior, psychosis, apathy, and depression strain healthcare and support provider resources, and lead to accelerated cognitive decline and earlier death. It is thus critical to better understand NPS in AD and their underlying biological mechanisms. Human postmortem brain tissue research is pivotal for a deeper knowledge on the molecular mechanism of AD and NPS, however, postmortem studies typically focus on categorical phenotypes such as the presence of a disease. To overcome limitations of categorical study designs and the lack of extended molecular studies on NPS in AD, we build on our prior work and implemented a proven natural language processing (NLP)-based approach to profile NIMH Research Domain Criteria (RDoC) dimensions in n=101 brain donors across a broad spectrum of AD neuropathology. We then associated the derived dimensional phenotypes with bulk RNAseq data from the insula (Brodmann Area 16) and dorsal anterior cingulate cortex (Brodmann Area 32), two brain regions important for motivation, decision making, and affective processing among others. Our findings show distinct transcriptional profiles associated with RDoC criteria dimensions and provide insight into genes and gene pathways implicated in NPS in AD. The transcriptional profiles associated with RDoC criteria dimensions augment the mechanistic insight into AD derived from prior studies and provide novel approaches for much needed interventions.1

Data and analysis
=

Raw data: RNAseq fastq files can be accessed through GEO accession GSE261050
Processed data: The data that is used for the scripts is uploaded to the github page, the counts table is stored in the tximport-counts3.csv, metadata is stored in the AD_RNAseq_metadata3.xlsx, RDoC scores of five domains for each individual is stored in the RDoC.csv.
Analysis: The main analysis normalization, differential expression analysis and enrichment analysis is performed in the RNAseq_analayis.R and the deconvolution of the RNAseq data with MuSiC is performed in the Deconvolution_MuSiC.R. Additionally, for the reviewer we also explored other Deseq2 models and the code is available for those models in RNAseq_analysis_multiple_DEseq2_models.R. The functions that are used throughout the analysis can be found in the functions.R script
