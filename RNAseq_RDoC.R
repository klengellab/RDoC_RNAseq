rm(list =ls())
library("corrplot")
library('apeglm')
library(DESeq2)
library(dplyr)
library("factoextra")
library(ggplot2)
library(ggrepel)
library(ImpulseDE2)
.libPaths()
library("pheatmap", lib.loc=.libPaths()[1])
library("readxl")
library('rtracklayer')
library(sva)
library(tibble)
library(UpSetR)
options(ggrepel.max.overlaps=Inf)

# read gene count matrix
file_path = 'E:/project/sabina/result/'
file_name = 'tximport-counts3.csv'
data<-read.table(file = paste0(file_path,file_name),
                 header = TRUE,sep = ",")
table(duplicated(data$gene))

# read metadata and processing
file_path = 'E:/project/sabina/result/'
file_name2 = 'AD_RNAseq_metadata3.xlsx'
metadata0<-read_excel(path = paste0(file_path,file_name2),
                              col_names = TRUE) 
metadata0<-data.frame(metadata0)

for (i in 1:nrow(metadata0)) {
  metadata0[i,'Sample_id'] = strsplit(metadata0[i,2],'_')[[1]][1]
  metadata0[i,'Area'] = strsplit(metadata0[i,2],'_')[[1]][2]
  metadata0[i,'Plate_number'] = strsplit(
    metadata0[i,'Special.Comments'],' ')[[1]][2]
}

metadata0[,2]<-gsub("_",".",metadata0[,2])
rownames(metadata0)<-metadata0[,2]

metadata0 <- data.frame(metadata0) %>%
  filter(Sample.Name. %in% colnames(data[,2:ncol(data)]))

metadata0 <- metadata0 %>%
  rename(Sample_name = Sample.Name., Group = Dx, 
         Age = age, Sex = sex,
         Extraction_Batch = RNA_extraction_batch,
         Sequence_Batch = Plate_number,
         RIN = RIN_genewiz)

metadata_final <- metadata0[,c("Sample_name","Sample_id",
                               "Group","Area","Age","Sex","B.B",
                               "RIN","PMI",
                               "Sequence_Batch","Extraction_Batch")]

!colSums(is.na(metadata_final))

metadata_final$Area<-as.factor(metadata_final$Area)
metadata_final$Group<-as.factor(metadata_final$Group)
metadata_final$Sequence_Batch<-as.factor(metadata_final$Sequence_Batch)
metadata_final$Sex<-as.factor(metadata_final$Sex)
metadata_final$Extraction_Batch<-as.factor(metadata_final$Extraction_Batch)
metadata_final$B.B<-as.factor(metadata_final$B.B)
metadata_final$RIN <- as.numeric(metadata_final$RIN)

metadata_final[metadata_final$Age<60,'Age2'] <- '50-60'
metadata_final[(metadata_final$Age>=60 & metadata_final$Age<70),
               'Age2'] <- '60-70'
metadata_final[(metadata_final$Age>=70 & metadata_final$Age<80),
               'Age2'] <- '70-80'
metadata_final[(metadata_final$Age>=80 & metadata_final$Age<90),
               'Age2'] <- '80-90'
metadata_final[(metadata_final$Age>=90 & metadata_final$Age<100),
               'Age2'] <- '90-100'
metadata_final$Age2<-as.factor(metadata_final$Age2)

metadata_final <- metadata_final[order(
  metadata_final$Group,metadata_final$Area,metadata_final$Sex,
  metadata_final$Age2),]
summary(metadata_final)


# read RDoC domains data and merge with metadata
file_path = 'E:/project/sabina/original_data/'
file_name3 = 'RDoC.csv'
rdoc<-read.table(file = paste0(file_path,file_name3),
                 header = TRUE,sep = ",")
metadata_final_rdoc <- metadata_final[
  metadata_final$Sample_id %in% rdoc$AN,]

# here 'rdoc' column stands for cognition
metadata_final_rdoc$rdoc <- 0
for (i in c(1:nrow(metadata_final_rdoc))) {
  metadata_final_rdoc[i,'rdoc'] = rdoc[rdoc$AN==metadata_final_rdoc[i,'Sample_id'],'cognitive']
}
metadata_final_rdoc$positive <- 0
for (i in c(1:nrow(metadata_final_rdoc))) {
  metadata_final_rdoc[i,'positive'] = rdoc[rdoc$AN==metadata_final_rdoc[i,'Sample_id'],'positive']
}
metadata_final_rdoc$negative <- 0
for (i in c(1:nrow(metadata_final_rdoc))) {
  metadata_final_rdoc[i,'negative'] = rdoc[rdoc$AN==metadata_final_rdoc[i,'Sample_id'],'negative']
}
metadata_final_rdoc$arousal <- 0
for (i in c(1:nrow(metadata_final_rdoc))) {
  metadata_final_rdoc[i,'arousal'] = rdoc[rdoc$AN==metadata_final_rdoc[i,'Sample_id'],'arousal_regulatory']
}
metadata_final_rdoc$social <- 0
for (i in c(1:nrow(metadata_final_rdoc))) {
  metadata_final_rdoc[i,'social'] = rdoc[rdoc$AN==metadata_final_rdoc[i,'Sample_id'],'social']
}

summary(metadata_final_rdoc)

rdoc2 <- cut(metadata_final_rdoc$rdoc,breaks = seq(0,0.8,by=0.1),
             labels = c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),
             right = F)
positive2 <- cut(metadata_final_rdoc$positive,breaks = seq(0,0.7,by=0.1),
                 labels = c('0.1','0.2','0.3','0.4','0.5','0.6','0.7'),
                 right = F)
negative2 <- cut(metadata_final_rdoc$negative,breaks = seq(0,0.7,by=0.1),
                 labels = c('0.1','0.2','0.3','0.4','0.5','0.6','0.7'),
                 right = F)
arousal2 <- cut(metadata_final_rdoc$arousal,breaks = seq(0,0.6,by=0.1),
                labels = c('0.1','0.2','0.3','0.4','0.5','0.6'),
                right = F)
social2 <- cut(metadata_final_rdoc$social,breaks = seq(0,0.8,by=0.1),
               labels = c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),
               right = F)

metadata_final_rdoc<-cbind(metadata_final_rdoc,rdoc2,positive2,negative2,arousal2,social2)

metadata_final_rdoc$rdoc2 <- as.factor(metadata_final_rdoc$rdoc2)
metadata_final_rdoc$positive2 <- as.factor(metadata_final_rdoc$positive2)
metadata_final_rdoc$negative2 <- as.factor(metadata_final_rdoc$negative2)
metadata_final_rdoc$arousal2 <- as.factor(metadata_final_rdoc$arousal2)
metadata_final_rdoc$social2 <- as.factor(metadata_final_rdoc$social2)

metadata_final_rdoc <- metadata_final_rdoc %>%
  rename(cognition = rdoc)

metadata_final_rdoc$interaction1 <- factor(
  paste(metadata_final_rdoc$Group,metadata_final_rdoc$Area,sep = '-'))
metadata_final_rdoc$interaction2 <- factor(
  paste(metadata_final_rdoc$Group,metadata_final_rdoc$Sex,sep = '-'))
metadata_final_rdoc$interaction3 <- factor(
  paste(metadata_final_rdoc$Group,metadata_final_rdoc$Area,metadata_final_rdoc$Sex,sep = '-'))

summary(metadata_final_rdoc)

# gene annotation
gtf_data = import('E:/project/sabina/result/Homo_sapiens.GRCh38.105.gtf')
gtf_data = as.data.frame(gtf_data)
gtf_data2 <- data.frame(gtf_data) %>%
  dplyr::filter(gene_id %in% data$gene)
gtf_data3 <- gtf_data2[gtf_data2$type=='gene',
                       c("gene_id",'gene_name')]
rm(gtf_data,gtf_data2)
data <- data.frame(data) %>%
  dplyr::filter(gene %in% gtf_data3$gene_id)
df_original <- merge(gtf_data3,data, by.x="gene_id", by.y="gene")
rownames(df_original) = df_original[,1]

# filter low-expressing genes
Ex<-df_original[,3:ncol(df_original)]
idx <- rowSums((Ex)>=10) >= 10
idx2 <- rowSums((Ex)>0) > 1

data2 <- Ex[idx,]
Ex2 <- Ex[idx2,]

data_final<-data2[,order(match(colnames(data2),
                               metadata_final$Sample_name))]
data_final0<-Ex2[,order(match(colnames(Ex2),
                              metadata_final$Sample_name))]

data_final_rdoc <-data_final[,order(match(
  colnames(data_final),metadata_final_rdoc$Sample_name))]
data_final_rdoc <- data_final_rdoc[
  ,1:nrow(metadata_final_rdoc)]

data_final0_rdoc <-data_final0[,order(match(
  colnames(data_final0),metadata_final_rdoc$Sample_name))]
data_final0_rdoc <- data_final0_rdoc[
  ,1:nrow(metadata_final_rdoc)]

# PCA analysis
dds <- DESeqDataSetFromMatrix(data_final_rdoc,
                                   colData = metadata_final_rdoc,
                                   design= ~ Group)
dds$Group<-relevel(dds$Group,ref = 'control')

dds0 <- dds
dds <- estimateSizeFactors(dds)
vsd <- vst(object=dds,blind=FALSE)

source("E:/project/sabina/result/plot_PCA_deseq2_modified.R")
colData(vsd)
pcaData <- plotPCA.deseq2_modified(vsd, 
                                   ntop = 5000,
                                   intgroup=c("Group","Area",'Sex','B.B',
                                              'RIN','PMI','Age',
                                              'Sequence_Batch', 'Extraction_Batch'), 
                                   returnData=TRUE)
pca_PCs = pcaData[[1]]
pca_origin = pcaData[[2]]
percentVar <- round(100 * attr(pcaData[[1]], "percentVar"))
res.ind<-get_pca_ind(pca_origin)
fviz_eig(pca_origin, addlabels = TRUE, 
                     ylim = c(0, 30),ncp = 10)

ggplot(pca_PCs, aes(PC1, PC2, color=Group, shape=Group)) +
  geom_point(size = 2) +
  labs(
    x = paste0("PC1"), 
    y = paste0("PC2"))

ggplot(pca_PCs, aes(PC1, PC2, color=Area, shape=Area)) +
  geom_point(size = 2) +
  labs(
    x = paste0("PC1"), 
    y = paste0("PC2"))

ggplot(pca_PCs, aes(PC1, PC2, color=Sex, shape=Sex)) +
  geom_point(size = 2) +
  labs(
    x = paste0("PC1"), 
    y = paste0("PC2"))

ggplot(pca_PCs, aes(PC1, PC2, color=Sequence_Batch, shape=Sequence_Batch)) +
  geom_point(size = 2) +
  labs(
    x = paste0("PC1"), 
    y = paste0("PC2"))

# correlation analysis
metadata_final0 <- metadata_final_rdoc
colnames(metadata_final0)
identical(rownames(pca_PCs),metadata_final0$Sample_name)

cor_mat <- cbind(pca_PCs[,1:5],
                 metadata_final0[,c('Age','PMI','B.B','RIN',
                                        'Group',"Area",'Sex',
                                        'Sequence_Batch','Extraction_Batch'
                 )])
cor_mat<-sapply(cor_mat, as.double)
rownames(cor_mat)<-rownames(pca_PCs)
all(!is.na(cor_mat))

M <- cor(cor_mat)
res1 <- cor.mtest(cor_mat, conf.level = .95) 
col1 <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
test <- res1$p
test2 <- as.vector(test)
test3 <- p.adjust(test2,method = 'fdr')
test3 <- as.matrix(test3)
test3 <- matrix(test3,nrow = nrow(test),byrow = T)
rownames(test3) <- rownames(test)
colnames(test3) <- colnames(test)
res1$p <- test3

corrplot(M, col = col1(20),type = 'upper',
         tl.pos = 'tp', tl.srt = 45, tl.col = 'black',
         p.mat = res1$p, sig.level = .05, insig = 'blank')
corrplot(M, add = TRUE,method = "number",col = col1(20),type = 'lower',
         tl.pos = 'n', cl.pos = 'n',
         p.mat = res1$p, sig.level = .05, insig = 'blank')

# SVA analysis: iterating number of surrogate variables until representing the best of residual PCs other
# than those PCs already correlated with known variables
dds <- estimateSizeFactors(dds)
dat  <- counts(dds, normalized = TRUE)
all(!is.na(dat))
sum(colnames(dds)==rownames(metadata_final_rdoc)) == nrow(metadata_final_rdoc)

mod  <- model.matrix(~ interaction1 + Sex + Sequence_Batch + Age + RIN, 
                     colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
svseq <- sva::svaseq(dat, mod, mod0, n.sv = 2)

# redo correlation analysis
metadata_final0 <- metadata_final_rdoc
colnames(metadata_final0)
identical(rownames(pca_PCs),metadata_final0$Sample_name)

cor_mat <- cbind(pca_PCs[,1:5],
                 metadata_final0[,c('Age','PMI','B.B','RIN',
                                    'Group',"Area",'Sex',
                                    'Sequence_Batch','Extraction_Batch'
                 )])
cor_mat<-sapply(cor_mat, as.double)
rownames(cor_mat)<-rownames(pca_PCs)
all(!is.na(cor_mat))

test <- as.data.frame(svseq$sv)
for (i in c(1:ncol(test))) {
  colnames(test)[i] <- paste0('SV',i)
}

cor_mat <- data.frame(cor_mat)
cor_mat <- cbind(cor_mat,test)

M <- cor(cor_mat)
res1 <- cor.mtest(cor_mat, conf.level = .95) 
col1 <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
test <- res1$p
test2 <- as.vector(test)
test3 <- p.adjust(test2,method = 'fdr')
test3 <- as.matrix(test3)
test3 <- matrix(test3,nrow = nrow(test),byrow = T)
rownames(test3) <- rownames(test)
colnames(test3) <- colnames(test)
res1$p <- test3

corrplot(M, col = col1(20),type = 'upper',
         tl.pos = 'tp', tl.srt = 45, tl.col = 'black',
         p.mat = res1$p, sig.level = .05, insig = 'blank')
corrplot(M, add = TRUE,method = "number",col = col1(20),type = 'lower',
         tl.pos = 'n', cl.pos = 'n',
         p.mat = res1$p, sig.level = .05, insig = 'blank')

# finally merge SVs with metadata
test <- as.data.frame(svseq$sv)
for (i in c(1:ncol(test))) {
  colnames(test)[i] <- paste0('SV',i)
}
metadata_final_rdoc <- cbind(metadata_final_rdoc,test)

# DE analysis
ddssva2 <- dds0
ddssva2$SV1 <- svseq$sv[,1]
ddssva2$SV2 <- svseq$sv[,2]

ddssva2$interaction1 <- factor(paste(ddssva2$Group,ddssva2$Area,
                                     sep = '_'))
ddssva2$interaction2 <- factor(paste(ddssva2$Group,ddssva2$Sex,
                                     sep = '_'))
ddssva2$interaction3 <- factor(paste(ddssva2$Group,ddssva2$Area,
                                     ddssva2$Sex,sep = '_'))
colData(ddssva2)

design(ddssva2)<- ~ SV1+SV2+Sex+Sequence_Batch+Age+RIN+interaction1

table(ddssva2$interaction1)

# control_insula as reference
ddssva2$interaction1<-relevel(ddssva2$interaction1,
                              ref = 'control_insula')
# control_BA32 as reference
# ddssva2$interaction1<-relevel(ddssva2$interaction1,
#                               ref = 'control_BA32')

ddssva2 <- DESeq(ddssva2)
resultsNames(ddssva2)

contrast_fmt = c("interaction1", "AD_insula", "control_insula")
# contrast_fmt = c("interaction1", "AD_BA32", "control_BA32")

res <- results(ddssva2,
               contrast = contrast_fmt,
               pAdjustMethod = "fdr",
               alpha = 0.05,
               lfcThreshold=0)
resultsNames(ddssva2)

coef_fmt = "interaction1_AD_insula_vs_control_insula"
# coef_fmt = "interaction1_AD_BA32_vs_control_BA32"

res <- lfcShrink(ddssva2, 
                 coef=coef_fmt,
                 type = 'apeglm')
# check for NAs
!colSums(is.na(res))
colnames(res)
res <- res[ !is.na(res$padj), ]
res <- res[ !is.na(res$pvalue), ]
res <- res[ !is.na(res$log2FoldChange), ]
!colSums(is.na(res))

summary(res,alpha = 0.05)

# gene annotation
test2 <- data.frame(res) %>%
  rownames_to_column(var = "gene")
res_test <- merge(test2, gtf_data3,
                  by.x="gene", by.y="gene_id")
rownames(res_test) <- res_test$gene_id

naOrDup <- is.na(res_test$gene_name) | duplicated(res_test$gene_name)
res_test$annotation <- ifelse(naOrDup, res_test$gene, res_test$gene_name)
table(duplicated(res_test$annotation))

# identify DEGs
padj_cutoff = 0.05; FC_cutoff = 0

res_tbl <- res_test %>%
  as_tibble()
sig_res <- filter(res_tbl, (padj < padj_cutoff & abs(log2FoldChange) >= FC_cutoff)) %>%
  arrange(-log2FoldChange)
sig_res <- data.frame(sig_res)

res_insula_AD_ctrl_rdoc_samples <- res
dds_insula_AD_ctrl_rdoc_samples <- ddssva2
sig_res_insula_AD_ctrl_rdoc_samples <- sig_res
# res_BA32_AD_ctrl_rdoc_samples <- res
# dds_BA32_AD_ctrl_rdoc_samples <- ddssva2
# sig_res_BA32_AD_ctrl_rdoc_samples <- sig_res

# volcano plot
plotVocalno_DIY <- function(res_test,padj_cutoff,FC_cutoff,plot_range){
  tt2<-as.data.frame(res_test)
  tt2$log10PValue<- -log10(tt2$padj)
  tt2$Group ="NS"
  tt2$Group[which((tt2[,'log2FoldChange']>= FC_cutoff) & (tt2[,'padj'] < padj_cutoff) )]="Up"
  tt2$Group[which((tt2[,'log2FoldChange']<= -FC_cutoff) & (tt2[,'padj'] < padj_cutoff) )]="Down"
  tt2$Group <- factor(tt2$Group, levels = c("Up","Down","NS"))
  
  down <- tt2[which(tt2$padj < padj_cutoff & tt2$log2FoldChange <= -FC_cutoff),] %>%
    arrange(padj)
  top_genes_neg <- head(down,10)
  up <- tt2[which(tt2$padj < padj_cutoff & tt2$log2FoldChange >= FC_cutoff),] %>%
    arrange(padj)
  top_genes_pos <- head(up,10)
  
  DEG_20 <- c(top_genes_pos$annotation,
               top_genes_neg$annotation)
  tt2$label <- ifelse(tt2$annotation %in% DEG_20,
                      tt2$annotation,"")
  
  ggplot(data = tt2,
         aes(x = log2FoldChange,
             y = log10PValue,
             color = Group)) +
    geom_point(alpha=0.8, size=1.5) +
    xlim(c(-plot_range,plot_range)) +
    geom_vline(xintercept = -FC_cutoff, 
               linetype = "dashed", colour = "black") +
    geom_vline(xintercept = FC_cutoff,
               linetype = "dashed", colour = "black") +
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed", colour = "black") + 
    geom_label_repel(data = tt2, aes(label = label))
}

plotVocalno_DIY(res_test,
                padj_cutoff = padj_cutoff,
                FC_cutoff = FC_cutoff,
                plot_range=6)

# heatmap
vsd <- vst(object=dds_insula_AD_ctrl_rdoc_samples,blind=FALSE)
DE_Genes <- sig_res_insula_AD_ctrl_rdoc_samples$gene
mat  <- assay(vsd)[ DE_Genes, ]

metadata_final2 <- metadata_final_rdoc[order(metadata_final_rdoc$Group,
                                             metadata_final_rdoc$Area),]

mat2<-mat[,order(match(colnames(mat), metadata_final2$Sample_name))]
mat2 <- mat2[,1:nrow(metadata_final2)]

annotation_col = data.frame(
  metadata_final2[, c('Group','Area')]
)
rownames(annotation_col) <- metadata_final2$Sample_name

colors <- colorRampPalette(c("#3399FF",'white',"#FF3333"))(15)
pheatmap(mat2,
         color = colors,
         show_rownames = F,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA,
         na_col = 'white',
         fontsize = 12,
         scale = "row",
         annotation_col = annotation_col) 


# ImpulseDE2 analysis
# replace insula with BA32 to regress RDoC of BA32 samples

# inherit estimates of dispersion parameters from DESeq2
head(dispersions(dds_insula_AD_ctrl_rdoc_samples))
dispersions_deseq <- dispersions(dds_insula_AD_ctrl_rdoc_samples)
names(dispersions_deseq) <- rownames(assay(dds_insula_AD_ctrl_rdoc_samples))
head(dispersions_deseq)

metadata_rdoc_combined_insula <- metadata_final_rdoc[metadata_final_rdoc$Area=='insula',]
!colSums(is.na(metadata_rdoc_combined_insula))

test <- counts(dds_insula_AD_ctrl_rdoc_samples)
counts <- test[,colnames(test) %in% metadata_rdoc_combined_insula$Sample_name]

metadata_rdoc_combined_insula['Sample'] <- metadata_rdoc_combined_insula['Sample_name']
rownames(metadata_rdoc_combined_insula) <- metadata_rdoc_combined_insula$Sample
# ImpulseDE2:case-only model
metadata_rdoc_combined_insula['Condition'] <- 'case'
metadata_rdoc_combined_insula['Batch'] <- metadata_rdoc_combined_insula['Sequence_Batch']

metadata_rdoc_combined_insula$Sample <- as.character(metadata_rdoc_combined_insula$Sample)
metadata_rdoc_combined_insula$Condition <- as.character(metadata_rdoc_combined_insula$Condition)
metadata_rdoc_combined_insula$Batch <- as.character(metadata_rdoc_combined_insula$Batch)
metadata_rdoc_combined_insula$Sex <- as.character(metadata_rdoc_combined_insula$Sex)

metadata_rdoc_heat_all <- metadata_rdoc_combined_insula

# dimensional regression over cognition domain
metadata_rdoc_combined_insula['Time'] <- metadata_rdoc_combined_insula['rdoc2']
metadata_rdoc_combined_insula$Time <- as.numeric(metadata_rdoc_combined_insula$Time)
metadata_rdoc_combined_insula2 <- metadata_rdoc_combined_insula[
  ,c('Sample','Condition','Time','Batch','Sex')]

objectImpulseDE2_rdoc_insula_cognition_rdoc_samples <- runImpulseDE2(
  matCountData    = counts, 
  dfAnnotation    = metadata_rdoc_combined_insula2,
  boolCaseCtrl    = FALSE,
  vecDispersionsExternal = dispersions_deseq,
  vecConfounders  = c("Batch","Sex"),
  boolIdentifyTransients = TRUE)

# dimensional regression over positive domain
metadata_rdoc_combined_insula['Time'] <- metadata_rdoc_combined_insula['positive2']
metadata_rdoc_combined_insula$Time <- as.numeric(metadata_rdoc_combined_insula$Time)
metadata_rdoc_combined_insula2 <- metadata_rdoc_combined_insula[
  ,c('Sample','Condition','Time','Batch','Sex')]

objectImpulseDE2_rdoc_insula_positive_rdoc_samples <- runImpulseDE2(
  matCountData    = counts, 
  dfAnnotation    = metadata_rdoc_combined_insula2,
  boolCaseCtrl    = FALSE,
  vecDispersionsExternal = dispersions_deseq,
  vecConfounders  = c("Batch","Sex"),
  boolIdentifyTransients = TRUE)

# dimensional regression over negative domain
metadata_rdoc_combined_insula['Time'] <- metadata_rdoc_combined_insula['negative2']
metadata_rdoc_combined_insula$Time <- as.numeric(metadata_rdoc_combined_insula$Time)
metadata_rdoc_combined_insula2 <- metadata_rdoc_combined_insula[
  ,c('Sample','Condition','Time','Batch','Sex')]

objectImpulseDE2_rdoc_insula_negative_rdoc_samples <- runImpulseDE2(
  matCountData    = counts, 
  dfAnnotation    = metadata_rdoc_combined_insula2,
  boolCaseCtrl    = FALSE,
  vecDispersionsExternal = dispersions_deseq,
  vecConfounders  = c("Batch","Sex"),
  boolIdentifyTransients = TRUE)

# dimensional regression over arousal domain
metadata_rdoc_combined_insula['Time'] <- metadata_rdoc_combined_insula['arousal2']
metadata_rdoc_combined_insula$Time <- as.numeric(metadata_rdoc_combined_insula$Time)
metadata_rdoc_combined_insula2 <- metadata_rdoc_combined_insula[
  ,c('Sample','Condition','Time','Batch','Sex')]

objectImpulseDE2_rdoc_insula_arousal_rdoc_samples <- runImpulseDE2(
  matCountData    = counts, 
  dfAnnotation    = metadata_rdoc_combined_insula2,
  boolCaseCtrl    = FALSE,
  vecDispersionsExternal = dispersions_deseq,
  vecConfounders  = c("Batch","Sex"),
  boolIdentifyTransients = TRUE)

# dimensional regression over social domain
metadata_rdoc_combined_insula['Time'] <- metadata_rdoc_combined_insula['social2']
metadata_rdoc_combined_insula$Time <- as.numeric(metadata_rdoc_combined_insula$Time)
metadata_rdoc_combined_insula2 <- metadata_rdoc_combined_insula[
  ,c('Sample','Condition','Time','Batch','Sex')]

objectImpulseDE2_rdoc_insula_social_rdoc_samples <- runImpulseDE2(
  matCountData    = counts, 
  dfAnnotation    = metadata_rdoc_combined_insula2,
  boolCaseCtrl    = FALSE,
  vecDispersionsExternal = dispersions_deseq,
  vecConfounders  = c("Batch","Sex"),
  boolIdentifyTransients = TRUE)

sig_res = objectImpulseDE2_rdoc_insula_cognition_rdoc_samples@dfImpulseDE2Results
