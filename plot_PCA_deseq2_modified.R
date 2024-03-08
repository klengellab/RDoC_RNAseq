library(genefilter)
library(ggplot2)
library(ggrepel)
library("factoextra")

plotPCA.deseq2_modified <- function (object, intgroup = "condition", 
                                     ntop = 500, returnData = FALSE) 
{
  # 求每个基因的样本间方差，并按方差降序排序，选取前5000个样本间方差最大的那些基因行
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  # 取方差最大的前5000个基因，转置后，行为观测样本，列为各种参数如基因A，做主成分分析
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  
  ## Select the PCAs and percentVar that you like instead of 1 and 2
  n_pc = ncol(pca[["x"]])
  d0 <- data.frame(pca$x)
  d <- cbind(d0,group,intgroup.df)
  # d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  #                 PC4 = pca$x[, 4], PC5 = pca$x[, 5], PC6 = pca$x[, 6],
  #                 PC7 = pca$x[, 7], PC8 = pca$x[, 8], PC9 = pca$x[, 9],
  #                 PC10 = pca$x[, 10], PC11 = pca$x[, 11], 
  #                 PC12 = pca$x[, 12],
  #                 PC13 = pca$x[, 13], PC14 = pca$x[, 14], 
  #                 PC15 = pca$x[, 15],
  #                 group = group, 
  #                 intgroup.df, name=colnames(object))
  
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:n_pc]
    return(list(d,pca))
  }
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group"))
  + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed()
  
}