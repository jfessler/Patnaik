
library(preprocessCore)
library(gplots)

setwd("~/Patnaik")

RNAseq_exp <- read.table("dv_prostate_RNAseq_exp.csv", sep = ',', header = TRUE)

mets_exp <- cbind(RNAseq_exp[,1:3],RNAseq_exp[,32:129])
mets_t <- as.data.frame(t(mets_exp[,4:ncol(mets_exp)]))
colnames(mets_t) <- mets_exp$gene_short_name

gene_list <- mets_exp[,"gene_short_name"]

rm(mets_exp)
## Genes expressed in less than 80% of samples removed
mets_new <- mets_t[, sapply(mets_t, function(x) { sum(x == 0) <  (0.8*(nrow(mets_t)))} )]

rm(mets_t)
## Upper-quantile normalized 
mets_norm <- normalize.quantiles(as.matrix(mets_new))
colnames(mets_norm) <- colnames(mets_new)
rownames(mets_norm) <- rownames(mets_new)

rm(mets_new)

## Make heatmap
heatmap.2(mets_norm)


### oTHERS 
  norm_exp <- RNAseq_exp[,1:6]
  norm_t <- as.data.frame(t(norm_exp[,4:ncol(norm_exp)]))
  colnames(norm_t) <- norm_exp$gene_short_name
  
  primary_exp <- cbind(RNAseq_exp[,1:3],RNAseq_exp[,7:31])
  primary_t <- as.data.frame(t(primary_exp[,4:ncol(primary_exp)]))
  colnames(primary_t) <- primary_exp$gene_short_name
  