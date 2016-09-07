
RNAseq_exp <- read.table("dv_prostate_RNAseq_exp.csv", sep = ',', header = TRUE)

norm_exp <- RNAseq_exp[,1:6]
norm_t <- as.data.frame(t(norm_exp[,4:ncol(norm_exp)]))
colnames(norm_t) <- norm_exp$gene_short_name

primary_exp <- cbind(RNAseq_exp[,1:3],RNAseq_exp[,7:31])
primary_t <- as.data.frame(t(primary_exp[,4:ncol(primary_exp)]))
colnames(primary_t) <- primary_exp$gene_short_name

mets_exp <- cbind(RNAseq_exp[,1:3],RNAseq_exp[,32:129])
mets_t <- as.data.frame(t(mets_exp[,4:ncol(mets_exp)]))
colnames(mets_t) <- mets_exp$gene_short_name

gene_list <- norm_exp[,"gene_short_name"]

wt_table <- data.frame()

for (gene in gene_list){
  wt_prim_pvalue <- NA
  if (sd(primary_t[[gene]]) != 0 & sd(norm_t[[gene]]) != 0){
    wt_prim <- wilcox.test(norm_t[[gene]], primary_t[[gene]], paired=FALSE)
    wt_prim_pvalue <- wt_prim$p.value    
  }

  wt_mets_pvalue <- NA
  if (sd(mets_t[[gene]]) != 0 & sd(norm_t[[gene]]) != 0){
    wt_mets <- wilcox.test(norm_t[[gene]], mets_t[[gene]], paired=FALSE)
    wt_mets_pvalue <- wt_mets$p.value
  }
  wt_table <- rbind(wt_table, data.frame(gene, wt_prim_pvalue, wt_mets_pvalue))
}

# Account for multiple hypothesis testing
prim_adj_pvalue <- p.adjust(wt_table$wt_prim_pvalue, method = "BH")
mets_adj_pvalue <- p.adjust(wt_table$wt_mets_pvalue, method = "BH")

wt_table <- cbind(wt_table, prim_adj_pvalue, mets_adj_pvalue)



met_lt_05_wt <- filter(wt_table, mets_adj_pvalue < 0.05)

