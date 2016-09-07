library(ConsensusClusterPlus)
library(gplots)
library(LPE)
library(plyr)
library(dplyr)
library(RColorBrewer)

selection_method <- "polyA" # or "hybrid"

setwd("~/Documents/fessler/prostate_bioinformatics/")
RNAseq_exp <- read.table("dv_prostate_RNAseq_exp.csv", sep = ',', header = TRUE)
mets_exp <- RNAseq_exp[,32:129] #mets

mets_t <- as.data.frame(t(mets_exp))
colnames(mets_t) <- RNAseq_exp$gene_short_name
rownames(mets_t) <- colnames(mets_exp)

libSelect <- read.table(file = "librarySelection.csv", header = TRUE, 
                        stringsAsFactors = FALSE, sep = ",", row.names = "Run_s")

hybrid_mets <- libSelect[libSelect$LibrarySelection_s == "Hybrid Selection",]
polyA_mets <-  libSelect[libSelect$LibrarySelection_s == "PolyA",]

if (selection_method == "polyA"){
  mets_t <- merge(polyA_mets, mets_t, by = 0)
}
if (selection_method == "hybrid"){
  mets_t <- merge(hybrid_mets, mets_t, by = 0)
}

rownames(mets_t) <- mets_t$Row.names
mets_t <- mets_t[4:ncol(mets_t)]

## Genes expressed in less than 80% of samples removed
mets_new <- mets_t[, sapply(mets_t, function(x) { sum(x == 0) <  (0.5*(nrow(mets_t)))} )]

## Genes are rows, samples are columns
mets_new <- t(mets_new)

# Upper-quartile normalization was based on the constitutive gene counts, 
# described above, but excluding any gene that had zero counts for all of the lanes.
# Then for each lane, di was the upper-quartile (75 percentile) of all the gene counts in lane i

# upper quartile normalize
mets_uqn <- quartile.normalize(mets_new, percent = 75)

# Log2 of gene expression
mets_log2 <- log2(mets_uqn+1)
mets_clean <- mets_log2

## Remove outlier patient sample
#mets_clean$SRR3018194 <- NULL

# Clean up
rm(mets_exp, mets_t, mets_new, mets_uqn, mets_log2)
rm(RNAseq_exp)

# scale rows and columns
d <- scale(mets_clean, scale = FALSE) # scale columns
d <- t(scale(t(d), scale = FALSE)) # scale rows

jpeg(filename = "scaleF_polyA_mets_823.jpg")
## Generate heatmap
hm <- heatmap.2(d, Rowv = TRUE, Colv = TRUE,
                trace = "none",
                col = bluered(100),
                main = paste0(selection_method," mets"),
                breaks = seq(-4,4,length.out = 101),
                labRow = ""
)
dev.off()

hc <- as.hclust(hm$rowDendrogram)
