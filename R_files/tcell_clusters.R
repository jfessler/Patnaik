library(ConsensusClusterPlus)
library(gplots)
library(LPE)
library(dplyr)
library(RColorBrewer)
set.seed(20)

setwd("~/Documents/fessler/prostate_bioinformatics/")

RNAseq_exp <- read.table("dv_prostate_RNAseq_exp.csv", sep = ',', header = TRUE)

mets_exp <- RNAseq_exp[,32:129] #mets

#primary_exp <- RNAseq_exp[,7:31] #primary

mets_t <- as.data.frame(t(mets_exp))

colnames(mets_t) <- RNAseq_exp$gene_short_name
rownames(mets_t) <- colnames(mets_exp)

## Genes expressed in less than 80% of samples removed
mets_new <- mets_t[, sapply(mets_t, function(x) { sum(x == 0) <  (0.8*(nrow(mets_t)))} )]

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

# Clean up
rm(mets_exp, mets_t, mets_new, mets_uqn, mets_log2)
rm(RNAseq_exp)

d <- scale(mets_clean)
sig_gene_list <- c('CD8A', 'CCL2', 'CCL3', 'CCL4', 'CXCL9', 
                   'CXCL10', 'ICOS', 'GZMK', 'IRF1',
                   'HLA-DMA', 'HLA-DOA', 'HLA-DOB', "HLA-DMB")

tcell_clust <- d[sig_gene_list,]
my.breaks <- c(seq(-4, -.51, length.out=60),seq(-.50, .5, length.out=30),seq(.51,4, length.out=60))

## Generate heatmap
hm <- heatmap.2(tcell_clust, 
                Rowv = TRUE, Colv = TRUE,
                trace = "none",
                col = bluered(149),
                main = "polyA mets, Tcell genes",
                breaks = my.breaks,
                density.info = "none"
                )

## T-Cell Signature Transcipts
#CD8A, CCL2, CCL3, CCL4, CXCL9, CXCL10, ICOS, 
#GZMK, IRF1, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB

tcell_clust <- d[sig_gene_list,]

results = ConsensusClusterPlus(tcell_clust,
                               maxK = 20, reps = 2000, 
                               pItem = 0.8, 
                               pFeature = 1, 
                               title = "primaryTcells_ONLY", 
                               clusterAlg = "hc", 
                               distance = 'euclidean',
                               seed = 20,
                               plot="png")

numClusters <- 8

# Cluster different samples fall into
cOrder <- as.data.frame(results[[numClusters]]["consensusClass"])
t_myclusters <- as.data.frame(t(tcell_clust))

#merge by row names
c_cluster <- merge(cOrder, t_myclusters, by = 0)
rownames(c_cluster) <- c_cluster$Row.names #add sampleID as row name

#sort by cluster and transpose
c_cluster <-arrange(c_cluster, consensusClass)
rownames(c_cluster) <- c_cluster$Row.names
new_cluster <- t(c_cluster[,3:ncol(c_cluster)]) #remove class and sampleID

## Better colors
my_colors <- colorspace::diverge_hsv(numClusters)
sidebarcolors <- unlist(lapply(c_cluster$consensusClass, color.map))

heatmap.2(new_cluster,
          Rowv = TRUE, 
          Colv = NA, dendrogram = "row", 
          trace = "none",
          col = rev(redblue(100)),
          main = "primary CC, Tcells",
          ColSideColors = sidebarcolors)




# Add side bar colors
color.map<-function(CATEGORY) {if(CATEGORY=="1") "blue" 
  else if(CATEGORY=="2") "green" 
  else if(CATEGORY=="3") "yellow"
  else if(CATEGORY=="4") "orange"
  else if(CATEGORY=="5") "red"
  else if(CATEGORY=="6") "purple"
  else if(CATEGORY=="7") "black"
  else if(CATEGORY=="8") "brown"
  else if(CATEGORY=="9") "gray"
}
sidebarcolors <- unlist(lapply(c_cluster$consensusClass, color.map))
