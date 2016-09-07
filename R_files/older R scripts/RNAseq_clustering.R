library(ConsensusClusterPlus)
library(gplots)
library(LPE)
library(dplyr)
set.seed(20)

setwd("~/Documents/fessler")

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

# upper quartile normalize, 2 = columns
mets_uqn <- quartile.normalize(mets_new,percent = 75)

# Log2 of gene expression
mets_log2 <- log2(mets_uqn+1)
mets_clean <- mets_log2

# Clean up
rm(mets_exp, mets_t, mets_new, mets_uqn, mets_log2)
rm(RNAseq_exp)

d <- scale(mets_clean)

## Generate heatmap
hm <- heatmap.2(d, Rowv = TRUE, Colv = TRUE, dendrogram = "none",
                trace = "none",
                col = rev(redblue(50)),
                main = "all mets")

hc <- as.hclust(hm$rowDendrogram)

cluster_map <- cutree(hc, k = 12)

## T-Cell Signature Transcipts
#CD8A, CCL2, CCL3, CCL4, CXCL9, CXCL10, ICOS, 
#GZMK, IRF1, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB

sig_gene_list <- c('CD8A', 'CCL2', 'CCL3', 'CCL4', 'CXCL9', 
                   'CXCL10', 'ICOS', 'GZMK', 'IRF1', 
                   'HLA-DMA', 'HLA-DOA', 'HLA-DOB', "HLA-DMB")

#For testing with small_mets
#sig_gene_list <- c("A1BG", "BOP1", "A1CF", "A2M", "BMX")

tcell_clusters <- list()
for (gene in sig_gene_list){
  tcell_clusters <- c(tcell_clusters, cluster_map[gene])
}

uniq_clusters <- unique(tcell_clusters)

myclusters <- data.frame()
for (cluster in uniq_clusters){
  myclusters <- rbind(myclusters, d[cluster_map == cluster,])
}

results = ConsensusClusterPlus(as.matrix(myclusters),
                               maxK = 20, reps = 2000, 
                               pItem = 0.8, 
                               pFeature = 1, 
                               title = "cluster_tcells_k12", 
                               clusterAlg = "hc", 
                               distance = 'euclidean',
                               seed = 20,
                               plot="png")

#2nd heatmap with samples ordered by Consensus Clustering

# Cluster different samples fall into
cOrder <- as.data.frame(results[[15]]["consensusClass"])
t_myclusters <- as.data.frame(t(myclusters))

#merge by row names
c_cluster <- merge(cOrder, t_myclusters, by = 0)
rownames(c_cluster) <- c_cluster$Row.names #add sampleID as row name

#sort by cluster and transpose
c_cluster <-arrange(c_cluster, consensusClass)
rownames(c_cluster) <- c_cluster$Row.names
c_cluster <- t(c_cluster[,3:ncol(c_cluster)]) #remove class and sampleID

#remove consensus class
#cOrder <- cbind(sampleID <- row.names(cOrder), cOrder)
#arrange(cOrder, consensusClass)

#dendrogram 
#cTree <- results[[2]]["consensusTree"]

heatmap.2(scale(c_cluster, center = -1, scale = FALSE),
          Rowv = TRUE, 
          Colv = FALSE, dendrogram = "row", 
          trace = "none", 
          col = rev(redblue(50)), 
          main = "mets, CC hAve")

########################################
#######################################


colors = c(seq(-2,-0.75,length=100),
           seq(-0.75,-1,length=100),
           seq(-1,1.5,length=100))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

heatmap.2(c_cluster,
         Rowv = TRUE, 
         Colv = FALSE, dendrogram = "row", 
         trace = "none", 
         col = my_palette, 
         main = "mets, CC hAve")

updated_cluster <- mets_clean[rownames(myclusters),]

heatmap.2(as.matrix(updated_cluster),
          Rowv = TRUE, 
          Colv = FALSE, dendrogram = "row", 
          trace = "none", 
          col = rev(redblue(50)), 
          main = "mets, CC")


collist <- c("D1","D2","D3","D4")
sel <- apply(data[,collist],1,function(row) "E" %in% row)
data[sel,]

