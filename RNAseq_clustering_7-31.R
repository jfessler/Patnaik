library(ConsensusClusterPlus)
library(gplots)
library(LPE)
library(dplyr)
library(RColorBrewer)
#set.seed(20)

setwd("~/Documents/fessler/prostate_bioinformatics/")

RNAseq_exp <- read.table("dv_prostate_RNAseq_exp.csv", sep = ',', header = TRUE)

mets_exp <- RNAseq_exp[,32:129] #mets

#primary_exp <- RNAseq_exp[,7:31] #primary
#mets_exp <- primary_exp

mets_t <- as.data.frame(t(mets_exp))

colnames(mets_t) <- RNAseq_exp$gene_short_name
rownames(mets_t) <- colnames(mets_exp)

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
mets_clean$SRR3018194 <- NULL

# Clean up
rm(mets_exp, mets_t, mets_new, mets_uqn, mets_log2)
rm(RNAseq_exp)

# scale rows and columns
d <- scale(mets_clean) # scale columns
d <- t(scale(t(d))) # scale rows

## Generate heatmap
hm <- heatmap.2(d, Rowv = TRUE, Colv = TRUE,
                trace = "none",
                col = rev(redblue(101)),
                #scale = "row",
                main = "all mets, rows+cols scaled"
                )

hc <- as.hclust(hm$rowDendrogram)
cluster_map <- cutree(hc, k = 60)

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
  if (gene %in% names(cluster_map)){
    tcell_clusters <- c(tcell_clusters, cluster_map[gene])
  }
  else (print(paste0(gene," filtered out")))
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
                               title = paste0("rmOUT_pearson_SCALED_k60_ctree_",
                                              format(Sys.time(),"%m-%e-%y_%H%M")), 
                               clusterAlg = "hc", 
                               distance = 'pearson',
                               seed = 20,
                               plot="png")

#2nd heatmap with samples ordered by Consensus Clustering

numClusters <- 12

# Cluster different samples fall into
cOrder <- as.data.frame(results[[numClusters]]["consensusClass"])
updated_cluster <- d[rownames(myclusters),]
t_myclusters <- as.data.frame(t(updated_cluster))

#merge by row names
c_cluster <- merge(cOrder, t_myclusters, by = 0)
rownames(c_cluster) <- c_cluster$Row.names #add sampleID as row name

#sort by cluster and transpose
c_cluster <-arrange(c_cluster, consensusClass)
rownames(c_cluster) <- c_cluster$Row.names
new_cluster <- t(c_cluster[,3:ncol(c_cluster)]) #remove class and sampleID

#remove consensus class
#cOrder <- cbind(sampleID <- row.names(cOrder), cOrder)
#arrange(cOrder, consensusClass)

my_colors <- colorspace::diverge_hsv(numClusters)
# Add side bar colors
color.map<-function(CATEGORY) {my_colors[CATEGORY]}
sidebarcolors <- unlist(lapply(c_cluster$consensusClass, color.map))

more_colors <- colorspace::rainbow_hcl(50)
secondColor.map <- function(CATEGORY) {more_colors[CATEGORY]}
allColors <- unlist(lapply(cluster_map, secondColor.map))
rowcolors <- allColors[rownames(new_cluster)]

heatmap.2(new_cluster,
          Rowv = TRUE, 
          Colv = NA, 
          dendrogram = "row", 
          trace = "none",
          col = rev(redblue(100)),
          main = paste0("Mets, ConsensusClusters k=",numClusters),
          ColSideColors = sidebarcolors,
          RowSideColors = rowcolors,
          breaks = seq(-3, 3, length.out = 101)
          )

########################################
#######################################
heatmap.2(as.matrix(edited_myclusters),
          Rowv = TRUE, 
          Colv = TRUE, 
          dendrogram = "row", 
          trace = "none",
          col = rev(redblue(100)),
          main = paste0("NOT EDITED Mets, k = 2000"),
          #ColSideColors = sidebarcolors,
          #RowSideColors = rowcolors,
          cexRow = 0.2,
          cexCol = 0.2,
          breaks = seq(-4, 4, length.out = 101)
)
