
library(tidyr)

setwd("~/Documents/fessler/prostate_bioinformatics/")
TCGA_data <- read.table("../TCGA_data/gdac.broadinstitute.org_PRAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/PRAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
                        sep = '\t', header = TRUE, stringsAsFactors = FALSE)
colnames(TCGA_data)[1] <- "Gene_ID"
TCGA_data <- TCGA_data[2:nrow(TCGA_data),]

TCGA_data <- separate(data = TCGA_data, col = Gene_ID, into = c("gene_id","gene_num"), sep = "\\|")

library(ConsensusClusterPlus)
library(gplots)
library(LPE)
library(dplyr)
library(RColorBrewer)
set.seed(20)

rownames(TCGA_data) = make.names(TCGA_data$gene_id, unique=TRUE) # rename rows gene_id
TCGA_data <- TCGA_data[,3:ncol(TCGA_data)] # remove gene info
TCGA_data2 <- sapply(TCGA_data, as.numeric)
rownames(TCGA_data2) <- rownames(TCGA_data)
TCGA_data <- TCGA_data2
rm(TCGA_data2)
t_data <- t(TCGA_data)
t_data <- as.data.frame(t_data)

## Genes expressed in less than 80% of samples removed
t_new <- t_data[, sapply(t_data, function(x) { sum(x == 0) <  (0.8*(nrow(t_data)))} )]

## Genes are rows, samples are columns
t_new <- t(t_new)

# Upper-quartile normalization was based on the constitutive gene counts, 
# described above, but excluding any gene that had zero counts for all of the lanes.
# Then for each lane, di was the upper-quartile (75 percentile) of all the gene counts in lane i

# upper quartile normalize
t_uqn <- quartile.normalize(t_new, percent = 75)

# Log2 of gene expression
t_log2 <- log2(t_uqn+1)
tcga_clean <- t_log2

# Clean up
rm(t_data, t_log2, t_new, t_uqn)
rm(TCGA_data)

# scale rows and columns
d <- scale(tcga_clean) # scale columns
d <- t(scale(t(d))) # scale rows

file_name <- paste0("TCGA_prostatePrimary_",format(Sys.time(), "%m-%e-%y_%H%M"),".tiff")
tiff(file_name, width = 6, height = 8, units = 'in', res = 300, compression = 'none')

## Generate heatmap
hm <- heatmap.2(d, Rowv = TRUE, Colv = TRUE,
                trace = "none",
                col = rev(redblue(100)),
                #scale = "row",
                cexCol = 0.2,
                cexRow = 0.2,
                main = "TCGA primary",
                breaks = seq(-5,5,length.out = 101)
)

dev.off()

hc <- as.hclust(hm$rowDendrogram)

kay = 5000
cluster_map <- cutree(hc, k = kay)

## T-Cell Signature Transcipts
#CD8A, CCL2, CCL3, CCL4, CXCL9, CXCL10, ICOS, 
#GZMK, IRF1, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB

sig_gene_list <- c('CD8A', 'CCL2', 'CCL3', 'CCL4', 'CXCL9', 
                   'CXCL10', 'ICOS', 'GZMK', 'IRF1', 
                   'HLA.DMA', 'HLA.DOA', 'HLA.DOB', "HLA.DMB")

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

more_colors <- colorspace::rainbow_hcl(kay)
secondColor.map <- function(CATEGORY) {more_colors[CATEGORY]}
allColors <- unlist(lapply(cluster_map, secondColor.map))
rowcolors <- allColors[rownames(myclusters)]

file_name <- paste0("TCGA_prostatePrimary_k=",kay,format(Sys.time(), "%m-%e-%y_%H%M"),".tiff")
tiff(file_name, width = 6, height = 8, units = 'in', res = 300, compression = 'none')

## Generate heatmap
hm <- heatmap.2(as.matrix(myclusters), 
                Rowv = TRUE, Colv = TRUE,
                trace = "none",
                col = rev(redblue(100)),
                labCol = "",
                cexRow = 0.2,
                main = paste0("TCGA primary, k=",kay),
                RowSideColors = rowcolors,
                breaks = seq(-5,5,length.out = 101)
)

dev.off()

results = ConsensusClusterPlus(as.matrix(myclusters),
                               maxK = 20, reps = 2000, 
                               pItem = 0.8, 
                               pFeature = 1, 
                               title = "TCGA_primary_k5000", 
                               clusterAlg = "hc", 
                               distance = 'euclidean',
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

###########
file_name <- paste0("TCGA_h25_tcellGenes",format(Sys.time(), "%m-%e-%y_%H%M"),".tiff")
tiff(file_name, width = 6, height = 8, units = 'in', res = 300, compression = 'none')

## Generate heatmap
hm <- heatmap.2(as.matrix(myclusters), Rowv = TRUE, Colv = TRUE,
                trace = "none",
                col = rev(redblue(100)),
                #scale = "row",
                labCol = "",
                cexRow = 0.15,
                main = "TCGA prim, h=25 Tcell Sig",
                breaks = seq(-5,5,length.out = 101)
)

dev.off()

