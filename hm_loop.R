for (klust in seq(5, 110, 5)){

cluster_map <- cutree(hc, k = klust)

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

more_colors <- colorspace::rainbow_hcl(1000)
secondColor.map <- function(CATEGORY) {more_colors[CATEGORY]}
allColors <- unlist(lapply(cluster_map, secondColor.map))
rowcolors <- allColors[rownames(myclusters)]

filename = paste0("Mets_k",klust,"_",format(Sys.time(), "%m-%e-%y_%H%M"),".jpg")

jpeg(filename, width = 750, height = 1200)

heatmap.2(as.matrix(myclusters),
          Rowv = TRUE, 
          Colv = TRUE, 
          dendrogram = "row", 
          trace = "none",
          col = rev(redblue(100)),
          main = paste0("Mets unsupervised clustering, k=",klust),
          #ColSideColors = sidebarcolors,
          RowSideColors = rowcolors,
          breaks = seq(-3, 3, length.out = 101)
)
dev.off()
print(klust)
}


#####

