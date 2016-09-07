results = ConsensusClusterPlus(tcell_prim_d,
                               maxK = 20, reps = 2000, 
                               pItem = 0.8, 
                               pFeature = 1, 
                               title = "primary_k2000Mets_TcellGenes", 
                               clusterAlg = "hc", 
                               distance = 'euclidean',
                               seed = 20,
                               plot="png")

for (i in 2:20){
numClusters <- i

# Cluster different samples fall into
cOrder <- as.data.frame(results[[numClusters]]["consensusClass"])
updated_cluster <- primary_d[rownames(tcell_prim_d),]
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

filename <- paste0('primary_',format(Sys.time(),"%m-%e-%y_%H%M"),".jpg")
jpeg(filename, width = 750, height = 1200)

heatmap.2(new_cluster,
          Rowv = TRUE, 
          Colv = NA, 
          dendrogram = "row", 
          trace = "none",
          col = rev(redblue(100)),
          main = paste0("Primary, genes from k=2000 Mets, unsupervised"))
          ColSideColors = sidebarcolors)
          #breaks = seq(-3, 3, length.out = 101)

dev.off()
}


sorted_d <- tcell_prim_d[,order(colSums(tcell_prim_d))]

jpeg("sortedPrim_k2000MetsGenes.jpg", width = 750, height = 1200)
heatmap.2(sorted_d,
          Rowv = TRUE, 
          Colv = NA, 
          dendrogram = "row", 
          trace = "none",
          col = rev(redblue(100)),
          main = "Primary, T cell gene signature: Sorted")
dev.off()
