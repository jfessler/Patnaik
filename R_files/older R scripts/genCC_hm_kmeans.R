library(randomcoloR)
n = 20
clusterMap <- new_d$fit.cluster
names(clusterMap) <- rownames(new_d)

rand_colors <- distinctColorPalette(n)
more_colors <- list()
for (color in rand_colors){
  more_colors <- c(more_colors, color)
}

secondColor.map <- function(CATEGORY) {more_colors[[CATEGORY]]}
allColors <- unlist(lapply(clusterMap, secondColor.map))
rowcolors <- allColors[rownames(new_cluster)]


cluster_assignment <- character()
for (g in sig_gene_list){
  cluster_assignment <- c(cluster_assignment, paste0(g,": ",new_d[g,"fit.cluster"]))
  print(paste0(g,": ",new_d[g,"fit.cluster"]))
}


for (num in c(2:20)){
  
  numClusters <- num
  
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
  
  filename = paste0("metsK12-k",numClusters,"_",format(Sys.time(), "%m-%e-%y_%H%M"),".jpg")
  
  jpeg(filename, width = 750, height = 1200)
  heatmap.2(as.matrix(new_cluster),
            Rowv = TRUE, 
            Colv = NA, 
            dendrogram = "row", 
            trace = "none",
            col = rev(redblue(100)),
            main = paste0("Mets,kmean:k12, ConsensusClusters k=",numClusters),
            ColSideColors = sidebarcolors,
            RowSideColors = rowcolors,
            breaks = seq(-5, 5, length.out = 101)
  )
  
  legend("topright",      
         legend = match(unique(rowcolors),more_colors),
         col = unique(rowcolors), 
         lty= 1,             
         lwd = 5,           
         cex=.7
  )
  
  legend("top",
         legend = cluster_assignment)
  
  dev.off()
  
}
########################################
#######################################

beep()
