library(randomcoloR)

#2nd heatmap with samples ordered by Consensus Clustering

numClusters <- 13
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

my_colors <- colorspace::diverge_hsv(numClusters)
# Add side bar colors
color.map<-function(CATEGORY) {my_colors[CATEGORY]}
sidebarcolors <- unlist(lapply(c_cluster$consensusClass, color.map))

secondColor.map <- function(CATEGORY) {more_colors[[CATEGORY]]}
allColors <- unlist(lapply(clusterMap, secondColor.map))
rowcolors <- allColors[rownames(new_cluster)]

my.breaks <- c(seq(-4.4, -1, length.out=60),seq(-1.000001, 2, length.out=30),seq(2.00001,4.4, length.out=60))


filename = paste0(folder,"/metsK60",format(Sys.time(), "%m-%e-%y_%H%M"),".jpg")
jpeg(filename, width = 750, height = 1200)

heatmap.2(new_cluster,
          Rowv = TRUE, 
          Colv = NA, 
          dendrogram = "row", 
          trace = "none",
          col=bluered(149),
          main = paste0("Mets, ConsensusClusters k=",numClusters),
          ColSideColors = sidebarcolors,
          RowSideColors = rowcolors,
          breaks = my.breaks
)

