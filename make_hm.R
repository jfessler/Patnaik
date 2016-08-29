numClusters <- 16

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

more_colors <- colorspace::rainbow_hcl(2000)
secondColor.map <- function(CATEGORY) {more_colors[CATEGORY]}
allColors <- unlist(lapply(cluster_map, secondColor.map))
rowcolors <- allColors[rownames(new_cluster)]

#filename = paste0("MetsCC_K100_k",numClusters,"_",format(Sys.time(), "%m-%e-%y_%H%M"),".jpg")

#jpeg(filename, width = 750, height = 1200)
heatmap.2(new_cluster,
          Rowv = TRUE, 
          Colv = NA, 
          dendrogram = "row", 
          trace = "none",
          col = rev(redblue(100)),
          main = paste0("Mets k=2000, ConsensusClusters k=",numClusters),
          ColSideColors = sidebarcolors,
          RowSideColors = rowcolors,
          breaks = seq(-4, 4, length.out = 101),
          keysize = 1.5,
          density.info = "none"
          )
#dev.off()


tiff("hm_sorted.tiff", width = 6, height = 8, units = 'in', res = 300, compression = 'none')

heatmap.2(new_cluster,
          Rowv = TRUE, 
          Colv = NA, 
          dendrogram = "row", 
          trace = "none",
          col = rev(redblue(100)),
          #main = paste0("Mets k=2000, ConsensusClusters k=",numClusters),
          #ColSideColors = sidebarcolors,
          #RowSideColors = rowcolors,
          breaks = seq(-4, 4, length.out = 101),
          keysize = 1.5,
          density.info = "none",
          labRow = "",
          labCol = ""
          )


sorted_cluster <- myclusters[,order(colSums(myclusters))]
file_name <- paste0("TCGA_k5000_tcellSORTED",format(Sys.time(), "%m-%e-%y_%H%M"),".tiff")
tiff(file_name, width = 8, height = 8, units = 'in', res = 400, compression = 'none')

heatmap.2(as.matrix(sorted_cluster),
          Rowv = TRUE, 
          Colv = NA, 
          dendrogram = "row", 
          trace = "none",
          col = rev(redblue(100)),
          main = paste0("TCGA k=5000, sorted"), 
          breaks = seq(-5, 5, length.out = 101),
          RowSideColors = rowcolors,
          keysize = 1.5,
          density.info = "none",
          labCol = "",
          cexRow = 0.2
)

dev.off()

####'

### Compare primary to genes in mets new_cluster ###
## Create heatmap with scaled data (d) from primary samples
tcell_prim_d <- primary_d[(rownames(primary_d) %in% rownames(myclusters)),]


heatmap.2(new_d,
          Rowv = TRUE, 
          Colv = TRUE, 
          dendrogram = "row", 
          trace = "none",
          col = rev(redblue(100)),
          main = "Primary, T cell gene signature"
          #ColSideColors = sidebarcolors,
          #RowSideColors = rowcolors,
)

sorted_d <- myclusters[,order(colSums(myclusters))]

file_name <- paste0("metsk60_tcellSORTED",format(Sys.time(), "%m-%e-%y_%H%M"),".tiff")
tiff(file_name, width = 6, height = 8, units = 'in', res = 300, compression = 'none')

my.breaks <- c(seq(-4, -.41, length.out=60),seq(-.40, .4, length.out=30),seq(.41,4, length.out=60))

heatmap.2(as.matrix(sorted_d),
          Rowv = TRUE, 
          Colv = NA, 
          dendrogram = "row", 
          trace = "none",
          col = bluered(149),
          #main = "TCGA k2000 Tcell gene signature: Sorted",
          breaks = my.breaks,
          #cexRow = 0.25,
          labCol = "",
          labRow = "",
          density.info = "none"
)

dev.off()
