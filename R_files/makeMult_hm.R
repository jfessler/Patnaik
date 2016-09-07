library(randomcoloR)

Kmeans.cut <- function(K){
  # K-Means Cluster Analysis
  set.seed(42)
  fit <- kmeans(d, K, iter.max = 30) # 20 cluster solution
  # get cluster means 
  #aggregate(d,by=list(fit$cluster),FUN=mean)
  # append cluster assignment
  new_d <- data.frame(d, fit$cluster)
  return(new_d)
}
 

tcellCusters.pick <- function(k ,new_d){
  sig_gene_list <<- c('CD8A', 'CCL2', 'CCL3', 'CCL4', 'CXCL9', 
                     'CXCL10', 'ICOS', 'GZMK', 'IRF1', 
                     'HLA-DMA', 'HLA-DOA', 'HLA-DOB', "HLA-DMB")
  x <- character()
  print(paste0("k = ",k))
  for (g in sig_gene_list){
    x <- c(x, paste0(g,": ",new_d[g,"fit.cluster"]))
    print(paste0(g,": ",new_d[g,"fit.cluster"]))
  }
  
  tcell_clusters <- list()
  for (gene in sig_gene_list){
    if (gene %in% rownames(new_d)){
      tcell_clusters <- c(tcell_clusters, new_d[gene,"fit.cluster"])
    }
  }  
  
  uniq_clusters <- unique(tcell_clusters)
  
  myclusters <- data.frame()
  for (cluster in uniq_clusters){
    myclusters <- rbind(myclusters, new_d[new_d$fit.cluster == cluster,])
  }
  
  myclusters <- myclusters[,1:(ncol(d))]
  return(myclusters)
  
}

genCC <- function(k,myclusters){
  folder <<- paste0("k",k,"_",
                   format(Sys.time(),"%m-%e-%y_%H%M"))
  results = ConsensusClusterPlus(as.matrix(myclusters),
                                 maxK = 20, reps = 2000, 
                                 pItem = 0.8, 
                                 pFeature = 1, 
                                 title = folder, 
                                 clusterAlg = "hc", 
                                 distance = 'euclidean',
                                 seed=1262118388.71279,
                                 plot="png")
  return(results)
}

plotCC <- function(new_d, k_init, results, myclusters){
  n = k_init
  
  clusterMap <- new_d$fit.cluster
  names(clusterMap) <- rownames(new_d)
  
  rand_colors <- distinctColorPalette(n)
  more_colors <- list()
  for (color in rand_colors){
    more_colors <- c(more_colors, color)
  }
  
  cluster_assignment <- character()
  for (g in sig_gene_list){
    cluster_assignment <- c(cluster_assignment, paste0(g,": ",new_d[g,"fit.cluster"]))
    #print(paste0(g,": ",new_d[g,"fit.cluster"]))
  }
  
  
  for (num in seq(8,14,2)){
  
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
    
    my_colors <- colorspace::diverge_hsv(numClusters)
    # Add side bar colors
    color.map<-function(CATEGORY) {my_colors[CATEGORY]}
    sidebarcolors <- unlist(lapply(c_cluster$consensusClass, color.map))
    
    secondColor.map <- function(CATEGORY) {more_colors[[CATEGORY]]}
    allColors <- unlist(lapply(clusterMap, secondColor.map))
    rowcolors <- allColors[rownames(new_cluster)]
    
    filename = paste0(folder,"/metsK",k_init,"_k",numClusters,"_",format(Sys.time(), "%m-%e-%y_%H%M"),".jpg")
    
    jpeg(filename, width = 750, height = 1200)
    heatmap.2(as.matrix(new_cluster),
              Rowv = TRUE, 
              Colv = NA, 
              dendrogram = "row", 
              trace = "none",
              col = rev(redblue(100)),
              main = paste0("Mets,kmean:,",k_init,"ConsensusClusters k=",numClusters),
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
  
}

for (k_num in seq(25,50,25)){
  NEW_D <- Kmeans.cut(K = k_num)
  MY_CLUS <- tcellCusters.pick(k = k_num, new_d = NEW_D)
  print(paste0("Number of genes: ",dim(MY_CLUS)[1]))
  RESULTS <- genCC(k = k_num, myclusters = MY_CLUS)
  plotCC(new_d = NEW_D, k_init = k_num, results = RESULTS, myclusters = MY_CLUS)
}

beep()