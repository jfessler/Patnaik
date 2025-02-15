if(exists("hc") == FALSE){
  source("R_files/makeHC.R")
  }

library(randomcoloR)
library(beepr)

####################################################################
cut_k <<- 500
cc_k <<- 13

coldClus <<- 7
hotClus <<- 2
####################################################################


cutree.cut <- function(K){
  # K-Means Cluster Analysis
  set.seed(42)
  cluster_map <- cutree(hc, k = K)
  return(cluster_map)
}


tcellCusters.pick <- function(k, cluster_map){
  sig_gene_list <<- c('CD8A', 'CCL2', 'CCL3', 'CCL4', 'CXCL9', 
                      'CXCL10', 'ICOS', 'GZMK', 'IRF1') 
#                      'HLA-DMA', 'HLA-DOA', 'HLA-DOB', "HLA-DMB")
  
  tcell_clusters <- list()
  for (gene in sig_gene_list){
    if (gene %in% names(cluster_map)){
      tcell_clusters <- c(tcell_clusters, cluster_map[gene])
    }
  }  
  
  uniq_clusters <- unique(tcell_clusters)
  
  myclusters <- data.frame()
  for (cluster in uniq_clusters){
    myclusters <- rbind(myclusters, d[cluster_map == cluster,])
  }
  
  return(myclusters)
  
}

genCC <- function(k,myclusters){
  folder <<- paste0("resultImages/k",k,"_",
                    format(Sys.time(),"%m-%e-%y_%H%M"))
  results = ConsensusClusterPlus(as.matrix(myclusters),
                                 maxK = 20, reps = 2000, 
                                 pItem = 0.8, 
                                 pFeature = 1, 
                                 title = folder, 
                                 clusterAlg = "hc", 
                                 distance = 'pearson',
                                 seed=1262118388.71279,
                                 plot="pdf",
                                 verbose = FALSE)
  return(results)
}

plotCC <- function(k_init, results, myclusters, cluster_map){
  n = k_init
  
  rand_colors <- distinctColorPalette(n)
  more_colors <- list()
  for (color in rand_colors){
    more_colors <- c(more_colors, color)
  }
  
  cluster_assignment <- character()
  for (g in sig_gene_list){
    print(paste0(g,": ",cluster_map[g]))
    cluster_assignment <- c(cluster_assignment, paste0(g,": ",cluster_map[g]))
  }
  
  for (num in c(cc_k)){
    
    numClusters <- num
    # Cluster different samples fall into
    cOrder <- as.data.frame(results[[numClusters]]["consensusClass"])
    cOrder$heatScore <- NA
    cOrder$sampleMed <- NA
    
    t_myclusters <- as.data.frame(t(myclusters))
    
    #merge by row names
    c_cluster <- merge(cOrder, t_myclusters, by = 0)
    rownames(c_cluster) <- c_cluster$Row.names #add sampleID as row name
    
    for (cclus in 1:num){
      median_sig <- median(as.matrix(c_cluster[c_cluster$consensusClass == cclus,5:ncol(c_cluster)]))
      c_cluster$heatScore[c_cluster$consensusClass == cclus] <- median_sig
    }
    
    for (SRR in rownames(c_cluster)){
      median_row <- median(as.matrix(c_cluster[SRR,5:ncol(c_cluster)]))
      c_cluster[SRR,"sampleMed"] <- median_row
    }
    
    #sort by cluster and transpose
    heatOrder <- rank(c_cluster$heatScore, ties.method = "max")
    ranks <- sort(unique(c_cluster$heatScore))
    for (pos in 1:num){
      c_cluster$heatOrder[c_cluster$heatScore == ranks[pos]] <- pos
    }
    
    c_cluster <-arrange(c_cluster, heatScore, sampleMed)
    rownames(c_cluster) <- c_cluster$Row.names
    new_cluster <- t(c_cluster[,4:ncol(c_cluster)]) #remove class and sampleID
    
    #browser()
    
    coldSRRs <<- rownames(subset(c_cluster, heatOrder <= coldClus))
    hotSRRs <<- rownames(subset(c_cluster, heatOrder > (numClusters - hotClus)))

    my_colors <- colorspace::diverge_hsv(numClusters)
    # Add side bar colors
    color.map<-function(CATEGORY) {my_colors[CATEGORY]}
    sidebarcolors <- unlist(lapply(c_cluster$heatOrder, color.map))
    
    secondColor.map <- function(CATEGORY) {more_colors[[CATEGORY]]}
    allColors <- unlist(lapply(cluster_map, secondColor.map))
    rowcolors <- allColors[rownames(new_cluster)]
    
    my.breaks <- c(seq(-4, -.31, length.out=60),seq(-.30, .3, length.out=30),seq(.31,4, length.out=60))
    
    filename = paste0(folder,"/metsK",k_init,"_k",numClusters,"_",format(Sys.time(), "%m-%e-%y_%H%M"),".tiff")
    
    tiff(filename, width = 6, height = 8, units = 'in', res = 300, compression = 'none')
    heatmap.2(as.matrix(new_cluster),
              Rowv = TRUE, 
              Colv = NA, 
              dendrogram = "none", 
              trace = "none",
              col = bluered(149),
              #main = paste0("Mets w/ outl:",k_init," CC k=",numClusters,"#: ",numGENES),
              main = "Metastatic Prostate Tumors (43)",
              ColSideColors = sidebarcolors,
              #RowSideColors = rowcolors,
              breaks = my.breaks,
              density.info = "none",
              labRow = "",
              cexCol = 0.7
              #labCol = ""
    )
    
  
    dev.off()
    
  }
  ########################################
  #######################################
  
}

for (k_num in c(cut_k)){
  print(paste0("k = ",k_num))
  CLUS_MAP <- cutree.cut(K = k_num)
  MY_CLUS <- tcellCusters.pick(k = k_num, cluster_map = CLUS_MAP)
  numGENES <<- dim(MY_CLUS)[1]
  print(paste0("Number of genes: ",numGENES))
  RESULTS <- genCC(k = k_num, myclusters = MY_CLUS)
  plotCC(k_init = k_num, results = RESULTS, myclusters = MY_CLUS, cluster_map = CLUS_MAP)
}

source('R_files/hotCold_CNVmutAnalysis_presentationVersion.R')

beep()