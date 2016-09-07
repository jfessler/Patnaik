# Determine number of clusters
wss <- (nrow(d)-1)*sum(apply(d,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(d, iter.max = 30,
                                     centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


# K-Means Cluster Analysis
fit <- kmeans(d, 12, iter.max = 30) # 20 cluster solution
# get cluster means 
#aggregate(d,by=list(fit$cluster),FUN=mean)
# append cluster assignment
new_d <- data.frame(d, fit$cluster)

x <- character()
for (g in sig_gene_list){
  x <- c(x, paste0(g,": ",new_d[g,"fit.cluster"]))
  print(paste0(g,": ",new_d[g,"fit.cluster"]))
}

sig_gene_list <- c('CD8A', 'CCL2', 'CCL3', 'CCL4', 'CXCL9', 
                   'CXCL10', 'ICOS', 'GZMK', 'IRF1', 
                   'HLA-DMA', 'HLA-DOA', 'HLA-DOB', "HLA-DMB")

#For testing with small_mets
#sig_gene_list <- c("A1BG", "BOP1", "A1CF", "A2M", "BMX")

tcell_clusters <- list()
for (gene in sig_gene_list){
  if (gene %in% rownames(new_d)){
    tcell_clusters <- c(tcell_clusters, new_d[gene,"fit.cluster"])
  }
  else (print(paste0(gene," filtered out")))
}  

uniq_clusters <- unique(tcell_clusters)

myclusters <- data.frame()
for (cluster in uniq_clusters){
  myclusters <- rbind(myclusters, new_d[new_d$fit.cluster == cluster,])
}

myclusters <- myclusters[,1:(ncol(d))]

results = ConsensusClusterPlus(as.matrix(myclusters),
                               maxK = 20, reps = 2000, 
                               pItem = 0.8, 
                               pFeature = 1, 
                               title = paste0("inclOUT-pF1-0.5_SCALED_k12_kmeans",
                                              format(Sys.time(),"%m-%e-%y_%H%M")), 
                               clusterAlg = "hc", 
                               distance = 'euclidean',
                               seed = 20,
                               plot="png")


#####

hc <- as.hclust(hm$rowDendrogram)
cluster_map <- cutree(hc, k = 20)

for (g in sig_gene_list){
  print(paste0(g,": ",cluster_map[g]))
  }

