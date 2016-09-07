sum_clusters <- rbind(colSums(myclusters),myclusters)

rownames(sum_clusters)[1] <- "heatScore"

sum_clusters <- as.data.frame(t(sum_clusters))

allData <- as.data.frame(t(d))
genes <- colnames(allData)

corTab <- data.frame()

for (gene in genes){
  PearsonCorrelation <- cor(sum_clusters$heatScore, allData[[gene]])
  corTab <- rbind(corTab, data.frame(gene, PearsonCorrelation))
}