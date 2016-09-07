log2.cnv <- read.table(file = "ptenPIK3C_log2cnv.csv", sep = ",", 
                       header = TRUE, stringsAsFactors = FALSE, 
                       row.names = "COMMON")

id_tab <- read.table("ID_key.csv", header = TRUE, stringsAsFactors = FALSE, 
                     sep = ",", row.names = "SRR_ID")

t_clusters <- t(myclusters)
merged_tab <- merge(id_tab, t_clusters, by = 0)
rownames(merged_tab) <- merged_tab$sample_ID

merged_tab$sample_ID <- NULL
merged_tab$Row.names <- NULL

merged_tab <- cbind(merged_tab, rowSums(merged_tab))
colnames(merged_tab)[ncol(merged_tab)] <- "heatScore"

merged_tab <- merge(merged_tab, log2.cnv, by = 0)

plot(merged_tab$heatScore, merged_tab$PTEN)
