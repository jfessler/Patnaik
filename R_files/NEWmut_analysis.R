mut_table <- read.table("ptenPI3K_04thresholdCNV.csv", header = TRUE, sep = ",", 
                        stringsAsFactors = FALSE, row.names = "COMMON")
 
id_tab <- read.table("ID_key.csv", header = TRUE, stringsAsFactors = FALSE, 
                     sep = ",", row.names = "SRR_ID")

t_clusters <- t(myclusters)

merged_tab <- merge(id_tab, t_clusters, by = 0)

rownames(merged_tab) <- merged_tab$sample_ID

merged_tab$sample_ID <- NULL
merged_tab$Row.names <- NULL

merged_tab <- cbind(merged_tab, rowSums(merged_tab))
colnames(merged_tab)[ncol(merged_tab)] <- "heatScore"

merged_tab <- merge(merged_tab, mut_table, by = 0)

p1 <-ggplot(merged_tab, aes(factor(PTEN), heatScore))
p1 <- p1 + geom_boxplot() + geom_jitter() + ggtitle("PTEN, biallelic loss (< 0.4)")

p2 <-ggplot(merged_tab, aes(factor(PIK3CA), heatScore))
p2 <- p2 + geom_boxplot() + geom_jitter() + ggtitle("PIK3CA, biallelic loss (< 0.4)")

p3 <-ggplot(merged_tab, aes(factor(PIK3CB), heatScore))
p3 <- p3 + geom_boxplot() + geom_jitter() + ggtitle("PIK3CB, biallelic loss (< 0.4)")

p4 <-ggplot(merged_tab, aes(factor(PIK3CG.y), heatScore))
p4 <- p4 + geom_boxplot() + geom_jitter() + ggtitle("PIK3CG, biallelic loss (< 0.4)")

p5 <-ggplot(merged_tab, aes(factor(PIK3CD), heatScore))
p5 <- p5 + geom_boxplot() + geom_jitter() + ggtitle("PIK3CD, biallelic loss (< 0.4)")



wilcox.test(merged_tab$heatScore[merged_tab$PTEN == -2], merged_tab$heatScore[merged_tab$PTEN == 0])
#p-value = 0.3836