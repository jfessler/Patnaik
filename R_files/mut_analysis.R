mut_table <- read.table("ptenPIK3_mut.csv", header = TRUE, sep = ",", 
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
p1 <- p1 + geom_boxplot() + geom_jitter() + ggtitle("PTEN mut")

p2 <-ggplot(merged_tab, aes(factor(PIK3CA), heatScore))
p2 <- p2 + geom_boxplot() + geom_jitter() + ggtitle("PIK3CA mut")

p3 <-ggplot(merged_tab, aes(factor(PIK3CB), heatScore))
p3 <- p3 + geom_boxplot() + geom_jitter() + ggtitle("PIK3CB mut")

p4 <-ggplot(merged_tab, aes(factor(PIK3CG), heatScore))
p4 <- p4 + geom_boxplot() + geom_jitter() + ggtitle("PIK3CG mut")

p5 <-ggplot(merged_tab, aes(factor(PIK3CD.y), heatScore))
p5 <- p5 + geom_boxplot() + geom_jitter() + ggtitle("PIK3CD mut")
