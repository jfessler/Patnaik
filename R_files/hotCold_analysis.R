library(ggplot2)

hotSamples <- read.table("HOTtumors.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE, row.names = "x")
coldSamples <- read.table("COLDtumors.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE, row.names = "x")

id_tab <- read.table("ID_key.csv", header = TRUE, stringsAsFactors = FALSE, 
                     sep = ",", row.names = "SRR_ID")

hotSamples <- merge(hotSamples, id_tab, by = 0)
row.names(hotSamples) <- hotSamples$sample_ID
hotSamples$Row.names <- NULL
numHot <- dim(hotSamples)[1]

coldSamples <- merge(coldSamples, id_tab, by= 0)
row.names(coldSamples) <- coldSamples$sample_ID
coldSamples$Row.names <- NULL
numCold <- dim(coldSamples)[1]

cnv04_tab <- read.table("ptenPI3K_04thresholdCNV.csv", sep = ",", header = TRUE, 
                        stringsAsFactors = FALSE, row.names = "COMMON")
cnvG_tab <- read.table("ptenPIK3_CNV.csv", sep = ",", header = TRUE, 
                       stringsAsFactors = FALSE, row.names = "COMMON")
mut_table <- read.table("ptenPIK3_mut.csv", sep = ",", header = TRUE, 
                        stringsAsFactors = FALSE, row.names = "COMMON")

cnv04_hot <- merge(hotSamples, cnv04_tab, by = 0)
cnv04_hot$temp <- "HOT"
cnv04_cold <- merge(coldSamples, cnv04_tab, by = 0)
cnv04_cold$temp <- "COLD"

plot.data.1 <- rbind(cnv04_hot, cnv04_cold)
print("% Cold tumors with biallelic loss of PTEN")
numPTEN.cold <- dim(plot.data.1[plot.data.1$temp == "COLD" & plot.data.1$PTEN == -2, ])[1]
print(numPTEN.cold/numCold)
print("% Hot tumors with biallelic loss of PTEN")
numPTEN.hot <- dim(plot.data.1[plot.data.1$temp == "HOT" & plot.data.1$PTEN == -2, ])[1]
print(numPTEN.hot/numHot)



comp_cnv04 <- data.frame(row.names = c("HOT", "COLD"), 
                         normal.cn = c((numHot- numPTEN.hot),(numCold-numPTEN.cold)), 
                         biallel.loss = c(numPTEN.hot, numPTEN.cold))

fisher.test(comp_cnv04)
  
cnvG_hot <- merge(hotSamples, cnvG_tab, by = 0)
cnvG_hot$temp <- "HOT"
cnvG_cold <- merge(coldSamples, cnvG_tab, by = 0)
cnvG_cold$temp <- "COLD"
plot.data.2 <- rbind(cnvG_hot, cnvG_cold)
print("% Cold tumors with biallelic loss of PTEN")
numPTEN <- dim(plot.data.2[plot.data.2$temp == "COLD" & plot.data.2$PTEN == -2, ])[1]
print(numPTEN/numCold)
print("% Hot tumors with biallelic loss of PTEN")
numPTEN <- dim(plot.data.2[plot.data.2$temp == "HOT" & plot.data.2$PTEN == -2, ])[1]
print(numPTEN/numHot)

write.table(plot.data.2, file = "plot2.csv", sep = ",", quote = FALSE, row.names = F)

mut_hot <- merge(hotSamples, mut_table, by = 0)
mut_hot$temp <- "HOT"
mut_cold <- merge(coldSamples, mut_table, by = 0)
mut_cold$temp <- "COLD"
plot.data.3 <- rbind(mut_hot, mut_cold)
write.table(plot.data.3, file = "plot3.csv", sep = ",", quote = FALSE, row.names = F)
print("% Cold tumors with LOF PTEN mut")
numPTEN <- dim(plot.data.3[plot.data.3$temp == "COLD" & plot.data.3$PTEN == 1, ])[1]
print(numPTEN/numCold)
print("% Hot tumors with LOF PTEN mut")
numPTEN <- dim(plot.data.3[plot.data.3$temp == "HOT" & plot.data.3$PTEN == 1, ])[1]
print(numPTEN/numHot)

print("% Cold tumors with activating PIK3CA mut")
numPTEN <- dim(plot.data.3[plot.data.3$temp == "COLD" & plot.data.3$PIK3CA == 1, ])[1]
print(numPTEN/numCold)
print("% Hot tumors activating PIK3CA mut")
numPTEN <- dim(plot.data.3[plot.data.3$temp == "HOT" & plot.data.3$PIK3CA == 1, ])[1]
print(numPTEN/numHot)



