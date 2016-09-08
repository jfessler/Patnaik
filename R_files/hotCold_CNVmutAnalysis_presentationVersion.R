library(ggplot2)
library(grid)
library(plyr)
library(dplyr)

#############################################
PIK3CA_CN_upperLim <- 100
#############################################

id_tab <- read.table("ID_key.csv", header = TRUE, stringsAsFactors = FALSE, 
                     sep = ",")#row.names = "SRR_ID")

hotSamples <- id_tab[id_tab$SRR_ID %in% hotSRRs,]
rownames(hotSamples) <- hotSamples$sample_ID
hotSamples$SRR_ID = NULL
hotSamples$sample_ID = NULL
hotSamples$temp <- "HOT"
numHot <- dim(hotSamples)[1]

coldSamples <- id_tab[id_tab$SRR_ID %in% coldSRRs,]
rownames(coldSamples) <- coldSamples$sample_ID
coldSamples$SRR_ID = NULL
coldSamples$sample_ID = NULL
coldSamples$temp <- "COLD"
numCold <- dim(coldSamples)[1]

allSamples <- rbind(hotSamples, coldSamples)

mut_table <- read.table("ptenPIK3_mut.csv", sep = ",", header = TRUE, 
                        stringsAsFactors = FALSE, row.names = "COMMON")
colnames(mut_table) <- paste0("mut_",colnames(mut_table))

cn_table <- read.table("ptenPIK3_absCN.csv", sep = ",", header = TRUE, 
                       stringsAsFactors = FALSE, row.names = "COMMON")
colnames(cn_table) <- paste0("cn_",colnames(cn_table))

cn_table$cn_PTEN <- ifelse(cn_table$cn_PTEN < 0.4, 1, 0)
cn_table$cn_PIK3CA <- ifelse(cn_table$cn_PIK3CA > PIK3CA_CN_upperLim, 1, 0)


cn_table <- merge(allSamples, cn_table, by = 0)
rownames(cn_table) <- cn_table$Row.names
cn_table <- cn_table[,c("temp","cn_PTEN","cn_PIK3CA")]

#browser()
allData <- transform(merge(cn_table, mut_table, by = 0), row.names=Row.names, Row.names=NULL)
freqTable <- plyr::count(allData,vars = colnames(allData))

for (col in colnames(freqTable[2:(ncol(freqTable)-2)])){
  if (sum(freqTable[[col]]) == 0){freqTable[[col]] <- NULL}
}

freqTable$Name <- NA
for (row in 1:nrow(freqTable)){
  rowdf <- freqTable[row,2:(ncol(freqTable)-2)]
  name <- paste(colnames(rowdf[which(!rowdf == 0)]),collapse = "+")
  freqTable[row,"Name"] <- name
}
freqTable <- freqTable[freqTable$Name!="",]


resultsList <- list()
statusList <- colnames(freqTable)[2:(ncol(freqTable)-2)]

for (status in statusList){
  countTab <- dplyr::summarise(group_by(freqTable[freqTable[[status]]==1,], temp), 
                               changed = sum(freq) )
  if (dim(countTab)[1] == 2) { # HOT and COLD
    countTab <- as.data.frame(mutate(countTab, 
                                     total = c(numCold, numHot), 
                                     unchanged = (total - changed)))
    } 
  if (dim(countTab)[1] == 1){ # HOT or COLD ONLY
    if (countTab$temp == "HOT"){
      countTab <- rbind(as.data.frame(countTab), data.frame(temp = "COLD", changed = 0))
      countTab <- mutate(countTab, 
                          total = c(numHot, numCold), 
                          unchanged = (total - changed))

      }
    else {
      countTab <- rbind(as.data.frame(countTab), data.frame(temp = "HOT", changed = 0))
      countTab <- mutate(countTab, 
                         total = c(numCold, numHot), 
                         unchanged = (total - changed))
      }
   }
  rownames(countTab) <- countTab$temp
  countTab$temp <- NULL
  countTab$total <- NULL
  ftest <- fisher.test(countTab)
  resultsList <- c(resultsList, paste0(status,": ",round(ftest$p.value,3)))
  print(paste0(status,": ",round(ftest$p.value,3)))
  print(countTab)
}



for (row in 1:nrow(freqTable)){
  if (freqTable[row,"temp"] == "COLD"){
    freqTable[row, "percent"] <- (freqTable[row, "freq"])/numCold
    freqTable[row, "total"] <- numCold
  }
  if (freqTable[row,"temp"] == "HOT"){
    freqTable[row, "percent"] <- (freqTable[row, "freq"])/numHot
    freqTable[row, "total"] <- numHot
  }
}

####

pik3Active_COLD <- sum(freqTable[freqTable$temp == "COLD","freq"])
pik3Active_HOT <- sum(freqTable[freqTable$temp == "HOT","freq"])

notActive_COLD <- numCold - pik3Active_COLD
notActive_HOT <- numHot - pik3Active_HOT

allChanges <- data.frame(c(pik3Active_COLD, pik3Active_HOT), c(notActive_COLD, notActive_HOT))
allChanges <- transform(allChanges, row.names = c("COLD","HOT"))
colnames(allChanges) <- c("activated", "notActivated")
ft <-fisher.test(allChanges)
resultsList <- c(resultsList, paste0("For all activating mutations, p-value=",round(ft$p.value,3)))
print(paste0("for all activating mutations, p-value=",round(ft$p.value,3)))

# Format the labels and calculate their positions
 
freqTable = ddply(freqTable, .(temp), transform, pos = (cumsum(percent) - 0.5 * percent))
freqTable$label = paste0(round(freqTable$percent*100), "%")

p <- ggplot(freqTable, aes(factor(temp), y = percent, fill = Name)) +
  geom_bar(stat = "identity") 

#TITLE <- paste0("K=",cut_k,", cc=",cc_k,", cold:",coldClus," hot:",hotClus,"\n",
#                "# Cold = ",numCold," # Hot = ",numHot)

#p <- p + labs(title = TITLE) + 
#  theme(plot.title = element_text(hjust = 0.5))
  
p <- p+ geom_text( aes(x = temp, y = pos, label = label), size = 4)

totals <- freqTable %>% group_by(temp) %>% dplyr::summarise(total = sum(percent))
totals$t_label<- paste0(round(totals$total*100), "%")

p <- p + geom_text(aes(temp, total + .02, label = t_label, fill = NULL), size = 6, data = totals) + 
  labs(x = "", y = "") + theme_classic()


# Using a manual scale instead of hue
p <- p + scale_fill_manual(values=c("#1E90FF", "#CD8500", "#CDCD00", "#32CD32"), 
                       name="",
                       breaks=c("cn_PTEN", "cn_PTEN+mut_PIK3CA", "cn_PTEN+mut_PIK3CB","cn_PTEN+mut_PTEN"),
                       labels=c("PTEN deletion", "PTEN deletion and PIK3CA activation mutation", "PTEN deletion and PIK3CB activation mutation", "PTEN deletion and PTEN LOF mutation"))

#tiff(paste0(folder,"/",cut_k,"_cc",cc_k,"_Cold",coldClus,"-Hot",hotClus,".tiff"), width = 500, height = 700, compression = "none")
tiff(paste0(folder,"/k500_forPresentation.tiff"), width = 600, height = 700, compression = "none")
plot(p)
dev.off()

sink(file  = paste0(folder,"/pvalueResults.txt"))
print(resultsList)
sink()

source('R_files/tumorSitePlot.R')

