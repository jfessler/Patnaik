freqTable <- transform(freqTable, row.names=temp, temp=NULL)

rownames(freqTable) <- 
  fisher.test(freqTable[,1:2])

mfreqTable <- melt(data = freqTable, id.vars = "temp", variable_name = "Name")



############
############
############

byTemp <- group_by(allData, temp)
ptenResults <- (summarize(byTempPTEN, count = n()))

byTempPIK3CA <- group_by(cn_table, temp, cn_PIK3CA)
pik3caResults <- summarize(byTempPIK3CA, count(n))


summaryTab <- summarize(group_by(cn_table, temp, PTEN, PIK3CA), count = n())


qplot(factor(temp), data=cn_table, geom = "bar", fill = factor(PIK3CA))

k <- ggplot(cn_table, aes(temp, fill=PIK3CA))
k + geom_bar()

library(scales)
(p <- ggplot(mdfr, aes(category, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = percent)
)

p <- ggplot(cn_table, aes(temp, fill = mut_PTEN)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = percent)


cn_hot <- merge(hotSamples, cn_table, by = 0)
cn_cold <- merge(coldSamples, cn_table, by = 0)

comp_cnv04 <- data.frame(row.names = c("HOT", "COLD"), 
                         normal.cn = c((numHot- numPTEN.hot),(numCold-numPTEN.cold)), 
                         biallel.loss = c(numPTEN.hot, numPTEN.cold))

fisher.test(comp_cnv04)

