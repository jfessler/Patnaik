
tumorSite <- read.table("tumorSite.csv", sep = ",", header = T, row.names = "Patient.ID")

hotTumors <- merge(hotSamples, tumorSite, by = 0)
hotTumors$temp <- "HOT"
numHot <- nrow(hotTumors)

coldTumors <- merge(coldSamples, tumorSite, by = 0)
coldTumors$temp <- "COLD"
numCold <- nrow(coldTumors)

allTumors <- transform(rbind(hotTumors, coldTumors), row.names = Row.names, Row.names = NULL)
tumorFreq <- plyr::count(allTumors)

for (row in 1:nrow(tumorFreq)){
  if (tumorFreq[row,"temp"] == "COLD"){
    tumorFreq[row, "percent"] <- (tumorFreq[row, "freq"])/numCold
  }
  if (tumorFreq[row,"temp"] == "HOT"){
    tumorFreq[row, "percent"] <- (tumorFreq[row, "freq"])/numHot
  }
}

tumorFreq = ddply(tumorFreq, .(temp), transform, pos = (cumsum(percent) - 0.5 * percent))
tumorFreq$label = paste0(round(tumorFreq$percent*100), "%")

p <- ggplot(tumorFreq, aes(factor(temp), y = percent, fill = Tumor.Site)) +          
      geom_bar(stat = "identity") + scale_y_continuous() +
      ggtitle(paste0("# Cold = ",numCold," # Hot = ",numHot))
      
p <- p+ geom_text( aes(x = temp, y = pos, label = label), size = 4)
p <- p + geom_text( aes(x = temp, y = pos + .03, label = freq), size = 2)


