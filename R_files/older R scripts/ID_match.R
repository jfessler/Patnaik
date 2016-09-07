id_key <- read.table("SRR_sampleID.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
rnaseq_id <- colnames(mets_clean)

new_table <- data.frame()
final_key <- data.frame()
count = 0
count2 = 0
count_y = 0
for (id in rnaseq_id){
  if (id %in% id_key$Run_s){
    IDs <- data.frame(id, id_key$SAMPLE.ID[id_key$Run_s == id])
    if (id_key$SAMPLE.ID[id_key$Run_s == id] == ''){
      print(paste0(count2,": Error! ",id))
      count2 = count2 + 1
    }
    else{count_y = count_y + 1
    keep <- data.frame(id, id_key$SAMPLE.ID[id_key$Run_s == id])
    final_key <- rbind(final_key, keep)
    }
    new_table <- rbind(new_table, IDs)
  }
  else {
    print(paste0(count,": ",id," not found."))
    count = count + 1 
    }
}
print(count_y)


