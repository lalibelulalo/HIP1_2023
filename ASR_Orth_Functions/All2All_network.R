AA <- unique(AACountsM1$AA)
AA[AA == "----"] <- "****"
#AACountsM1$AA
LINKS <- c()
for(i in 1:length(AA)){
  for (j in 1:length(AA)){
    #print(paste0(AA[i]," -- ",AA[j]))
    LINKS <- rbind(LINKS,c(AA[i],AA[j]))
  }
}


library(Biostrings)
data(BLOSUM62)

for (i in 1:length(LINKS[,1])){
  pep1 = seqinr::s2c(LINKS[i,][1])
  pep2 = seqinr::s2c(LINKS[i,][2])
  for (j in 1:4){
    Score <- BLOSUM62[pep1[j],pep2[j]]
    while (Score>=0) {
      print(Score)
    }
  }
}



