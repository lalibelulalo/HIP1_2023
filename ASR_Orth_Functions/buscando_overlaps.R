setwd("/home/lalibelulalo/New_Pipe_test/Calothrix_D/PALINDROMES/GCGATCGC/")
Sites336 = read.table("Markov_count_GCGATCGC_Octanuc_.txt",header = TRUE,sep = "\t",row.names = NULL)%>%
  filter(spp == 'Calothrix_sp_336/3')%>%# 'NIES-267'
  filter(obs >= 1)%>% distinct()
Sites3974 = read.table("Markov_count_GCGATCGC_Octanuc_.txt",header = TRUE,sep = "\t",row.names = NULL)%>%
  filter(spp == 'Calothrix_sp_NIES-3974')%>%# 'NIES-267'
  filter(obs >= 1)%>% distinct()
Sites6303 = read.table("Markov_count_GCGATCGC_Octanuc_.txt",header = TRUE,sep = "\t",row.names = NULL)%>%
  filter(spp == 'Calothrix_sp_PCC_6303')%>%# 'NIES-267'
  filter(obs >= 1)%>% distinct()

sum(Sites336$obs)
sum(Sites3974$obs)
sum(Sites6303$obs)



setwd("/home/lalibelulalo/New_Pipe_test/Calothrix_C/PALINDROMES/GCGATCGC/SUBCLADE_ALL/")

Sites336.1291 = read.table("336-6.test.RF1.1291",header = FALSE,sep = "\t",row.names = NULL)
Sites336.1291 = Sites336.1291[,c(1,7,3,4)]
colnames(Sites336.1291) <- c("FILE","PAL","START","END")
write.table(Sites336.1291,file="336-6.test.RF1.1291.format",sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)

Sites336.1290 = read.table("336-6.test.RF1.1290",header = FALSE,sep = "\t",row.names = NULL)
colnames(Sites336.1290) <- c("FILE","PAL","START","END")
write.table(Sites336.1290,file="336-6.test.RF1.1290.format",sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
