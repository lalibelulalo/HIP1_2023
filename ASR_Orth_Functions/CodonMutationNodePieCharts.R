library(lemon)
library(dplyr)
library(ggplot2)
library(ggtree)
Codon_Mutation_Node_Pie_Charts <- function(Tree,Root,SppPath,HIP1NodeStatus,Palindrome,TreePlot,SPP,Outgroup,Exclude,TaxonNumber,Width,Height,Vjust,VTitle,A,B,C,D) {
  TreeTipsA <- Tree[["tip.label"]]
  #TreeTipsA <- TreeTipsA [ !TreeTipsA  == 'NIES-267']
  TreeTipsB <-  as.data.frame(tidytree::as_tibble(Tree))[,c(2,4)]
  TreeParents <- as.data.frame(tidytree::as_tibble(Tree))[,c(1,2,4)]
  #Root=7
  for (i in Root:length(TreeParents[,1])) {
    TreeParents[i,3] <- TreeParents[i,2]
  }
  TreeParents <- TreeParents[,c(1,3)]
  
  if(is.null(TreePlot)==TRUE) { 
    p <- ggtree(Tree,branch.length='none') +
      geom_tiplab(as_ylab=TRUE, color='black')+
      theme(plot.title = element_text(color="black", size=14, face="bold",hjust = 0.5))+
    #ggplot2::expand_limits(x = 5.5)
      ggtree::geom_label2(ggplot2::aes(label=node, 
                                       subset = !is.na(as.numeric(node)) & as.numeric(node) > TaxonNumber))
  }else(
    p <- TreePlot
  )
  ##-----------------------------------------------------------------------##
  RF_list <- list()
  RF_listALL <- list()
  #ls *rooted.txt | sed 's/codon_mutations_RF//g' | sed 's/.rooted.txt//g'
  #RFS <- sort(unique(system('awk \'{if(NR!=1) {print \$5}}\' Orthologues_Palindrome_sites.AllFrames.SECOND.txt |uniq',intern = TRUE)))
  RFFiles <- system('ls *rooted.txt | sed \'s/codon_mutations_RF//g\' | sed \'s/.rooted.txt//g\'', intern = TRUE)
  #RF =3
  for (RF in RFFiles){
    FileRF = paste0(SppPath,"codon_mutations_RF",RF,".rooted.txt")
    if (Exclude==TRUE){
      tableRF = read.table(FileRF,header = TRUE,sep = "\t",row.names = NULL)%>%
        filter(Spp!=Outgroup)# 'NIES-267'
    }
    if (Exclude==FALSE){
      tableRF = read.table(FileRF,header = TRUE,sep = "\t",row.names = NULL)
    }
    if(HIP1NodeStatus == "Ancestor"){
      ##Filtro dejando solo aquellos casos en los que HIP fue el palíndromo parental
      tableRFALL=tableRF
      tableRF<-tableRF%>%
        filter(AncestorType=='SITE')
      PlotsTitle = paste0("Subset: PARENT.\nTransitions with \"",Palindrome,"\" sequences in the PARENT node.\n")
      #PlotsTitle = paste0("Nodos con sitios ",Palindrome," en el nodo parental")
    }
    if(HIP1NodeStatus == "AncestorO"){
      ##Filtro dejando solo aquellos casos en los que HIP fue el palíndromo parental
      tableRFALL=tableRF
      tableRF<-tableRF%>%
        filter(AncestorType=='SITE')
      tableRF<-tableRF%>%
        filter(ActualType=='NoSITE')
      PlotsTitle = paste0("Subset: ONLY PARENT\nTransitions with \"",Palindrome,"\" sequences ONLY in the PARENT node.\n")
      #PlotsTitle = paste0("Nodos con sitios ",Palindrome," en el nodo parental")
    }
    if(HIP1NodeStatus == "ActualO"){
      ##Filtro dejando solo aquellos casos en los que HIP fue el palíndromo parental
      tableRFALL=tableRF
      tableRF<-tableRF%>%
        filter(ActualType=='SITE')
      tableRF<-tableRF%>%
        filter(AncestorType=='NoSITE')
      PlotsTitle = paste0("Subset: ONLY CHILD\nTransitions with \"",Palindrome,"\" sequences ONLY in the CHILD node.\n")
      #PlotsTitle = paste0("Nodos con sitios ",Palindrome," en el nodo Actual")
    }
    if(HIP1NodeStatus == "Actual"){
      ##Filtro dejando solo aquellos casos en los que HIP fue el palíndromo parental
      tableRFALL=tableRF
      tableRF<-tableRF%>%
        filter(ActualType=='SITE')
      PlotsTitle = paste0("Subset: CHILD.\nTransitions with \"",Palindrome,"\" sequences in the CHILD node.\n")
      #PlotsTitle = paste0("Nodos con sitios ",Palindrome," en el nodo Actual")
    }
    if(HIP1NodeStatus == "All"){
      tableRFALL=tableRF
      tableRF=tableRF
      #Transitions at the "HIP" sites of species 336-3.
      PlotsTitle = paste0("All changes between node pairs at the \"",Palindrome,"\" sites.")
      #PlotsTitle = paste0("Mutaciones en sitios ",Palindrome," con un octamero CUALESQUIERA del nodo parental al nodo actual")
    }
    
    if (RF==3){
      tableRF = tableRF[ , c(2,22)]
      tableRFALL=tableRFALL[ , c(2,22)]
    }else{
      tableRF = tableRF[ , c(2,20)]
      tableRFALL=tableRFALL[ , c(2,20)]
    }
    
    if (nrow(tableRF) == 0) stop("FALTAN DATOS.")

    RowNames <- (unique(tableRF[ , c(1)]))
    RowNamesALL <- (unique(tableRFALL[ , c(1)]))
    #ColNames <- c('Conservative','ConservativeNoSiteMut','Deletion','NoMutation','NoSynonym','NoSynonymNoSiteMut','Synonym','node')
    #ColNames <- c('NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','node')
    ColNames <- c('NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','SynonymNoHIPMutation','node')
    ancstats2 <- matrix(0,
                        nrow = length(RowNames),
                        ncol = length(ColNames))
    colnames(ancstats2) <- ColNames
    
    ancstats2ALL <- matrix(0,
                        nrow = length(RowNamesALL),
                        ncol = length(ColNames))
    colnames(ancstats2ALL) <- ColNames
    
    ancstats5 <- matrix(0,
                           nrow = length(RowNames),
                           ncol = length(ColNames))
    colnames(ancstats5) <- ColNames
    #---
    for(j in 1:length(RowNamesALL)){
      tableM.STypeALL <- filter(tableRFALL,Spp==RowNamesALL[j])
      for (i in 1:(length(ColNames)-1)){
        ancstats2ALL[j,i] <- length(filter(tableM.STypeALL,SType==ColNames[i])[,1])#/length(filter(tableM.SType,Spp==RowNames[j])[,1])
      }
      ancstats2ALL[j,length(ColNames)] <- RowNamesALL[j]
    }
    ancstats2ALL <- as.data.frame(ancstats2ALL)
    dataALL <- data.frame(
      #group=c('Conservative','ConservativeNoSiteMut','Deletion','NoMutation','NoSynonym','NoSynonymNoSiteMut','Synonym'),
      #group=c('NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym'),
      group=c('NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','SynonymNoHIPMutation'),
      value=as.numeric(as.vector(ancstats2ALL[1:8][1,])))
    TotalAll = sum(dataALL$value)
    
    #---
    for(j in 1:length(RowNames)){
      tableM.SType <- filter(tableRF,Spp==RowNames[j])
      for (i in 1:(length(ColNames)-1)){
        ancstats2[j,i] <- length(filter(tableM.SType,SType==ColNames[i])[,1])#/length(filter(tableM.SType,Spp==RowNames[j])[,1])
        ancstats5[j,i] <- round((length(filter(tableM.SType,SType==ColNames[i])[,1])/TotalAll)*100,digits = 2) #signif(7.5678, digits = 3)
      }
      ancstats2[j,length(ColNames)] <- RowNames[j]
      ancstats5[j,length(ColNames)] <- RowNames[j]
    }
    ancstats2 <- as.data.frame(ancstats2)
    ancstats5 <- as.data.frame(ancstats5)
    RF_list[[RF]] <- ancstats2
   #----------- 
    

    #----------- 
    Parent <- c() 
    for (k in 1:length(ancstats2[,1])) {
      for (l in 1:length(TreeParents[,1])) {
        if (ancstats2[k,9] == TreeParents[l,2]) {
          Parent <- rbind(Parent,TreeParents[l,1])
        }
      }
    }
    ReFr <- vector( "numeric" , length = length(ancstats2[,1]) )
    ReFr <- c(rep(paste0("RF",RF), length(ancstats2[,1])))
    TotalRF <- c(rep(TotalAll,length(ancstats2[,1])))
    #---
    TotalVec <- c()
    for (count in 1:length(ancstats2[,1])) {
      #print(ancstats3[count,])
      datax <- data.frame(
      #group=c('Conservative','ConservativeNoSiteMut','Deletion','NoMutation','NoSynonym','NoSynonymNoSiteMut','Synonym'),
      #group=c('NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym'),
      group=c('NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','SynonymNoHIPMutation'),
      value=as.numeric(as.vector(ancstats2[1:8][count,])))
      TotalAllx = sum(datax$value)
      TotalVec <- append(TotalVec,TotalAllx)
    }
    #---
    TotalVec2 <- round((TotalVec/TotalRF[1])*100,digits = 2)
    ancstats3 <- cbind(ancstats2,Parent,ReFr,TotalVec,TotalRF)
    ancstats7 <- cbind(ancstats5,Parent,ReFr,TotalVec2,TotalRF)
    
    #ancstats3 <- ancstats3[,c(10,9,8,1,2,3,4,5,6,7,11,12)]
    #ancstats7 <- ancstats7[,c(10,9,8,1,2,3,4,5,6,7,11,12)]
    ancstats3 <- ancstats3[,c(11,10,9,1,2,3,4,5,6,7,8,12,13)]
    ancstats7 <- ancstats7[,c(11,10,9,1,2,3,4,5,6,7,8,12,13)]
    
    #colnames(ancstats3) <- c('RF','Parent','Node','NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','Total','SitesNumber')
    #colnames(ancstats7) <- c('RF','Parent','Node','NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','Total','SitesNumber')
    colnames(ancstats3) <- c('RF','Parent','Node','NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','SynonymNoHIPMutation','Total','SitesNumber')
    colnames(ancstats7) <- c('RF','Parent','Node','NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','SynonymNoHIPMutation','Total','SitesNumber')
    
    
    write.table(ancstats3,file=paste0("Ancstats.",HIP1NodeStatus,".",RF,".RF"),sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
    write.table(ancstats7,file=paste0("Ancstats.",HIP1NodeStatus,".",RF,".percent"),sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
  }
  output = paste0("Ancstats.",HIP1NodeStatus,".All_RF.txt")
  output7 = paste0("Ancstats.percent.",HIP1NodeStatus,".All_RF.txt")
  system(sprintf("cat *.RF >%s", output))
  system(sprintf("cat *.percent >%s", output7))
  system('rm *.RF')
  system('rm *.percent')
  ancstats4 <- read.table(output,header = FALSE,sep = "\t",row.names = NULL)
  ancstats6 <- read.table(output7,header = FALSE,sep = "\t",row.names = NULL)
  #colnames(ancstats4) <- c('RF','Parent','Node','NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','Total','SitesNumber')
  #colnames(ancstats6) <- c('RF','Parent','Node','NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','Total','SitesNumber')
  colnames(ancstats4) <- c('RF','Parent','Node','NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','SynonymNoHIPMutation','Total','SitesNumber')
  colnames(ancstats6) <- c('RF','Parent','Node','NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','SynonymNoHIPMutation','Total','SitesNumber')
  
  write.table(ancstats4,file=output,sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
  write.table(ancstats6,file=output7,sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
  
  ##-----------------------------------------------------------------------##
  plots <- list()
  for (RF in RFFiles) {
    pies = list()
    for (i in 1:length(RF_list[[RF]][,1])) {
      data <- data.frame(
        #group=c('Conservative','ConservativeNoSiteMut','Deletion','NoMutation','NoSynonym','NoSynonymNoSiteMut','Synonym'),
        #group=c('NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym'),
        group=c('NoSynonymConservative','ConservativeNoHIPMutation','Deletion','NoMutation','NoSynonym','NoSynonymNoHIPMutation','Synonym','SynonymNoHIPMutation'),
        value=as.numeric(as.vector(RF_list[[RF]][1:8][i,]))
      )
      if (length(RFFiles)==1){
        VJUST = -3
      }
      if (length(RFFiles)==2){
        VJUST = -4
      }
      if (length(RFFiles)==3){
        VJUST = -5
      }
      VALUES <- c('NoMutation'='chartreuse2',
                  'Synonym'='dodgerblue',
                  'NoSynonymConservative'='goldenrod1',
                  'NoSynonym'='magenta',
                  'Deletion'='red',
                  'ConservativeNoHIPMutation'='black',
                  'NoSynonymNoHIPMutation'='darkslategray',
                  'SynonymNoHIPMutation'='darkorange4')
      TYPES <- c('NoMutation',
                 'Synonym',
                 'NoSynonymConservative',
                 'NoSynonym',
                 'Deletion',
                 'ConservativeNoHIPMutation',
                 'NoSynonymNoHIPMutation',
                 'SynonymNoHIPMutation')
      pies[[i]] =  ggplot(data = data, aes(x = "", y = -value, 
                                           fill = reorder(group, -value))) + 
        geom_bar(stat = "identity") +
        coord_polar("y") +
        theme_void()+
        #theme(legend.position='right',legend.key.size = unit(0.4, 'cm'))+
        #guides(fill=guide_legend(title="Substitution Type",title.theme = element_text(color="gray21", size=10, face="bold",hjust = 0.5) ))+
        scale_fill_manual("",values=VALUES,
                          breaks=TYPES,
                          labels=TYPES)+
        #ggplot(data, aes(x="", y=value, fill=group)) +
        #geom_bar(stat="identity", width=0.1, color="white") +
        #coord_polar("y", start=0) +
        #theme_void() +
        theme(legend.position="none")+
        ggtitle(paste0(sum(data$value)))+
        theme(
          plot.title = element_text(color = "black",
                                    size = 8,
                                    face = "bold",
                                    hjust = ifelse(nchar(sum(data$value))==1,
                                                   A,#-0.2
                                                   ifelse(nchar(sum(data$value))==2,
                                                          B,#-0.6
                                                          ifelse(nchar(sum(data$value))==3,
                                                                 C,#-1.2
                                                                 D))), #-2.6
                                    vjust = VTitle))
    }
    PiesVector = RF_list[[RF]]$node
    
    for (k in 1:length(TreeTipsA)){
      PiesVector = replace(PiesVector,PiesVector==TreeTipsA[k],as.numeric(as.vector(TreeTipsB[TreeTipsB $label == TreeTipsA[k],][1])))
    }
    names(pies) = PiesVector
    
    # z <- ggplot(data, aes(x="", y=value, fill=group)) +
    #   geom_bar(stat="identity")+
    #   scale_fill_manual("",values=c('NoMutation'='green','Synonym'='cyan','NoSynonymConservative'='yellow','NoSynonym'='magenta','Deletion'='red','ConservativeNoHIPMutation'='black','NoSynonymNoHIPMutation'='darkslategray','SynonymNoHIPMutation'='darkorange4'),
    #                     breaks=c('NoMutation','Synonym','NoSynonymConservative','NoSynonym','Deletion','ConservativeNoHIPMutation','NoSynonymNoHIPMutation','SynonymNoHIPMutation'),
    #                     labels=c('NoMutation','Synonym','NoSynonymConservative','NoSynonym','Deletion','ConservativeNoHIPMutation','NoSynonymNoHIPMutation','SynonymNoHIPMutation'))+
    #   coord_polar("y", start=0) +
    #   theme_void()+
    #   theme(legend.position='right',legend.key.size = unit(0.4, 'cm'))+
    #   guides(fill=guide_legend(title="Substitution Type",title.theme = element_text(color="gray21", size=10, face="bold",hjust = 0.5) ))
    #legend <- gtable::gtable_filter(ggplotGrob(z), 'guide-box')

      

    z <- ggplot(data = data, aes(x = "", y = -value, 
                               fill = reorder(group, -value))) + 
      geom_bar(stat = "identity") +
      coord_polar("y") +
      theme_void()+
      theme(legend.position='right',legend.key.size = unit(0.4, 'cm'))+
      guides(fill=guide_legend(title="Substitution Type",title.theme = element_text(color="gray21", size=10, face="bold",hjust = 0.5) ))+
      scale_fill_manual("",values=VALUES,
                                     breaks=TYPES,
                                     labels=TYPES)
    legend <- gtable::gtable_filter(ggplotGrob(z), 'guide-box') 
      
    p <- p +
      ggtitle(paste0("Reading Frame ",RF))+
      #ggtitle(paste0("Marco de lectura ",RF))+
      theme(
        plot.title = element_text(color = "gray21", size = 8, face = "bold",vjust = -7))
    #p2 <- inset(p, pies, width=.4, height=.25, hjust=0)
    if(length(RF_list)==1){
      p2 <- inset(p, pies, x='branch', width=Width, height=Height, hjust=0,vjust = Vjust)
    }
    if(length(RF_list)==2){
      p2 <- inset(p, pies, x='branch', width=Width, height=Height, hjust=0,vjust = Vjust)
    }
    if(length(RF_list)==3){
      p2 <- inset(p, pies, x='branch', width=Width, height=Height, hjust=0,vjust = Vjust) ### <========================================
    }
    
    
    plots[[RF]] <- p2
  }
  
  g <- ggplotGrob(z + theme(legend.position="right"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  
  if (length(RF_list)==1){
    pl <- gridExtra::grid.arrange(plots[[1]],
                                  widths = 1,
                                  right=legend$grobs[[1]],
                                  top=grid::textGrob(PlotsTitle))
  }
  if (length(RF_list)==2){
    pl <- gridExtra::grid.arrange(plots[[1]],plots[[2]],
                                  widths = 1,
                                  right=legend$grobs[[1]],
                                  top=grid::textGrob(PlotsTitle))
  }
  if (length(RF_list)==3){
    pl <- gridExtra::grid.arrange(plots[[1]],plots[[2]],plots[[3]],
                                  widths = 1,
                                  right=legend$grobs[[1]],
                                  top=grid::textGrob(PlotsTitle))
  }
  return(pl)
}