library(lemon)
library(dplyr)
library(ggplot2)
library(ggtree)
Codon_Mutation_Node_Pie_Charts <- function(Tree,SppPath,HIP1NodeStatus,Palindrome,TreePlot) {
  TreeTipsA <- Tree[["tip.label"]]
  TreeTipsB <-  as.data.frame(tidytree::as_tibble(Tree))[,c(2,4)]
  
  if(is.null(TreePlot)==TRUE) { 
    p <- ggtree(Tree,branch.length='none') +
      geom_tiplab(as_ylab=TRUE, color='black')+
      theme(plot.title = element_text(color="black", size=14, face="bold",hjust = 0.5))#+
    #ggplot2::expand_limits(x = 5.5)
  }else(
    p <- TreePlot
  )
  ##-----------------------------------------------------------------------##
  RF_list <- list()
  #ls *rooted.txt | sed 's/codon_mutations_RF//g' | sed 's/.rooted.txt//g'
  #RFS <- sort(unique(system('awk \'{if(NR!=1) {print \$5}}\' Orthologues_Palindrome_sites.AllFrames.SECOND.txt |uniq',intern = TRUE)))
  RFFiles <- system('ls *rooted.txt | sed \'s/codon_mutations_RF//g\' | sed \'s/.rooted.txt//g\'', intern = TRUE)
  #RF =3
  for (RF in RFFiles){
    FileRF = paste0(SppPath,"codon_mutations_RF",RF,".rooted.txt")
    tableRF = read.table(FileRF,header = TRUE,sep = "\t",row.names = NULL)
    
    if(HIP1NodeStatus == "Ancestor"){
      ##Filtro dejando solo aquellos casos en los que HIP fue el palíndromo parental
      tableRF<-tableRF%>%
        filter(AncestorType=='SITE')
      
      PlotsTitle = paste0("Nodos con sitios ",Palindrome," en el nodo parental")
    }
    if(HIP1NodeStatus == "Actual"){
      ##Filtro dejando solo aquellos casos en los que HIP fue el palíndromo parental
      tableRF<-tableRF%>%
        filter(ActualType=='SITE')
      tableRF<-tableRF%>%
        filter(AncestorType=='NoSITE')
      PlotsTitle = paste0("Nodos con sitios ",Palindrome," en el nodo Actual")
    }
    if(HIP1NodeStatus == "All"){
      tableRF=tableRF
      PlotsTitle = paste0("Mutaciones en sitios ",Palindrome," con un octamero CUALESQUIERA del nodo parental al nodo actual")
    }
    
    if (RF==3){
      tableRF = tableRF[ , c(2,22)]
    }else{
      tableRF = tableRF[ , c(2,20)]
    }
    
    if (nrow(tableRF) == 0) stop("FALTAN DATOS.")

    RowNames <- (unique(tableRF[ , c(1)]))
    ColNames <- c('Conservative','ConservativeNoSiteMut','Deletion','NoMutation','NoSynonym','NoSynonymNoSiteMut','Synonym','node')
    ancstats2 <- matrix(0,
                        nrow = length(RowNames),
                        ncol = length(ColNames))
    colnames(ancstats2) <- ColNames
    
    for(j in 1:length(RowNames)){
      tableM.SType <- filter(tableRF,Spp==RowNames[j])
      for (i in 1:(length(ColNames)-1)){
        ancstats2[j,i] <- length(filter(tableM.SType,SType==ColNames[i])[,1])#/length(filter(tableM.SType,Spp==RowNames[j])[,1])
      }
      ancstats2[j,length(ColNames)] <- RowNames[j]
    }
    ancstats2 <- as.data.frame(ancstats2)
    RF_list[[RF]] <- ancstats2
  }
 
  ##-----------------------------------------------------------------------##
  plots <- list()
  for (RF in RFFiles) {
    pies = list()
    for (i in 1:length(RF_list[[RF]][,1])) {
      data <- data.frame(
        group=c('Conservative','ConservativeNoSiteMut','Deletion','NoMutation','NoSynonym','NoSynonymNoSiteMut','Synonym'),
        value=as.numeric(as.vector(RF_list[[RF]][1:7][i,]))
      )
      if (length(RFFiles)==1){
        VJUST = -3
      }
      if (length(RFFiles)==2){
        VJUST = -5
      }
      if (length(RFFiles)==3){
        VJUST = -7
      }
      pies[[i]] =  ggplot(data, aes(x="", y=value, fill=group)) +
        geom_bar(stat="identity", width=0.1, color="white") +
        coord_polar("y", start=0) +
        theme_void() +
        theme(legend.position="none")+
        ggtitle(paste0(sum(data$value)))+
        theme(
          plot.title = element_text(color = "black",
                                    size = 8,
                                    face = "bold",
                                    hjust = ifelse(nchar(sum(data$value))==1,
                                                   -0.1,
                                                   ifelse(nchar(sum(data$value))==2,
                                                          -0.4,
                                                          ifelse(nchar(sum(data$value))==3,
                                                                 -0.9,
                                                                 -1.9))),
                                    vjust = VJUST))
    }
    PiesVector = RF_list[[RF]]$node
    
    for (k in 1:length(TreeTipsA)){
      PiesVector = replace(PiesVector,PiesVector==TreeTipsA[k],as.numeric(as.vector(TreeTipsB[TreeTipsB $label == TreeTipsA[k],][1])))
    }
    names(pies) = PiesVector
    
    z <- ggplot(data, aes(x="", y=value, fill=group)) +
      geom_bar(stat="identity", width=0.1, color="white") +
      coord_polar("y", start=0) +
      theme_void()+
      theme(legend.position='right',legend.key.size = unit(0.4, 'cm'))+
      guides(fill=guide_legend(title="Substitution Type",title.theme = element_text(color="gray21", size=12, face="bold",hjust = 0.5) ))
    #legend <- gtable::gtable_filter(ggplotGrob(z), 'guide-box')
    p <- p +
      ggtitle(paste0("Marco de lectura ",RF))+
      theme(
        plot.title = element_text(color = "gray21", size = 8, face = "bold",vjust = -7))
    #p2 <- inset(p, pies, width=.4, height=.25, hjust=0)
    if(length(RF_list)==1){
      p2 <- inset(p, pies, width=(.6*0.65), height=(.15*0.65), hjust=0,vjust = (-0.3*0.3))
    }
    if(length(RF_list)==2){
      p2 <- inset(p, pies, width=(.6*0.7), height=(.3*0.7), hjust=0,vjust = (-0.3*0.4))
    }
    if(length(RF_list)==3){
      p2 <- inset(p, pies, width=.6, height=.3, hjust=0,vjust = -0.3)
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