library(lemon)
library(dplyr)
library(ggplot2)
library(ggtree)
Codon_Mutation_Node_Pie_Charts2 <- function(Tree,SppPath,HIP1NodeStatus,Palindrome,TreePlot) {
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
  RF_title <-list()
  #ls *rooted.txt | sed 's/codon_mutations_RF//g' | sed 's/.rooted.txt//g'
  #RFS <- sort(unique(system('awk \'{if(NR!=1) {print \$5}}\' Orthologues_Palindrome_sites.AllFrames.SECOND.txt |uniq',intern = TRUE)))
  RFFiles <- system('ls *rooted.txt | sed \'s/codon_mutations_RF//g\' | sed \'s/.rooted.txt//g\'', intern = TRUE)
  #RF =3
  for (RF in RFFiles){
    FileRF = paste0(SppPath,"codon_mutations_RF",RF,".rooted.txt")
    tableRF = read.table(FileRF,header = TRUE,sep = "\t",row.names = NULL)
    TitleRF = tableRF
    
    if(HIP1NodeStatus == "Ancestor"){
      ##Filtro dejando solo aquellos casos en los que HIP fue el palíndromo parental
      tableRF<-tableRF%>%
        filter(AncestorType=='SITE')
      TitleRF = tableRF
      PlotsTitle = paste0("Nodos con sitios ",Palindrome," en el nodo parental")
    }
    if(HIP1NodeStatus == "Actual"){
      ##Filtro dejando solo aquellos casos en los que HIP fue el palíndromo parental
      tableRF<-tableRF%>%
        filter(ActualType=='SITE')
      tableRF<-tableRF%>%
        filter(AncestorType=='NoSITE')
      TitleRF = tableRF
      PlotsTitle = paste0("Nodos con sitios ",Palindrome," en el nodo Actual")
    }
    if(HIP1NodeStatus == "All"){
      tableRF=tableRF
      TitleRF = tableRF
      PlotsTitle = paste0("Mutaciones en sitios ",Palindrome," con un octamero CUALESQUIERA del nodo parental al nodo actual")
    }
    
    if (RF==3){
      tableRF = tableRF[ , c(2,27,28,29,30,31,32,33,34)]
    }else{
      tableRF = tableRF[ , c(2,25,26,27,28,29,30,31,32)]
    }
    
    if (nrow(tableRF) == 0) stop("FALTAN DATOS.")
    
    Spps <-unique(tableRF[ , c(1)])
    HIP.Muts <- c()
    for (m in Spps) {
      df<-filter(tableRF,Spp==m)
      NewRow <- c(colSums(df[, -1]),m)
      HIP.Muts <- rbind(HIP.Muts,NewRow)
    }
    ColNames  <- c("Nuc1","Nuc2","Nuc3","Nuc4","Nuc5","Nuc6","Nuc7","Nuc8","node")
    colnames(HIP.Muts) <- ColNames
    rownames(HIP.Muts) <- c(1:length(Spps))
    #ColNames <- unlist(strsplit("GCGATCGC", split=""))
    ancstats2 <- as.data.frame(HIP.Muts)
    RF_list[[RF]] <- ancstats2
    RF_title[[RF]] <- TitleRF
  }
  
  ##-----------------------------------------------------------------------##
  plots <- list()
  for (RF in RFFiles) {
    pies = list()
    for (i in 1:length(RF_list[[RF]][,1])) {
      data <- data.frame(
        group=c(ColNames[1:8]),
        value=as.numeric(as.vector(RF_list[[RF]][1:8][i,]))
      )
      AcualNode <- RF_list[[RF]][9][i,]
      PieTitle <- length(filter(RF_title[[RF]],Spp==AcualNode)[,1])
      PieTitle <- paste0(PieTitle,"|",sum(data$value))
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
        ggplot2::scale_fill_brewer(palette="Set2")+
        ggtitle(PieTitle)+
        theme(
          plot.title = element_text(color = "black",
                                    size = 8,
                                    face = "bold",
                                    hjust = ifelse(nchar(PieTitle)==3,
                                                   -0.8,
                                                   ifelse(nchar(PieTitle)==4,
                                                          -2.0,
                                                          ifelse(nchar(PieTitle)==5,
                                                                 -4.0,
                                                                 ifelse(nchar(PieTitle)==6,
                                                                        -20.5,
                                                                        ifelse(nchar(PieTitle)==7,
                                                                               9.8,
                                                                               ifelse(nchar(PieTitle)==8,
                                                                                      4.7,
                                                                                      ifelse(nchar(PieTitle)==9,
                                                                                             3.3,-10))))))),
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
      ggplot2::scale_fill_brewer(palette="Set2")+
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