library(lemon)
library(dplyr)
library(ggplot2)
library(ggtree)
##-----------------------------------------------------------------------##
#Tree = ggtree::read.tree("/home/lalibelulalo/TESIS/Clados/Callothrix_clade/SpeciesTree_rooted.txt")
Tree = ggtree::read.tree("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/SpeciesTree_rooted.txt")
SppPath <- "/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/TGGCGCCA/NIES-4071/"
setwd(SppPath)
HIP1NodeStatus = "Ancestor"

ggtree::ggtree(Tree,branch.length="none",) +
  ggtree::geom_tiplab(color='firebrick', offset = .14)+
  ggtree::geom_label(ggplot2::aes(label = node),show.legend = FALSE)+
  ggplot2::expand_limits(x = 5.5)  

tree2 <- ggtree::groupClade(Tree, c(8,11,12))
tree3 <- as.data.frame(tidytree::as_tibble(tree2))
tree3[4,5]<- 1
tree3[11,5]<- 1
tree3[12,5]<- 2
tree3[2,5]<- 2
p <- ggtree::ggtree(tidytree::as.treedata(tree3), ggplot2::aes(color=group),branch.length='none') + 
  ggplot2::theme(legend.position='none')+
  ggtree::scale_color_manual(values=c("red","gold1","steelblue"))+
  ggtree::geom_cladelab(node=1, label="336-3", align=TRUE, 
                        geom='label', fill='red',textcolor='white')+
  ggtree::geom_cladelab(node=2, label="NIES-267", align=TRUE, 
                        geom='label', fill='lightblue4',textcolor='white')+
  ggtree::geom_cladelab(node=3, label="NIES-3974", align=TRUE, 
                        geom='label', fill='red',textcolor='white')+
  ggtree::geom_cladelab(node=4, label="PCC_6303", align=TRUE, 
                        geom='label', fill='red',textcolor='white')+
  ggtree::geom_cladelab(node=5, label="PCC_7716", align=TRUE, 
                        geom='label', fill='steelblue',textcolor='white')+
  ggtree::geom_cladelab(node=6, label="NIES-4071", align=TRUE, 
                        geom='label', fill='steelblue',textcolor='white')+
  ggtree::geom_cladelab(node=7, label="NIES-4105", align=TRUE, 
                        geom='label', fill='steelblue',textcolor='white')+
  ggtree::theme(plot.title = ggplot2::element_text(color="black", size=14, face="bold",hjust = 0.5))+
  ggplot2::expand_limits(x = 5.5)

TreePlot <- p
#TreePlot <- NULL
{
TreeTipsA <- Tree[["tip.label"]]
TreeTipsB <-  as.data.frame(tidytree::as_tibble(Tree))[,c(2,4)]
#TreePlot = c()
if(is.null(TreePlot)==TRUE) { 
  p <- ggtree::ggtree(Tree,branch.length='none') +
    geom_tiplab(as_ylab=TRUE, color='black')+
    theme(plot.title = element_text(color="black", size=14, face="bold",hjust = 0.5))#+
    #ggplot2::expand_limits(x = 5.5)
}else(
  p <- TreePlot
)
##-----------------------------------------------------------------------##
RF_list <- list()
RFFiles <- system('ls *rooted.txt | sed \'s/codon_mutations_RF//g\' | sed \'s/.rooted.txt//g\'', intern = TRUE)
#RF ="3"
for (RF in RFFiles){
  FileRF = paste0(SppPath,"codon_mutations_RF",RF,".rooted.txt")
  tableRF = read.table(FileRF,header = TRUE,sep = "\t",row.names = NULL)
  
  if(HIP1NodeStatus == "Ancestor"){
    ##Filtro dejando solo aquellos casos en los que HIP fue el palíndromo parental
    tableRF<-tableRF%>%
      filter(AncestorType=='SITE')
    PlotsTitle = "Nodos con sitios HIP1 en el nodo parental"
  }
  if(HIP1NodeStatus == "Actual"){
    ##Filtro dejando solo aquellos casos en los que HIP fue el palíndromo parental
    tableRF<-tableRF%>%
      filter(ActualType=='SITE')
    tableRF<-tableRF%>%
      filter(AncestorType=='NoSITE')
    PlotsTitle = "Nodos con sitios HIP1 en el nodo Actual"
  }
  
  if (RF==3){
    tableRF = tableRF[ , c(2,22)]
  }else{
    tableRF = tableRF[ , c(2,20)]
  }
  
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
for (RF in 1:3) {
  pies = list()
  for (i in 1:length(RF_list[[RF]][,1])) {
    data <- data.frame(
      group=c('Conservative','ConservativeNoSiteMut','Deletion','NoMutation','NoSynonym','NoSynonymNoSiteMut','Synonym'),
      value=as.numeric(as.vector(RF_list[[RF]][1:7][i,]))
    )
    pies[[i]] =  ggplot(data, aes(x="", y=value, fill=group)) +
      geom_bar(stat="identity", width=0.1, color="white") +
      coord_polar("y", start=0) +
      theme_void() +
      theme(legend.position="none")+
      ggtitle(paste0(sum(data$value)))+
      theme(
        plot.title = element_text(color = "black", size = 8, face = "bold",hjust = ifelse(sum(data$value)<10,-0.1,ifelse(sum(data$value)<=20,-0.5,-1.5)),vjust = -8))
  }
  
  PiesVector = RF_list[[RF]]$node
  
  #for (k in 1:length(TreeTips)) {
  #  PiesVector = replace(PiesVector,PiesVector==TreeTips[k],k)
  #}
  for (k in 1:length(TreeTipsA)){
    PiesVector = replace(PiesVector,PiesVector==TreeTipsA[k],as.numeric(as.vector(TreeTipsB[TreeTipsB $label == TreeTipsA[k],][1])))
  }
  names(pies) = PiesVector
  
  z <- ggplot(data, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=0.1, color="white") +
    coord_polar("y", start=0) +
    theme_void()+
    theme(legend.position='left')+
    guides(fill=guide_legend(title="Substitution Type",title.theme = element_text(color="gray21", size=12, face="bold",hjust = 0.5)))#+
    #ggtitle(paste0(sum(data$value)))+
    #theme(
    #  plot.title = element_text(color = "gray21", size = 22, face = "bold",hjust = 0.50,vjust = - 13))
  
  p<- p +
    ggtitle(paste0("Marco de lectura ",RF))+
    theme(
      plot.title = element_text(color = "gray21", size = 8, face = "bold",vjust = -7))
  p2 <- inset(p, pies, width=.6, height=.3, hjust=0,vjust = -0.3)
  
  plots[[RF]] <- p2
}

plots<-list(plots[[1]],plots[[2]],plots[[3]])
g <- ggplotGrob(z + theme(legend.position="right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#legend <- gtable::gtable_filter(ggplotGrob(z), 'guide-box')

PLOTS <- gridExtra::grid.arrange(plots[[1]],plots[[2]],plots[[3]],
                                 widths = 1,
                                 right=legend$grobs[[1]],
                                 top=grid::textGrob(PlotsTitle))

}
img <- png::readPNG("/home/lalibelulalo/TESIS/Clados/Calothrix_B/Venn/Venn_Calothrix_B-6RELEVANT.png", native = TRUE)

rasterImage(img,2,2,4,4)


source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/CodonMutationNodePieCharts.R")
Spps <- c("336-3",
          "NIES-267",
          "NIES-4105",
          "NIES-3974",
          "PCC_6303",
          "NIES-4071",
          "PCC_7716")
PALINDROME = "GCGATCGC"

#Spp=2
for (Spp in 1:length(Spps)){
  print(Spps[Spp])
  setwd(paste0("/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/",PALINDROME,"/",Spps[Spp],"/"))
  for(i in c("All","Ancestor","Actual")){
    print(paste0(" --", i))
    res <- try(PLOTS <- Codon_Mutation_Node_Pie_Charts(Tree = ggtree::read.tree("/home/lalibelulalo/TESIS/Clados/Calothrix_B/SpeciesTree_rooted.txt"),
                                              SppPath = paste0("/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/",PALINDROME,"/",Spps[Spp],"/"),
                                              HIP1NodeStatus = i,
                                              Palindrome = PALINDROME,
                                              TreePlot = p),silent=TRUE)
    if(inherits(res, "try-error")){
      print("  Error. Faltan Datos")
      next
    }
  ggsave(PLOTS, file=paste0(Spps[Spp],"_",i,"_codon_mutations_tree.png"),width=8, height=6, units="in", scale=1.5)
  }
}



source("/home/lalibelulalo/HIP1_2023/ASR_Orth_Functions/CodonMutationNodePieCharts.R")
PALINDROME = "GCGATCGC"
for(i in c("All","Ancestor","Actual")){
  print(paste0(" --", i))
  res <- try(PLOTS <- Codon_Mutation_Node_Pie_Charts(Tree = ggtree::read.tree("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/SpeciesTree_rooted.txt"),
                                                     SppPath = paste0("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/GCGATCGC/336-3/"),
                                                     HIP1NodeStatus = i,
                                                     Palindrome = "GCGATCGC",
                                                     TreePlot = p),silent=TRUE)
  if(inherits(res, "try-error")){
    print("  Error. Faltan Datos")
    next
  }
  ggsave(PLOTS, file=paste0("336-3","_",i,"_codon_mutations_tree.png"),width=8, height=6, units="in", scale=1.5)
}

