Tree = ggtree::read.tree("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/SpeciesTree_rooted.txt")
SppPath <- "/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/GCGATCGC/336-3/"
setwd(SppPath)
HIP1NodeStatus = "All"

RF_list <- list()
RF_title <-list()
RFFiles <- system('ls *rooted.txt | sed \'s/codon_mutations_RF//g\' | sed \'s/.rooted.txt//g\'', intern = TRUE)
RF ="1"
for (RF in RFFiles){
  FileRF = paste0(SppPath,"codon_mutations_RF",RF,".rooted.txt")
  tableRF = read.table(FileRF,header = TRUE,sep = "\t",row.names = NULL)
  TitleRF = tableRF
  #tableRF<-tableRF%>%
  #  filter(SType!='Deletion')
  
  if(HIP1NodeStatus == "Ancestor"){
    ##Filtro dejando solo aquellos casos en los que HIP fue el palíndromo parental
    tableRF<-tableRF%>%
      filter(AncestorType=='SITE')
    TitleRF = tableRF
    PlotsTitle = "Nodos con sitios HIP1 en el nodo parental"
  }
  if(HIP1NodeStatus == "Actual"){
    ##Filtro dejando solo aquellos casos en los que HIP fue el palíndromo parental
    tableRF<-tableRF%>%
      filter(ActualType=='SITE')
    tableRF<-tableRF%>%
      filter(AncestorType=='NoSITE')
    TitleRF = tableRF
    PlotsTitle = "Nodos con sitios HIP1 en el nodo Actual"
  }
  
  
  if (RF==3){
    tableRF = tableRF[ , c(2,43,44)]#,45
  }else{
    tableRF = tableRF[ , c(2,41,42)]#,43
  }
  
  Spps <-unique(tableRF[ , c(1)])
  HIP.Muts <- c()
  for (m in Spps) {
    df<-filter(tableRF,Spp==m)
    NewRow <- c(colSums(df[, -1]),m)
    HIP.Muts <- rbind(HIP.Muts,NewRow)
  }
  ColNames  <- c("Transitions","Transversions","node") #"NoMutations",
  colnames(HIP.Muts) <- ColNames
  rownames(HIP.Muts) <- c(1:length(Spps))
  #ColNames <- unlist(strsplit("GCGATCGC", split=""))
  ancstats2 <- as.data.frame(HIP.Muts)
  RF_list[[RF]] <- ancstats2
  RF_title[[RF]] <- TitleRF
}

plots <- list()
#RF=1
for (RF in 1:3) {
  pies = list()
  for (i in 1:length(RF_list[[RF]][,1])) {
    data <- data.frame(
      group=c(ColNames[1:2]),
      value=as.numeric(as.vector(RF_list[[RF]][1:2][i,]))
    )
    AcualNode <- RF_list[[RF]][3][i,]
    PieTitle <- length(filter(RF_title[[RF]],Spp==AcualNode)[,1])
    pies[[i]] =  ggplot(data, aes(x="", y=value, fill=group)) +
      geom_bar(stat="identity", width=0.1, color="white") +
      coord_polar("y", start=0) +
      theme_void() +
      theme(legend.position="none")+
      ggplot2::scale_fill_brewer(palette="Set2")+
      ggtitle(paste0(PieTitle,"|",PieTitle*8,"|",sum(data$value)))+
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
    ggplot2::scale_fill_brewer(palette="Set2")+
    theme(legend.position='left')+
    guides(fill=guide_legend(title="SMutation Type",title.theme = element_text(color="gray21", size=12, face="bold",hjust = 0.5)))#+
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

#plots<-list(plots[[1]],plots[[2]],plots[[3]])
g <- ggplotGrob(z + theme(legend.position="right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#legend <- gtable::gtable_filter(ggplotGrob(z), 'guide-box')

PLOTS <- gridExtra::grid.arrange(plots[[1]],plots[[2]],plots[[3]],
                                 widths = 1,
                                 right=legend$grobs[[1]],
                                 top=grid::textGrob(PlotsTitle))

i="All"
source("/home/lalibelulalo/HIP1_2023/ASR_Orth_Functions/CodonMutationNodePieCharts3.R")
for(i in c("All","Ancestor","Actual")){
  print(paste0(" --", i))
  res <- try(PLOTS <- Codon_Mutation_Node_Pie_Charts3(Tree = ggtree::read.tree("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/SpeciesTree_rooted.txt"),
                                                      SppPath = paste0("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/GCGATCGC/336-3/"),
                                                      HIP1NodeStatus = i,
                                                      Palindrome = "GCGATCGC",
                                                      TreePlot = p),silent=TRUE)
  if(inherits(res, "try-error")){
    print("  Error. Faltan Datos")
    next
  }
  ggsave(PLOTS, file=paste0("336-3","_",i,"_HIP_mutation_type_tree.png"),width=8, height=6, units="in", scale=1.5)
}




print(paste0(" --", i))
res <- try(PLOTS <- Codon_Mutation_Node_Pie_Charts3(Tree = ggtree::read.tree("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/SpeciesTree_rooted.txt"),
                                                    SppPath = paste0("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/GCGATCGC/336-3/"),
                                                    HIP1NodeStatus = i,
                                                    Palindrome = "GCGATCGC",
                                                    TreePlot = p),silent=TRUE)
if(inherits(res, "try-error")){
  print("  Error. Faltan Datos")
  next
}
ggsave(PLOTS, file=paste0("336-3","_",i,"_HIP_mutation_type_tree.png"),width=8, height=6, units="in", scale=1.5)

