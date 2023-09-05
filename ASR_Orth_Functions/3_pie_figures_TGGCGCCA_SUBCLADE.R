tree2 <- ggtree::groupClade(ggtree::read.tree("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/SpeciesTree_rooted.txt"), c(8,11,12))
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

source("/home/lalibelulalo/HIP1_2023/ASR_Orth_Functions/CodonMutationNodePieCharts.R")
source("/home/lalibelulalo/HIP1_2023/ASR_Orth_Functions/CodonMutationNodePieCharts2.R")
source("/home/lalibelulalo/HIP1_2023/ASR_Orth_Functions/CodonMutationNodePieCharts3.R")

Spps <- c("SUBCLADE")
#PALINDROME = "GCGATCGC"
PALINDROME = "TGGCGCCA"
Output1 = "codon_mutations"
Output2 = "HIP_nuc_mutations"
Output3 = "HIP_mutation_type"

#Spp=1
#i= "Ancestor"
for (Spp in 1:length(Spps)){
  print(Spps[Spp])
  setwd(paste0("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/",PALINDROME,"/",Spps[Spp],"/"))
  print('  PIE 1')
  for(i in c("All","Ancestor","Actual")){
    print(paste0("    --", i))
    res <- try(PLOTS <- Codon_Mutation_Node_Pie_Charts(Tree = ggtree::read.tree("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/SpeciesTree_rooted.txt"),
                                                       SppPath = paste0("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/",PALINDROME,"/",Spps[Spp],"/"),
                                                       HIP1NodeStatus = i,
                                                       Palindrome = PALINDROME,
                                                       TreePlot = p,
                                                       SPP = Spps[Spp]),silent=TRUE)
    if(inherits(res, "try-error")){
      print("      Error. Faltan Datos")
      next
    }
    ggsave(PLOTS, file=paste0(Spps[Spp],"_",i,"_",Output1,"_tree.png"),width=8, height=6, units="in", scale=1.5)
  }
  
  setwd(paste0("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/",PALINDROME,"/",Spps[Spp],"/"))
  print('  PIE 2')
  for(i in c("All","Ancestor","Actual")){
    print(paste0("    --", i))
    res <- try(PLOTS <- Codon_Mutation_Node_Pie_Charts2(Tree = ggtree::read.tree("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/SpeciesTree_rooted.txt"),
                                                        SppPath = paste0("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/",PALINDROME,"/",Spps[Spp],"/"),
                                                        HIP1NodeStatus = i,
                                                        Palindrome = PALINDROME,
                                                        TreePlot = p,
                                                        SPP = Spps[Spp]),silent=TRUE)
    if(inherits(res, "try-error")){
      print("      Error. Faltan Datos")
      next
    }
    ggsave(PLOTS, file=paste0(Spps[Spp],"_",i,"_",Output2,"_tree.png"),width=8, height=6, units="in", scale=1.5)
  }
  
  setwd(paste0("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/",PALINDROME,"/",Spps[Spp],"/"))
  print('  PIE 3')
  for(i in c("All","Ancestor","Actual")){
    print(paste0("    --", i))
    res <- try(PLOTS <- Codon_Mutation_Node_Pie_Charts3(Tree = ggtree::read.tree("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/SpeciesTree_rooted.txt"),
                                                        SppPath = paste0("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/",PALINDROME,"/",Spps[Spp],"/"),
                                                        HIP1NodeStatus = i,
                                                        Palindrome = PALINDROME,
                                                        TreePlot = p,
                                                        SPP = Spps[Spp]),silent=TRUE)
    if(inherits(res, "try-error")){
      print("      Error. Faltan Datos")
      next
    }
    ggsave(PLOTS, file=paste0(Spps[Spp],"_",i,"_",Output3,"_tree.png"),width=8, height=6, units="in", scale=1.5)
  }
}
