setwd("/home/lalibelulalo/HIP1_2023/figures/")
Tree = ggtree::read.tree("/home/lalibelulalo/HIP1_2023/Clados/Callothrix_clade/SpeciesTree_rooted.txt")
tree2 <- ggtree::groupClade(Tree, c(7,10,11))
tree3 <- as.data.frame(tidytree::as_tibble(tree2))

tree3[4,5]<- 3

p<- ggtree::ggtree(tidytree::as.treedata(tree3), ggplot2::aes(color=group),branch.length='none') + 
  ggplot2::theme(legend.position='none')+
  ggtree::scale_color_manual(values=c("red","gold1","steelblue"))+
  ggtree::geom_label(ggplot2::aes(label = node),show.legend = TRUE)+
  ggtree::geom_cladelab(node=1, label="336-3", align=TRUE, 
                        geom='label', fill='red')+
  ggtree::geom_cladelab(node=2, label="NIES-3974", align=TRUE, 
                        geom='label', fill='red')+
  ggtree::geom_cladelab(node=3, label="PCC_6303", align=TRUE, 
                        geom='label', fill='red',)+
  ggtree::geom_cladelab(node=4, label="PCC_7716", align=TRUE, 
                        geom='label', fill='steelblue')+
  ggtree::geom_cladelab(node=5, label="NIES-4071", align=TRUE, 
                        geom='label', fill='steelblue')+
  ggtree::geom_cladelab(node=6, label="NIES-4105", align=TRUE, 
                        geom='label', fill='steelblue')+
  ggplot2::expand_limits(x = 6)

ggplot2::ggsave(p, file=paste0("Calothrix_6_species_node transition_yellow_tree.png"),width=6, height=4, units="in", scale=1.5)

getwd()
