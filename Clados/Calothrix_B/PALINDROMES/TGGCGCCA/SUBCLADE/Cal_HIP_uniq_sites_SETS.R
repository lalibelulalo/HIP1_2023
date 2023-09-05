library(dplyr)
source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R")
setwd("/home/lalibelulalo/TESIS/Clados/Geminocystis_clade/")


clades <- c(#"Geminocystis_clade",
            #"Callothrix_clade",
            "Pseudoanabaena_clade",
            #"Clado_A18-40",
            #"SynechococcusCyanobium_clade",
            #"Cyanobacterium_clade",
            "Thermosynechococcus_clade")
for (j in clades){
  #print(j)
  spps <- system(paste0("ls /home/lalibelulalo/TESIS/Clados/",j,"/PALINDROMES/GCGATCGC/"), intern = TRUE)
  for (i in 1:length(spps)){
    setwd(paste0("/home/lalibelulalo/TESIS/Clados/",j,"/PALINDROMES/GCGATCGC/",spps[i],"/"))
    Create_Transition_Table_No_Fit(SitesTable = "Orthologues_Palindrome_sites.txt",
                                   EvolutionModel = "F81",
                                   Method = "bayes",
                                   Phylogeny = "../../../SpeciesTree_rooted.txt",
                                   OrthoPath = "Only_ORTHOLOGUES/")
  }
}

j = "Thermosynechococcus_clade"
spps <- system(paste0("ls /home/lalibelulalo/TESIS/Clados/",j,"/PALINDROMES/CAGGCCTG/"), intern = TRUE)
for (i in 1:length(spps)){
  setwd(paste0("/home/lalibelulalo/TESIS/Clados/",j,"/PALINDROMES/CAGGCCTG/",spps[i],"/"))
  Create_Transition_Table_No_Fit(SitesTable = "Orthologues_Palindrome_sites.txt",
                                 EvolutionModel = "F81",
                                 Method = "bayes",
                                 Phylogeny = "../../../SpeciesTree_rooted.txt",
                                 OrthoPath = "Only_ORTHOLOGUES/")
}

#------------------------------------
setwd("/home/lalibelulalo/HIP1_2023/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/SUBLCLADE/")
SitesTableA = "/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/GCGATCGC/336-3/Orthologues_Palindrome_sites.txt"
SitesTableB = "/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/GCGATCGC/NIES-3974/Orthologues_Palindrome_sites.txt"
SitesTableC = "/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/GCGATCGC/PCC_6303/Orthologues_Palindrome_sites.txt"


SitesA <- read.table(SitesTableA, sep = "\t", header = TRUE)
SitesA$MergedCoordsOrth <- paste(SitesA$START,"-",SitesA$END,"-",SitesA$FILE)

SitesB <- read.table(SitesTableB, sep = "\t", header = TRUE)
SitesB$MergedCoordsOrth <- paste(SitesB$START,"-",SitesB$END,"-",SitesB$FILE)

SitesC <- read.table(SitesTableC, sep = "\t", header = TRUE)
SitesC$MergedCoordsOrth <- paste(SitesC$START,"-",SitesC$END,"-",SitesC$FILE)


Full_AB <- full_join(SitesA, SitesB, by = "MergedCoordsOrth")
Full_AB <- Full_AB %>% distinct()
Full_AB <- Full_AB[ , c(1,2,3,4,5)]

Full_BC <- full_join(SitesB, SitesC, by = "MergedCoordsOrth")
Full_BC <- Full_BC %>% distinct()
Full_BC <- Full_BC[ , c(1,2,3,4,5)]

Full_AC <- full_join(SitesA, SitesC, by = "MergedCoordsOrth")
Full_AC <- Full_AC %>% distinct()
Full_AC <- Full_AC[ , c(1,2,3,4,5)]

Cal_336 <- SitesA[!(SitesA$MergedCoordsOrth %in% Full_BC$MergedCoordsOrth),]
Cal_3974 <- SitesB[!(SitesB$MergedCoordsOrth %in% Full_AC$MergedCoordsOrth),]
Cal_6303 <- SitesC[!(SitesC$MergedCoordsOrth %in% Full_AB$MergedCoordsOrth),]

Cal_HIP_uniq_sites <- Cal_336 %>%
  bind_rows(Cal_3974) %>%
  bind_rows(Cal_6303)

Cal_HIP_uniq_sites <- Cal_HIP_uniq_sites[ , c(1,2,3,4)]
write.table(Cal_HIP_uniq_sites,file="Cal_HIP_uniq_sites.txt",sep="\t",row.names = FALSE,col.names = TRUE)


source("/home/lalibelulalo/HIP1_2023/ASR_Orth_Functions/NodeAndEdges.R")
setwd("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/GCGATCGC/SUBLCLADE/")
Nodes.Edges <- Create_Transition_Table_No_Fit(SitesTable = "Cal_HIP_uniq_sites.txt",
                                              EvolutionModel = "F81",
                                              Method = "bayes",
                                              Phylogeny = "../../../SpeciesTree_rooted.txt",
                                              OrthoPath = "Only_ORTHOLOGUES/",
                                              OutName = "All_RF")

#------------------------------------------------------------------------------------------------------------------------------------------------
setwd("/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/TGGCGCCA/SUBCLADE/")
SecondTableA = "/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/TGGCGCCA/NIES-4105/Orthologues_Palindrome_sites.AllFrames.SECOND.txt"
SecondTableB = "/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/TGGCGCCA/NIES-4071/Orthologues_Palindrome_sites.AllFrames.SECOND.txt"
SecondTableC = "/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/TGGCGCCA/PCC_7716/Orthologues_Palindrome_sites.AllFrames.SECOND.txt"

SecondA <- read.table(SecondTableA, sep = "\t", header = TRUE)
SecondA$MergedCoordsOrth <- paste(SecondA$FILE,"-",SecondA$Spp,"-",SecondA$START)

SecondB <- read.table(SecondTableB, sep = "\t", header = TRUE)
SecondB$MergedCoordsOrth <- paste(SecondB$FILE,"-",SecondB$Spp,"-",SecondB$START)

SecondC <- read.table(SecondTableC, sep = "\t", header = TRUE)
SecondC$MergedCoordsOrth <- paste(SecondC$FILE,"-",SecondC$Spp,"-",SecondC$START)


Full_AB <- full_join(SecondA, SecondB, by = "MergedCoordsOrth")
Full_AB <- Full_AB %>% distinct()

Full_BC <- full_join(SecondB, SecondC, by = "MergedCoordsOrth")
Full_BC <- Full_BC %>% distinct()

Full_AC <- full_join(SecondA, SecondC, by = "MergedCoordsOrth")
Full_AC <- Full_AC %>% distinct()


Cal_336 <- SecondA[!(SecondA$MergedCoordsOrth %in% Full_BC$MergedCoordsOrth),]
Cal_336 <- Cal_336[ , c(1:8)]

Cal_3974 <- SecondB[!(SecondB$MergedCoordsOrth %in% Full_AC$MergedCoordsOrth),]
Cal_3974 <- Cal_3974[ , c(1:8)]

Cal_6303 <- SecondC[!(SecondC$MergedCoordsOrth %in% Full_AB$MergedCoordsOrth),]
Cal_6303 <- Cal_6303[ , c(1:8)]

Cal_HIP_uniq_sites <- Cal_336 %>%
  bind_rows(Cal_3974) %>%
  bind_rows(Cal_6303)

#Cal_HIP_uniq_sites <- Cal_HIP_uniq_sites[ , c(1,2,3,4)]
write.table(Cal_HIP_uniq_sites,file="Cal_HIP_uniq_sites.AllFrames.SECOND.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)


RFS <- length(unique(system('awk \'{if(NR!=1) {print $5}}\' Cal_HIP_uniq_sites.AllFrames.SECOND.txt |uniq',intern = TRUE)))

for(RF in 1:RFS){ ## Para los tres marcos de lectura (RF)
  CreateCodonMutationsTable(Tree = '../../../SpeciesTree_rooted.txt',
                            Second = 'Cal_HIP_uniq_sites.AllFrames.SECOND.txt',
                            TaxonNumber = 7,
                            Root = 8,
                            RF = RF)
  system(paste0('python3 ../../../../CodonMutations.py RF',RF,'.RECONSTRUCTION.rooted.parents.txt ',RF,' codon_mutations_RF',RF,'.rooted.txt TGGCGCCA'))
}

system('cat *.rooted.txt >codon_mutations_AllF.rooted.txt')

#---------------------------------------------------------
setwd("/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/")
PalsitesA <- read.table("/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Orthologues_Palindrome_sites.txt", sep = "\t", header = TRUE)
PalsitesB <- read.table("/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Orthologues_Palindrome_sites.AllFrames.FIRST.txt", sep = "\t", header = TRUE)

PalsitesA$MergedCoordsOrth <- paste(PalsitesA$FILE,"-",PalsitesA$START,"-",PalsitesA$END)
PalsitesB$MergedCoordsOrth <- paste(PalsitesB$FILE,"-",PalsitesB$START,"-",PalsitesB$END)


AllFrames <- PalsitesB[!(PalsitesB$MergedCoordsOrth %in% PalsitesA$MergedCoordsOrth),]

AllFrames <- dplyr::inner_join(PalsitesA, PalsitesB, by = "MergedCoordsOrth")
AllFrames <- AllFrames[ , c(1,2,3,4,7,10)]
colnames(AllFrames)<- c("FILE","PAL","START","END","SPP","RF")

AllFrames <- AllFrames%>%
  filter(SPP=="336-3")
AllFrames <- AllFrames[ , c(1,2,3,4,6)]

for (i in 1:3){
  RF <- AllFrames%>%
    filter(RF==paste0(i))
  write.table(RF,file=paste0("Orthologues_Palindrome_sites.RF",i),sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
}

RF1 <- AllFrames%>%
  filter(RF=="1")
RF2 <- AllFrames%>%
  filter(RF=="2")
RF3 <- AllFrames%>%
  filter(RF=="3")


source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R")
source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/CreateTransitionTableNoFitCodon.R")
setwd("/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/")
for (RF in 1:3){
  Nodes.Edges <- Create_Transition_Table_No_Fit(SitesTable = paste0("/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Orthologues_Palindrome_sites.RF",RF),
                                                EvolutionModel = "F81",
                                                Method = "bayes",
                                                Phylogeny = "/home/lalibelulalo/TESIS/Clados/Callothrix_clade/SpeciesTree_rooted.txt",
                                                OrthoPath = "/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Only_ORTHOLOGUES/",
                                                OutName = paste0("RF",RF))
}

RF= 3
Nodes.Edges <- Create_Transition_Table_No_Fit_Codon(SitesTable = paste0("/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Orthologues_Palindrome_sites.RF",RF),
                                              EvolutionModel = "F81",
                                              Method = "bayes",
                                              Phylogeny = "/home/lalibelulalo/TESIS/Clados/Callothrix_clade/SpeciesTree_rooted.txt",
                                              OrthoPath = "/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Only_ORTHOLOGUES/",
                                              ReadingFrame = RF,
                                              OutName = paste0("RF",RF))
#---

Transition_table <- read.table("/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/RF3.TRANSICIONES.AA", sep = "\t", header = TRUE )
Network <- Create_Network(Transitions = Transition_table,
                          MinimumCounts = 7,
                          Loops = TRUE)

nodes = Network[[1]]
edges = Network[[2]]

## Este es otro grafo de fuerzas con peso en sus vertices y te muestra las conceciones de cada nodo de manera mas visual
nodes_d3 <- mutate(nodes, id = id - 1)
edges_d3 <- mutate(edges, from = from - 1, to = to - 1)

#Sys.setenv("OPENSSL_CONF"="/dev/null")
## Este ultimo grafo muestra el grafo anterior pero de una forma mas analizable.
networkD3::sankeyNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
                         NodeID = "label", Value = "weight", fontSize = 16, unit = "TIME(s)")

#---
Directed_Transition_table <- read.table("/home/lalibelulalo/TESIS/Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/RF2.TRANSICIONES_DIRECCION.AA", sep = "\t", header = TRUE )
Directed_Network <- Create_Directed_Network(DirectedTransitions = Directed_Transition_table,
                                            Direction = "10--9",
                                            Weight = 2)
nodesL = Directed_Network[[1]]
edgesL = Directed_Network[[2]]

nodes_d3L <- mutate(nodesL, id = id - 1)
edges_d3L <- mutate(edgesL, from = from - 1, to = to - 1)

networkD3::sankeyNetwork(Links = edges_d3L, Nodes = nodes_d3L, Source = "from", Target = "to", 
                         NodeID = "label", Value = "weight", fontSize = 16, unit = "TIME(s)")
