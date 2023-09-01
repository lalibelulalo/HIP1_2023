	suppressMessages(source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R"))
	setwd('/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/GCGATCGC/NIES-4071/')
	Nodes.Edges <- Create_Transition_Table_No_Fit(SitesTable = "/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/GCGATCGC/NIES-4071/Orthologues_Palindrome_sites.txt",
		                        EvolutionModel = "F81",
		                        Method = "bayes",
		                        Phylogeny = "/home/lalibelulalo/TESIS/Clados/Calothrix_B/SpeciesTree_rooted.txt",
		                        OrthoPath = "/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/GCGATCGC/NIES-4071/Only_ORTHOLOGUES/",
		                        OutName = "All_RF")
