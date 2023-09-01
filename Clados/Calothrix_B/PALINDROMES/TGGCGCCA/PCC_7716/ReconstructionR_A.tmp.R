	suppressMessages(source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R"))
	setwd('/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/TGGCGCCA/PCC_7716/')
	Create_Reconstruction_Files(SitesTable = "/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/TGGCGCCA/PCC_7716/Orthologues_Palindrome_sites.AllFrames.FIRST.txt",
		                    EvolutionModel = "F81",
		                    Method ="bayes",
		                    Phylogeny = "/home/lalibelulalo/TESIS/Clados/Calothrix_B/SpeciesTree_rooted.txt",
		                    OrthoPath = "/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/TGGCGCCA/PCC_7716/Only_ORTHOLOGUES/",
		                    TAXON = 'PCC_7716')
