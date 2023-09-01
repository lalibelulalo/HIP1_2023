	suppressMessages(source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R"))
	setwd('/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/TGGCGCCA/PCC_6303/')
	Create_Reconstruction_Files(SitesTable = "/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/TGGCGCCA/PCC_6303/Orthologues_Palindrome_sites.AllFrames.FIRST.txt",
		                    EvolutionModel = "F81",
		                    Method ="bayes",
		                    Phylogeny = "/home/lalibelulalo/TESIS/Clados/Calothrix_B/SpeciesTree_rooted.txt",
		                    OrthoPath = "/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/TGGCGCCA/PCC_6303/Only_ORTHOLOGUES/",
		                    TAXON = 'PCC_6303')
