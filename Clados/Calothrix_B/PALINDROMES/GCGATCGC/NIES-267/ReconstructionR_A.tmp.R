	suppressMessages(source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R"))
	setwd('/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/GCGATCGC/NIES-267/')
	Create_Reconstruction_Files(SitesTable = "/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/GCGATCGC/NIES-267/Orthologues_Palindrome_sites.AllFrames.FIRST.txt",
		                    EvolutionModel = "F81",
		                    Method ="bayes",
		                    Phylogeny = "/home/lalibelulalo/TESIS/Clados/Calothrix_B/SpeciesTree_rooted.txt",
		                    OrthoPath = "/home/lalibelulalo/TESIS/Clados/Calothrix_B/PALINDROMES/GCGATCGC/NIES-267/Only_ORTHOLOGUES/",
		                    TAXON = 'NIES-267')
