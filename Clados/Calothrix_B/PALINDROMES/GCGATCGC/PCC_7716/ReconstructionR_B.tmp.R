	suppressMessages(source("/home/lalibelulalo/HIP1_2023/ASR_Orth_Functions/NodeAndEdges.R"))
	setwd('/home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/PALINDROMES/GCGATCGC/PCC_7716/')

	RFS <- sort(unique(system('awk \'{if(NR!=1) {print $5}}\' Orthologues_Palindrome_sites.AllFrames.SECOND.txt |uniq',intern = TRUE)))

	for(RF in RFS){ ## Para los tres marcos de lectura (RF)
	  CreateCodonMutationsTable(Tree = '../../../SpeciesTree_rooted.txt',
		                    Second = 'Orthologues_Palindrome_sites.AllFrames.SECOND.txt',
		                    TaxonNumber = 7,
		                    Root = 8,
		                    RF = RF)
	  system(paste0('python3 ../../../../CodonMutations.py RF',RF,'.RECONSTRUCTION.rooted.parents.txt ',RF,' codon_mutations_RF',RF,'.rooted.txt GCGATCGC'))
	}

	#system('cat *.rooted.txt >codon_mutations_AllF.rooted.txt')
