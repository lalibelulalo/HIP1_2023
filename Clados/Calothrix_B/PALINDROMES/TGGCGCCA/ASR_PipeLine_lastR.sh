#!/bin/bash
#!/bin/bash

# USE:
# for f in $(cat Spps.txt); do sh ASR_PipeLine_lastR.sh $f GCGATCGC /home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/gbff_homologues/CalothrixspNIES-3974NIES-3974_f0_0taxa_algOMCL_e0_ Calothrix_B 7 8;done
#for f in $(cat Spps.txt); do sh ASR_PipeLine_lastR.sh $f TGGCGCCA /home/lalibelulalo/HIP1_2023/Clados/Calothrix_B/gbff_homologues/CalothrixspNIES-3974NIES-3974_f0_0taxa_algOMCL_e0_ Calothrix_B 7 8;done

SppA="$1"
PALINDROME="$2"
Path="$3"
Clade="$4"
TaxonNumber=$5
ROOT=$6
Phylogeny=$7

echo "\n"
echo "-------------------------------------"
echo "Reconstrucción Ancestral de Sitios.\n"
echo "PALINDROMO: $PALINDROME"
echo "ESPECIE:    $SppA"
echo "-------------------------------------\n"

startA=$(date +%s.%N)
SppA2=$(echo "$SppA" | sed -r 's/\//-/g')
SppB=$(echo "$SppA" | sed -r 's/ /_/g' | sed 's/\.//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/\//-/g' | sed 's/Calothrix_sp_//g' | sed 's/Calothrix_parasitica_//g')

if [ -z "$(ls -A $SppB/Only_ORTHOLOGUES)" ]; then
	echo "\nNo hay Ortólogos. Fin del codigo.                                           "
	end=$(date +%s.%N)
	runtime=$( echo "($end - $startA)/60" | bc -l )
	echo "Execution Time: $runtime mins."
else
	echo -n "\rAnálisis de mutaciones por codon y marco de lectura *"
	#---------------------------------------------------
	startB=$(date +%s.%N)
	cat <<EOF >$SppB/ReconstructionR_B.tmp.R
	suppressMessages(source("/home/lalibelulalo/HIP1_2023/ASR_Orth_Functions/NodeAndEdges.R"))
	setwd('/home/lalibelulalo/HIP1_2023/Clados/$Clade/PALINDROMES/$PALINDROME/$SppB/')

	RFS <- sort(unique(system('awk \'{if(NR!=1) {print \$5}}\' Orthologues_Palindrome_sites.AllFrames.SECOND.txt |uniq',intern = TRUE)))

	for(RF in RFS){ ## Para los tres marcos de lectura (RF)
	  CreateCodonMutationsTable(Tree = '../../../SpeciesTree_rooted.txt',
		                    Second = 'Orthologues_Palindrome_sites.AllFrames.SECOND.txt',
		                    TaxonNumber = $TaxonNumber,
		                    Root = $ROOT,
		                    RF = RF)
	  system(paste0('python3 ../../../../CodonMutations.py RF',RF,'.RECONSTRUCTION.rooted.parents.txt ',RF,' codon_mutations_RF',RF,'.rooted.txt $PALINDROME'))
	}

	#system('cat *.rooted.txt >codon_mutations_AllF.rooted.txt')
EOF
	Rscript $SppB/ReconstructionR_B.tmp.R > Analisis_Mutaciones_$SppB\_log.tmp
	end=$(date +%s.%N)
	runtime=$( echo "$end - $startB" | bc -l )
	#---------------------------------------------------
	echo "\rAnálisis de mutaciones por codon y marco de lectura ✓	Runtime: $runtime s.\n"
	end=$(date +%s.%N)
	runtime=$( echo "($end - $startA)/60" | bc -l )
	echo "Execution Time: $runtime mins."
	rm *.tmp
fi

