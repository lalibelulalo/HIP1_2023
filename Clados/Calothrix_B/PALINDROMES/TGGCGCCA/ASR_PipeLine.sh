#!/bin/bash
#!/bin/bash

# USE:
# for f in $(cat Spps.txt); do sh ASR_PipeLine.sh $f TGGCGCCA /home/lalibelulalo/TESIS/Clados/Calothrix_B/gbff_homologues/CalothrixspNIES-3974NIES-3974_f0_0taxa_algOMCL_e0_ Calothrix_B 7 8;done

#for f in $(cat Spps.txt); do sh ASR_PipeLine.sh $f GCGATCGC /home/lalibelulalo/TESIS/Clados/Calothrix_B/gbff_homologues/CalothrixspNIES-3974NIES-3974_f0_0taxa_algOMCL_e0_ Calothrix_B 7 8;done

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
echo -n "\rObteniendo ortólogos con sitios $PALINDROME            *"
#---------------------------------------------------
startB=$(date +%s.%N)
mkdir $SppB

awk -v r=$1 '{if ($5!=0 && $2==r){print $1}}' Markov_count_$PALINDROME\_Octanuc_.txt | sed 's/.fna//g' | uniq >$SppB/Ortologos_$PALINDROME\_$SppB.txt

mkdir $SppB/Orthologues_$PALINDROME\_$SppB

for word in $(cat $SppB/Ortologos_$PALINDROME\_$SppB.txt); do cp $Path/$word.fna $SppB/Orthologues_$PALINDROME\_$SppB; done
for word in $(cat $SppB/Ortologos_$PALINDROME\_$SppB.txt); do cp $Path/$word.faa $SppB/Orthologues_$PALINDROME\_$SppB; done

sleep 1
end=$(date +%s.%N)
runtime=$( echo "$end - $startB" | bc -l )
#---------------------------------------------------
echo "\rObteniendo ortólogos con sitios $PALINDROME            ✓	Runtime: $runtime s."



echo -n "\rFiltrando Parálogos                                 *"
#---------------------------------------------------
startB=$(date +%s.%N)
mkdir $SppB/Only_ORTHOLOGUES
mkdir $SppB/PARALOGUES

python3 ../../../FiltradoParalogosSHELL.py $SppB/Orthologues_$PALINDROME\_$SppB $TaxonNumber $SppB/Only_ORTHOLOGUES/ $SppB/PARALOGUES/

if [ -z "$(ls -A $SppB/Only_ORTHOLOGUES)" ]; then
	echo "\rFiltrando Parálogos                                 ✓	Runtime: $runtime s."
	echo "\nNo hay Ortólogos. Fin del codigo.                                           "
	end=$(date +%s.%N)
	runtime=$( echo "($end - $startA)/60" | bc -l )
	echo "Execution Time: $runtime mins."
else
	for f in $SppB/Only_ORTHOLOGUES/*.fna; do awk -F "|" '{if (NR%2){print ">"$2}else{print $1}}' $f | sed 's/ /_/g' | sed 's/\.//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/\//-/g' | sed 's/Calothrix_sp_//g' | sed 's/Calothrix_parasitica_//g' >$f.awk1;done
	for f in $SppB/Only_ORTHOLOGUES/*.faa; do awk -F "|" '{if (NR%2){print ">"$2}else{print $1}}' $f | sed 's/ /_/g' | sed 's/\.//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/\//-/g' | sed 's/Calothrix_sp_//g' | sed 's/Calothrix_parasitica_//g' >$f.awk1;done
	sleep 1
	end=$(date +%s.%N)
	runtime=$( echo "$end - $startB" | bc -l )
	#---------------------------------------------------
	echo "\rFiltrando Parálogos                                 ✓	Runtime: $runtime s."
	
	echo -n "\rAlineación multiple                                 *"
	#---------------------------------------------------
	startB=$(date +%s.%N)
	#sh Alineacion.sh $SppB > Alineacion_log.tmp > /dev/null 2>&1
	for f in $SppB/Only_ORTHOLOGUES/*faa.awk1; do mafft $f >$f.mafft;done > /dev/null 2>&1
	ls $SppB/Only_ORTHOLOGUES/*.fna | sed 's/.fna//g' >$SppB/only.orthologues.txt
	for f in $(cat $SppB/only.orthologues.txt); do perl ../../../pal2nal.pl $f.faa.awk1.mafft $f.fna.awk1 -output fasta >$f.codon.alg.fasta;done > /dev/null 2>&1
	python3 ../../../Fasta2Phylip.py $SppB/Only_ORTHOLOGUES/
	end=$(date +%s.%N)
	runtime=$( echo "$end - $startB" | bc -l )
	#---------------------------------------------------
	echo "\rAlineación multiple                                 ✓	Runtime: $runtime s."
	
	echo -n "\rObteniendo coordenadas de los sitios $PALINDROME       *"
	#---------------------------------------------------
	startB=$(date +%s.%N)
	python3 ../../../AlignmentPalindromeCoords.py $SppB/Only_ORTHOLOGUES/ $PALINDROME $SppB
	mv Orthologues_Palindrome_sites.txt $SppB/
	python3 ../../../AlignmentPalindromeCoordsAA.py $SppB/Only_ORTHOLOGUES/ $PALINDROME $SppB FIRST >FIRST.output
	mv Orthologues_Palindrome_sites.AllFrames.FIRST.txt $SppB/
	mv FIRST.output $SppB/
	mv $SppB\_Reading_Frames.counts $SppB/
	mv CodonErrors.$SppB.txt $SppB/
	sleep 1
	end=$(date +%s.%N)
	runtime=$( echo "$end - $startB" | bc -l )
	#---------------------------------------------------
	echo "\rObteniendo coordenadas de los sitios $PALINDROME       ✓	Runtime: $runtime s."
	
	echo -n "\rGrafo                                               *"
	#---------------------------------------------------
	startB=$(date +%s.%N)
	cat <<EOF >$SppB/ReconstructionR_C.tmp.R
	suppressMessages(source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R"))
	setwd('/home/lalibelulalo/TESIS/Clados/$Clade/PALINDROMES/$PALINDROME/$SppB/')
	Nodes.Edges <- Create_Transition_Table_No_Fit(SitesTable = "/home/lalibelulalo/TESIS/Clados/$Clade/PALINDROMES/$PALINDROME/$SppB/Orthologues_Palindrome_sites.txt",
		                        EvolutionModel = "F81",
		                        Method = "bayes",
		                        Phylogeny = "/home/lalibelulalo/TESIS/Clados/$Clade/SpeciesTree_rooted.txt",
		                        OrthoPath = "/home/lalibelulalo/TESIS/Clados/$Clade/PALINDROMES/$PALINDROME/$SppB/Only_ORTHOLOGUES/",
		                        OutName = "All_RF")
EOF
	Rscript $SppB/ReconstructionR_C.tmp.R > Grafo_$SppB\_log.tmp
	end=$(date +%s.%N)
	runtime=$( echo "$end - $startB" | bc -l )
	#---------------------------------------------------
	echo "\rGrafo                                               ✓	Runtime: $runtime s."
	
	echo -n "\rReconstrucción Ancestral de sitios $PALINDROME         *"
	#---------------------------------------------------
	startB=$(date +%s.%N)
	cat <<EOF >$SppB/ReconstructionR_A.tmp.R
	suppressMessages(source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R"))
	setwd('/home/lalibelulalo/TESIS/Clados/$Clade/PALINDROMES/$PALINDROME/$SppB/')
	Create_Reconstruction_Files(SitesTable = "/home/lalibelulalo/TESIS/Clados/$Clade/PALINDROMES/$PALINDROME/$SppB/Orthologues_Palindrome_sites.AllFrames.FIRST.txt",
		                    EvolutionModel = "F81",
		                    Method ="bayes",
		                    Phylogeny = "/home/lalibelulalo/TESIS/Clados/$Clade/SpeciesTree_rooted.txt",
		                    OrthoPath = "/home/lalibelulalo/TESIS/Clados/$Clade/PALINDROMES/$PALINDROME/$SppB/Only_ORTHOLOGUES/",
		                    TAXON = '$SppB')
EOF
	Rscript $SppB/ReconstructionR_A.tmp.R > Reconstrucción_Ancestral_$SppB\_log.tmp
	mkdir $SppB/RECONSTRUCCIONES
	mv $SppB/Only_ORTHOLOGUES/*.RECONSTRUCTION.fasta $SppB/RECONSTRUCCIONES
	python3 ../../../Fasta2Phylip.py $SppB/RECONSTRUCCIONES/
	python3 ../../../AlignmentPalindromeCoordsAA.py $SppB/RECONSTRUCCIONES/ $PALINDROME $SppB SECOND >SECOND.output
	mv Orthologues_Palindrome_sites.AllFrames.SECOND.txt $SppB/
	mv SECOND.output $SppB/
	mv $SppB\_Reading_Frames.counts $SppB/
	mv CodonErrors.$SppB.txt $SppB/
	end=$(date +%s.%N)
	runtime=$( echo "$end - $startB" | bc -l )
	#---------------------------------------------------
	echo "\rReconstrucción Ancestral de sitios $PALINDROME         ✓	Runtime: $runtime s."
	
	echo -n "\rAnálisis de mutaciones por codon y marco de lectura *"
	#---------------------------------------------------
	startB=$(date +%s.%N)
	cat <<EOF >$SppB/ReconstructionR_B.tmp.R
	suppressMessages(source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R"))
	setwd('/home/lalibelulalo/TESIS/Clados/$Clade/PALINDROMES/$PALINDROME/$SppB/')

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

