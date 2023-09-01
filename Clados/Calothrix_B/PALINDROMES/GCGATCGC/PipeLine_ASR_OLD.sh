#!/bin/bash

# USE:
# sh PipeLine_ASR.sh Calothrix_sp_336/3 GCGATCGC /home/lalibelulalo/TESIS/Clados/Calothrix_B/gbff_homologues/CalothrixspNIES-3974NIES-3974_f0_0taxa_algOMCL_e0_

SppA="$1"
PALINDROME="$2"
Path="$3"
#Path="/home/lalibelulalo/TESIS/Clados/Calothrix_B/gbff_homologues/CalothrixspNIES-3974NIES-3974_f0_0taxa_algOMCL_e0_"
Clade="$4"
TaxonNumber=$5

SppA2=$(echo "$SppA" | sed -r 's/\//-/g')
SppB=$(echo "$SppA" | sed -r 's/ /_/g' | sed 's/\.//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/\//-/g' | sed 's/Calothrix_sp_//g' | sed 's/Calothrix_parasitica_//g')

## Obtengo los ortólogos con conteos positivos del palíndromo de interés
mkdir $SppB

awk -v r=$1 '{if ($5!=0 && $2==r){print $1}}' Markov_count_$PALINDROME\_Octanuc_.txt | sed 's/.fna//g' | uniq >$SppB/Ortologos_$PALINDROME\_$SppB.txt

mkdir $SppB/Orthologues_$PALINDROME\_$SppB

for word in $(cat $SppB/Ortologos_$PALINDROME\_$SppB.txt); do cp $Path/$word.fna $SppB/Orthologues_$PALINDROME\_$SppB; done
for word in $(cat $SppB/Ortologos_$PALINDROME\_$SppB.txt); do cp $Path/$word.faa $SppB/Orthologues_$PALINDROME\_$SppB; done


## Filtro Ortólogos y Paralogos
mkdir $SppB/Only_ORTHOLOGUES
mkdir $SppB/PARALOGUES

python3 ../../../FiltradoParalogosSHELL.py $SppB/Orthologues_$PALINDROME\_$SppB $TaxonNumber $SppB/Only_ORTHOLOGUES/ $SppB/PARALOGUES/

for f in $SppB/Only_ORTHOLOGUES/*.fna; do awk -F "|" '{if (NR%2){print ">"$2}else{print $1}}' $f | sed 's/ /_/g' | sed 's/\.//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/\//-/g' | sed 's/Calothrix_sp_//g' | sed 's/Calothrix_parasitica_//g' >$f.awk1;done

for f in $SppB/Only_ORTHOLOGUES/*.faa; do awk -F "|" '{if (NR%2){print ">"$2}else{print $1}}' $f | sed 's/ /_/g' | sed 's/\.//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/\//-/g' | sed 's/Calothrix_sp_//g' | sed 's/Calothrix_parasitica_//g' >$f.awk1;done


## Alineación multiple de ortólogos
for f in $SppB/Only_ORTHOLOGUES/*faa.awk1; do mafft $f >$f.mafft;done

ls $SppB/Only_ORTHOLOGUES/*.fna | sed 's/.fna//g' >$SppB/only.orthologues.txt

for f in $(cat $SppB/only.orthologues.txt); do perl ../../../pal2nal.pl $f.faa.awk1.mafft $f.fna.awk1 -output fasta >$f.codon.alg.fasta;done

## Convierto a formato PHYLIP
python3 ../../../Fasta2Phylip.py $SppB/Only_ORTHOLOGUES/

## Obtengo las coordenadas de los palindromos
python3 ../../../AlignmentPalindromeCoords.py $SppB/Only_ORTHOLOGUES/ $PALINDROME $SppB
mv Orthologues_Palindrome_sites.txt $SppB/

## Obtengo las coordenadas de los palindromos por marco de lectura
python3 ../../../AlignmentPalindromeCoordsAA.py $SppB/Only_ORTHOLOGUES/ $PALINDROME $SppB FIRST >FIRST.output
mv Orthologues_Palindrome_sites.AllFrames.FIRST.txt $SppB/
mv FIRST.output $SppB/
mv $SppB\_Reading_Frames.counts $SppB/
mv CodonErrors.$SppB.txt $SppB/

## Hago la reconstrucción ancestral de sitios ancestrales
cat <<EOF >$SppB/ReconstructionR_A.tmp.R
source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R")
setwd('/home/lalibelulalo/TESIS/Clados/$Clade/PALINDROMES/$PALINDROME/$SppB/')
Create_Reconstruction_Files(SitesTable = "/home/lalibelulalo/TESIS/Clados/$Clade/PALINDROMES/$PALINDROME/$SppB/Orthologues_Palindrome_sites.AllFrames.FIRST.txt",
                            EvolutionModel = "F81",
                            Method ="bayes",
                            Phylogeny = "/home/lalibelulalo/TESIS/Clados/$Clade/SpeciesTree_rooted.txt",
                            OrthoPath = "/home/lalibelulalo/TESIS/Clados/$Clade/PALINDROMES/$PALINDROME/$SppB/Only_ORTHOLOGUES/",
                            TAXON = '$SppB')
EOF

Rscript $SppB/ReconstructionR_A.tmp.R

mkdir $SppB/RECONSTRUCCIONES
mv $SppB/Only_ORTHOLOGUES/*.RECONSTRUCTION.fasta $SppB/RECONSTRUCCIONES

python3 ../../../Fasta2Phylip.py $SppB/RECONSTRUCCIONES/
python3 ../../../AlignmentPalindromeCoordsAA.py $SppB/RECONSTRUCCIONES/ $PALINDROME $SppB SECOND >SECOND.output

mv Orthologues_Palindrome_sites.AllFrames.SECOND.txt $SppB/
mv SECOND.output $SppB/
mv $SppB\_Reading_Frames.counts $SppB/
mv CodonErrors.$SppB.txt $SppB/


## Análisis de mutaciones por codon y marco de lectura
cat <<EOF >$SppB/ReconstructionR_B.tmp.R
source("/home/lalibelulalo/TESIS/ASR_Orth_Functions/NodeAndEdges.R")
setwd('/home/lalibelulalo/TESIS/Clados/$Clade/PALINDROMES/$PALINDROME/$SppB/')

RFS <- length(unique(system('awk \'{if(NR!=1) {print \$5}}\' Orthologues_Palindrome_sites.AllFrames.SECOND.txt |uniq',intern = TRUE)))

for(RF in 1:RFS){ ## Para los tres marcos de lectura (RF)
  CreateCodonMutationsTable(Tree = '../../../SpeciesTree_rooted.txt',
                            Second = 'Orthologues_Palindrome_sites.AllFrames.SECOND.txt',
                            TaxonNumber = $TaxonNumber,
                            RF = RF)
  system(paste0('python3 ../../../../CodonMutations.py RF',RF,'.RECONSTRUCTION.rooted.parents.txt ',RF,' codon_mutations_RF',RF,'.rooted.txt $PALINDROME'))
}

#system('cat *.rooted.txt >codon_mutations_AllF.rooted.txt')
EOF

Rscript $SppB/ReconstructionR_B.tmp.R

