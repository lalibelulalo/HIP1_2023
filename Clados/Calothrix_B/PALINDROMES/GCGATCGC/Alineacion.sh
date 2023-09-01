#!/bin/bash
#!/bin/bash
SppB=$1

for f in $SppB/Only_ORTHOLOGUES/*faa.awk1; do mafft $f >$f.mafft;done > /dev/null

ls $SppB/Only_ORTHOLOGUES/*.fna | sed 's/.fna//g' >$SppB/only.orthologues.txt

for f in $(cat $SppB/only.orthologues.txt); do perl ../../../pal2nal.pl $f.faa.awk1.mafft $f.fna.awk1 -output fasta >$f.codon.alg.fasta;done > /dev/null

python3 ../../../Fasta2Phylip.py $SppB/Only_ORTHOLOGUES/
