from Bio import SeqIO
import re
from CountFunctions import OctanucGen
from CountFunctions import NucsGenLen
from CountFunctions import CountNucsAndKmers
from CountFunctions import PalKmerCounts
from CountFunctions import Markov
import time
import os
import sys
import numpy as np

GENOMES_PATH = sys.argv[1]
PalsFile = sys.argv[2]
GenomesFile = sys.argv[3]
ORDER = sys.argv[4]
KMERS = 8


model = str("".join(['M',str(ORDER)]))

if GENOMES_PATH.endswith("/"):
    genomepath = re.sub('/', '', GENOMES_PATH)
else:
    genomepath = GENOMES_PATH

st = time.time()
st2 = time.time()
DateTime = time.strftime("%Y,%m,%d,%H,%M,%S")
t = DateTime.split(',')
numbers = [ int(x) for x in t ]  
current_date = "-".join([str(numbers[0]),str(numbers[1]),str(numbers[2])])
current_time = "".join ([str(numbers[3]),'hrs',str(numbers[4]),'mins'])
date = str("_".join ([current_date,current_time]))

Nucs = ["N","A","G","C","T"]
with open(PalsFile) as f:
    for line in f:
        pals = [i.strip() for i in line.split(',')]
        
if KMERS == 8:
    kmer = 'Octanuc'


## Entramos al path proporcionado que contiene los genomas
if GENOMES_PATH.endswith("/"):
    GenomeDir = str(GENOMES_PATH)
else:
    GenomeDir = str("".join ([GENOMES_PATH,'/']))
    
#genomes = [x for x in os.listdir(GenomeDir) if x.endswith(".fna")]
genomes = np.loadtxt(GenomesFile,dtype="str")

contador = 0
PalCounter = 0

for pal in pals:
    PalCounter += 1
    GP = re.sub('/', '', GENOMES_PATH)
    output_file = str("_".join (['Markov_count',pal,kmer,'.txt']))#,date
    header = str("".join(['prot\tspp\tID\tpalindrome\tobs\tmarkov',str(ORDER),'\tgenomesize\tA\tTh\tC\tG\tN\n']))
    output = open (output_file, 'w')
    output.write(header)
    print ("The output file is: {}.\n".format(output_file))
    print ("Palindromo {} de {}: {}.".format(PalCounter,len(pals),pal))

    pal = pal.upper()
    for GenomeFile in genomes:
        GenomeFile2 = str("".join ([GenomeDir,GenomeFile]))
        KNUCSUM = 0
        contador += 1
        fasta_sequences = SeqIO.parse(open(GenomeFile2),'fasta')
        print ("\tArchivo {} de {}: {}".format(contador,len(genomes), GenomeFile))
        for fasta in fasta_sequences:
            name, sequence, description = fasta.id, str(fasta.seq), fasta.description
            DescriptionList = description.split("|")
            ID = re.sub('ID:', '', name)
            Spp = DescriptionList[1]
            Spp = re.sub(' ','_',Spp)
            Spp = re.sub(']','',Spp)
            Spp = re.sub('\[','',Spp)
            Spp = re.sub('\.','',Spp)
            # BUSCAMOS NUCLEOTIDOS EN TODO EL GENOMA Y OBTENEMOS EL TAMAÑO DEL GENOMA
            NUCLEOTIDES, genome_length = NucsGenLen(sequence)
            # k-meros del genoma
            NucsKmers = CountNucsAndKmers(sequence,KMERS)
            if KMERS == 8 or KMERS == 82 or KMERS == 83 or KMERS == 61:
                DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES,PENTANUCLEOTIDES,HEXANUCLEOTIDES,HEPTANUCLEOTIDES,OCTANUCLEOTIDES = PalKmerCounts(pal,KMERS,NucsKmers)
                KNUCS = OCTANUCLEOTIDES
                if KMERS == 4 or KMERS == 5 or KMERS == 6 or KMERS == 7 or KMERS == 73 or KMERS == 8 or KMERS == 82 or KMERS == 83:
                    markov = Markov(pal,genome_length,ORDER,NUCLEOTIDES,DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES)
                    output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(GenomeFile,Spp,ID,pal,KNUCS[pal],markov,genome_length,NUCLEOTIDES['A'],NUCLEOTIDES['T'],NUCLEOTIDES['C'],NUCLEOTIDES['G'],NUCLEOTIDES['N']))   
    # get the end time
    et = time.time()
    # get the execution time
    elapsed_time2 = et - st2
    st2 = time.time()
    contador = 0
    print('Execution time: {} seconds.'.format(elapsed_time2))
    print ("------------------------------------")
    output.close()
elapsed_time = et - st
print ('*** Execution time: {} mins. ***'.format(elapsed_time/60))
