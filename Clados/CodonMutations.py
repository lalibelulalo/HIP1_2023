import re
import sys
import blosum as bl

FILE = sys.argv[1] #'RF3.RECONSTRUCTION.rooted.parents.txt' #
ML = sys.argv[2] #'3' #
output_file = sys.argv[3] #'codon_mutations_RF3.rooted.txt' #
PALINDROME = sys.argv[4]

RF = open(FILE, 'r')
mat = bl.BLOSUM(62)
DEL = mat["-"]["-"]

if ML == '1':
    A,B,C, = range(0,3),range(3,6),range(6,8)
    SA,EA,SB,EB,SC,EC = 0,3,3,6,6,9
elif ML == '2':
    A,B,C=range(1,3),range(3,6),range(6,9)
    SA,EA,SB,EB,SC,EC = 0,3,3,6,6,9
elif ML == '3':
    A,B,C,D=range(2,3),range(3,6),range(6,9),range(9,10)
    SA,EA,SB,EB,SC,EC,SD,ED = 0,3,3,6,6,9,9,12

output = open (output_file, 'w') ## Abrimos el archivo de salida
if ML == '3':
    output.write('FILE\tSpp\tSTART\tEND\tReadingFrame\tOrthLength\tSequence\tAA\tparent\tParentSequence\tparentAA\tC1Mutations\tC2Mutations\tC3Mutations\tC4Mutations\tTotalMutations\tBL1\tBL2\tBL3\tBL4\tAASubs\tSType\tAncestorSite\tActualSite\tAncestorType\tActualType\tHIPMutPos1\tHIPMutPos2\tHIPMutPos3\tHIPMutPos4\tHIPMutPos5\tHIPMutPos6\tHIPMutPos7\tHIPMutPos8\tHIPMutTypePos1\tHIPMutTypePos2\tHIPMutTypePos3\tHIPMutTypePos4\tHIPMutTypePos5\tHIPMutTypePos6\tHIPMutTypePos7\tHIPMutTypePos8\tTotalTransitions\tTotalTransversions\tTotalNoMutations\tTotalInsertions\tTotalDeletions\n')
else:
    output.write('FILE\tSpp\tSTART\tEND\tReadingFrame\tOrthLength\tSequence\tAA\tparent\tParentSequence\tparentAA\tC1Mutations\tC2Mutations\tC3Mutations\tTotalMutations\tBL1\tBL2\tBL3\tAASubs\tSType\tAncestorSite\tActualSite\tAncestorType\tActualType\tHIPMutPos1\tHIPMutPos2\tHIPMutPos3\tHIPMutPos4\tHIPMutPos5\tHIPMutPos6\tHIPMutPos7\tHIPMutPos8\tHIPMutTypePos1\tHIPMutTypePos2\tHIPMutTypePos3\tHIPMutTypePos4\tHIPMutTypePos5\tHIPMutTypePos6\tHIPMutTypePos7\tHIPMutTypePos8\tTotalTransitions\tTotalTransversions\tTotalNoMutations\tTotalInsertions\tTotalDeletions\n')

Lines = RF.readlines()

for line in Lines:
    Mutations4,C4S,Codon4 = 0,'',''
    line = re.sub('\n', '', line)
    NodeRec = line.split('\t')
    ## CHILD & PARENTS
    From,To,AA,ParentAA = NodeRec[9],NodeRec[6],NodeRec[7],NodeRec[10]
    ## BLOSUM SCORES
    if ML == '3':
        BL1,BL2,BL3,BL4 = mat[ParentAA[0]][AA[0]],mat[ParentAA[1]][AA[1]],mat[ParentAA[2]][AA[2]],mat[ParentAA[3]][AA[3]]
    else:
        BL1,BL2,BL3 = mat[ParentAA[0]][AA[0]],mat[ParentAA[1]][AA[1]],mat[ParentAA[2]][AA[2]]
    
    HIPMut1,HIPMut2,HIPMut3,HIPMut4,HIPMut5,HIPMut6,HIPMut7,HIPMut8 = 0,0,0,0,0,0,0,0
    SMT1,SMT2,SMT3,SMT4,SMT5,SMT6,SMT7,SMT8 = 'NoMutation','NoMutation','NoMutation','NoMutation','NoMutation','NoMutation','NoMutation','NoMutation'
    ##-----------------------
    ## NUC
    Codon1,Mutations1,MutPosC1= '',0,0
    for i in A:
        MutPosC1+=1
        Codon1+=From[i]
        if From[i] != To[i]:
            Mutations1 += 1
            if ML == '1' and MutPosC1 == 1:
                HIPMut1 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT1 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT1 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT1 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT1 = 'Insertion'
                elif From[i] == To[i]:
                    SMT1 = 'NoMutation'
            elif ML == '1' and MutPosC1 == 2:
                HIPMut2 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT2 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT2 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT2 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT2 = 'Insertion'
                elif From[i] == To[i]:
                    SMT2 = 'NoMutation'
            elif ML == '1' and MutPosC1 == 3:
                HIPMut3 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT3 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT3 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT3 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT3 = 'Insertion'
                elif From[i] == To[i]:
                    SMT3 = 'NoMutation'
            #--------------------------------
            elif ML == '2' and MutPosC1 == 1:
                HIPMut1 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT1 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT1 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT1 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT1 = 'Insertion'
                elif From[i] == To[i]:
                    SMT1 = 'NoMutation'
            elif ML == '2' and MutPosC1 == 2:
                HIPMut2 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT2 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT2 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT2 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT2 = 'Insertion'
                elif From[i] == To[i]:
                    SMT2 = 'NoMutation'
            #--------------------------------
            elif ML == '3' and MutPosC1 == 1:
                HIPMut1 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT1 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT1 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT1 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT1 = 'Insertion'
                elif From[i] == To[i]:
                    SMT1 = 'NoMutation'
                

    ## AA
    # No hay Cambio de AA
    if ParentAA[0] == AA[0]:
        if ParentAA[0] != '-' or ParentAA[0] != '_':
            C1S = 'S'
    # Hay cambio de AA pero el cambio es conservativo 
    elif BL1>=0:
        C1S = 's'
    # Hay una delecion
    elif BL1 == DEL:
        #FROM[SA:EA],FROM[SB:EB],FROM[SC:EC],FROM[SD:ED]
        if ParentAA[0] == '_' and AA[0] == '_':
            if (From[SA:EA]==To[SA:EA]):
                C1S = 'S'
            if (From[SA:EA]=='TAA' and To[SA:EA]=='TGA') or (From[SA:EA]=='TAA' and To[SA:EA]=='TAG'):
                C1S = 'S'
            if (From[SA:EA]=='TGA' and To[SA:EA]=='TAA') or (From[SA:EA]=='TGA' and To[SA:EA]=='TAG'):
                C1S = 'S'
            if (From[SA:EA]=='TAG' and To[SA:EA]=='TAA') or (From[SA:EA]=='TAG' and To[SA:EA]=='TGA'):
                C1S = 'S'
        if (ParentAA[0] != '_' and AA[0] == '_') or (ParentAA[0] == '_' and AA[0] != '_'):
            C1S = 'N'
        if ParentAA[0] != '-' and AA[0] == '-':
            C1S = 'D'
        if ParentAA[0] == '-' and AA[0] != '-':
            C1S = 'I'
    # Hay un cambio de AA y no es conservativo 
    elif BL1 < 0:
        C1S = 'N'
    ##-----------------------
    ## NUC
    Codon2,Mutations2,MutPosC2 = '',0,0
    for i in B:
        MutPosC2+=1
        Codon2+=From[i]
        if From[i] != To[i]:
            Mutations2 += 1
            if ML == '1' and MutPosC2 == 1:
                HIPMut4 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT4 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT4 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT4 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT4 = 'Insertion'
                elif From[i] == To[i]:
                    SMT4 = 'NoMutation'
            elif ML == '1' and MutPosC2 == 2:
                HIPMut5 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT5 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT5 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT5 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT5 = 'Insertion'
                elif From[i] == To[i]:
                    SMT5 = 'NoMutation'
            elif ML == '1' and MutPosC2 == 3:
                HIPMut6 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT6 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT6 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT6 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT6 = 'Insertion'
                elif From[i] == To[i]:
                    SMT6 = 'NoMutation'
            #--------------------------------
            elif ML == '2' and MutPosC2 == 1:
                HIPMut3 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT3 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT3 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT3 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT3 = 'Insertion'
                elif From[i] == To[i]:
                    SMT3 = 'NoMutation'
            elif ML == '2' and MutPosC2 == 2:
                HIPMut4 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT4 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT4 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT4 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT4 = 'Insertion'
                elif From[i] == To[i]:
                    SMT4 = 'NoMutation'
            elif ML == '2' and MutPosC2 == 3:
                HIPMut5 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT5 = 'Transition'
                if (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT5 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT5 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT5 = 'Insertion'
                elif From[i] == To[i]:
                    SMT5 = 'NoMutation'
            #--------------------------------
            elif ML == '3' and MutPosC2 == 1:
                HIPMut2 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT2 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT2 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT2 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT2 = 'Insertion'
                elif From[i] == To[i]:
                    SMT2 = 'NoMutation'
            elif ML == '3' and MutPosC2 == 2:
                HIPMut3 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT3 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT3 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT3 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT3 = 'Insertion'
                elif From[i] == To[i]:
                    SMT3 = 'NoMutation'
            elif ML == '3' and MutPosC2 == 3:
                HIPMut4 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT4 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT4 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT4 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT4 = 'Insertion'
                elif From[i] == To[i]:
                    SMT4 = 'NoMutation'
   
    ## AA
    if ParentAA[1] == AA[1]:
        if ParentAA[1] != '-' or ParentAA[1] != '_':
            C2S = 'S'
    elif BL2>=0:
        C2S = 's'
    elif BL2 == DEL:
        if ParentAA[1] == '_' and AA[1] == '_':
            if (From[SB:EB]==To[SA:EA]):
                C2S = 'S'
            if (From[SB:EB]=='TAA' and To[SB:EB]=='TGA') or (From[SB:EB]=='TAA' and To[SB:EB]=='TAG'):
                C2S = 'S'
            if (From[SB:EB]=='TGA' and To[SB:EB]=='TAA') or (From[SB:EB]=='TGA' and To[SB:EB]=='TAG'):
                C2S = 'S'
            if (From[SB:EB]=='TAG' and To[SB:EB]=='TAA') or (From[SB:EB]=='TAG' and To[SB:EB]=='TGA'):
                C2S = 'S'
        if (ParentAA[1] != '_' and AA[1] == '_') or (ParentAA[1] == '_' and AA[1] != '_'):
            C2S = 'N'
        if ParentAA[1] != '-' and AA[1] == '-':
            C2S = 'D'
        if ParentAA[1] == '-' and AA[1] != '-':
            C2S = 'I'
    elif BL2 < 0:
        C2S = 'N'
    ##-----------------------
    ## NUC
    Codon3,Mutations3,MutPosC3 = '',0,0
    for i in C:
        MutPosC3+=1
        Codon3+=From[i]
        if From[i] != To[i]:
            Mutations3 += 1
            if ML == '1' and MutPosC3 == 1:
                HIPMut7 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT7 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT7 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT7 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT7 = 'Insertion'
                elif From[i] == To[i]:
                    SMT7 = 'NoMutation'
            elif ML == '1' and MutPosC3 == 2:
                HIPMut8 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT8 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT8 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT8 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT8 = 'Insertion'
                elif From[i] == To[i]:
                    SMT8 = 'NoMutation'
            #--------------------------------
            elif ML == '2' and MutPosC3 == 1:
                HIPMut6 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT6 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT6 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT6 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT6 = 'Insertion'
                elif From[i] == To[i]:
                    SMT6 = 'NoMutation'
            elif ML == '2' and MutPosC3 == 2:
                HIPMut7 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT7 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT7 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT7 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT7 = 'Insertion'
                elif From[i] == To[i]:
                    SMT7 = 'NoMutation'
            elif ML == '2' and MutPosC3 == 3:
                HIPMut8 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT8 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT8 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT8 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT8 = 'Insertion'
                elif From[i] == To[i]:
                    SMT8 = 'NoMutation'
            #--------------------------------
            elif ML == '3' and MutPosC3 == 1:
                HIPMut5 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT5 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT5 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT5 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT5 = 'Insertion'
                elif From[i] == To[i]:
                    SMT5 = 'NoMutation'
            elif ML == '3' and MutPosC3 == 2:
                HIPMut6 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT6 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT6 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT6 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT6 = 'Insertion'
                elif From[i] == To[i]:
                    SMT6 = 'NoMutation'
            elif ML == '3' and MutPosC3 == 3:
                HIPMut7 = 1
                if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                    SMT7 = 'Transition'
                elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                    SMT7 = 'Transversion'
                elif From[i] != '-' and To[i] == '-':
                    SMT7 = 'Deletion'
                elif From[i] == '-' and To[i] != '-':
                    SMT7 = 'Insertion'
                elif From[i] == To[i]:
                    SMT7 = 'NoMutation'
                    
    ## AA   
    if ParentAA[2] == AA[2]:
        if ParentAA[2] != '-' or ParentAA[2] != '_':
            C3S = 'S'
    elif BL3 >= 0:
        C3S = 's'
    elif BL3 == DEL:
        if ParentAA[2] == '_' and AA[2] == '_':
            if (From[SC:EC]==To[SC:EC]):
                C3S = 'S'
            if (From[SC:EC]=='TAA' and To[SC:EC]=='TGA') or (From[SC:EC]=='TAA' and To[SC:EC]=='TAG'):
                C3S = 'S'
            if (From[SC:EC]=='TGA' and To[SC:EC]=='TAA') or (From[SC:EC]=='TGA' and To[SC:EC]=='TAG'):
                C3S = 'S'
            if (From[SC:EC]=='TAG' and To[SC:EC]=='TAA') or (From[SC:EC]=='TAG' and To[SC:EC]=='TGA'):
                C3S = 'S'
        if (ParentAA[2] != '_' and AA[2] == '_') or (ParentAA[2] == '_' and AA[2] != '_'):
            C3S = 'N'
        if ParentAA[2] != '-' and AA[2] == '-':
            C3S = 'D'
        if ParentAA[2] == '-' and AA[2] != '-':
            C3S = 'I'
    elif BL3 < 0:
        C3S = 'N'
    ##-----------------------
    ## NUC
    if ML == '3':
        Codon4,Mutations4,MutPosC3 = '',0,0
        for i in D:
            MutPosC3+=1
            Codon4+=From[i]
            if From[i] != To[i]:
                Mutations4 += 1
                if ML == '3' and MutPosC3 == 1:
                    HIPMut8 = 1
                    if (From[i] == 'A' and To[i] == 'G') or (From[i] == 'G' and To[i] == 'A') or (From[i] == 'T' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'T'):
                        SMT8 = 'Transition'
                    elif (From[i] == 'A' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'G') and (From[i] == 'A' and To[i] == 'C') or (From[i] == 'C' and To[i] == 'A') or (From[i] == 'G' and To[i] == 'T') or (From[i] == 'T' and To[i] == 'G'):
                        SMT8 = 'Transversion'
                    elif From[i] != '-' and To[i] == '-':
                        SMT8 = 'Deletion'
                    elif From[i] == '-' and To[i] != '-':
                        SMT8 = 'Insertion'
                    elif From[i] == To[i]:
                        SMT8 = 'NoMutation'

        ## AA
        if ParentAA[3] == AA[3]:
            if ParentAA[3] != '-' or ParentAA[3] != '_':
                C4S = 'S'
        elif BL4 >= 0:
            C4S = 's'
        elif BL4 == DEL:
            if ParentAA[3] == '_' and AA[3] == '_':
                if (From[SD:ED]==To[SD:ED]):
                    C4S = 'S'
                if (From[SD:ED]=='TAA' and To[SD:ED]=='TGA') or (From[SD:ED]=='TAA' and To[SD:ED]=='TAG'):
                    C4S = 'S'
                if (From[SD:ED]=='TGA' and To[SD:ED]=='TAA') or (From[SD:ED]=='TGA' and To[SD:ED]=='TAG'):
                    C4S = 'S'
                if (From[SD:ED]=='TAG' and To[SD:ED]=='TAA') or (From[SD:ED]=='TAG' and To[SD:ED]=='TGA'):
                    C4S = 'S'
            if (ParentAA[3] != '_' and AA[3] == '_') or (ParentAA[3] == '_' and AA[3] != '_'):
                C4S = 'N'
            if ParentAA[3] != '-' and AA[3] == '-':
                C4S = 'D'
            if ParentAA[3] == '-' and AA[3] != '-':
                C4S = 'I'
        elif BL4 < 0:
            C4S = 'N'
        #TotalMutations = Mutations1+Mutations2+Mutations3+Mutations4
        #AASubs = C1S+C2S+C3S+C4S
    ##-----------------------
    if ML == '1':
        CHILD = To[0:8]
    if ML == '2':
        CHILD = To[1:9]
    if ML == '3':
        CHILD = To[2:10]
    
    ##-----------------------
    SMTs,TransitionCount,TransversionCount,NoMutationCount,InsertionCount,DeletionCount = [SMT1,SMT2,SMT3,SMT4,SMT5,SMT6,SMT7,SMT8],0,0,0,0,0
    for smt in SMTs:
        if smt == 'Transition':
            TransitionCount+=1
        elif smt == 'Transversion':
            TransversionCount+=1
        elif smt == 'Deletion':
            DeletionCount+=1
        elif smt == 'Insertion':
            InsertionCount+=1
        elif smt == 'NoMutation':
            NoMutationCount+=1
        
    ##-----------------------

    TotalMutations = Mutations1+Mutations2+Mutations3+Mutations4
    AASubs = C1S+C2S+C3S+C4S
    PalSite = Codon1+Codon2+Codon3+Codon4
    if PalSite == PALINDROME:
        Sequence = 'SITE'
    else:
        Sequence = 'NoSITE'
    ##-----------------------
    if CHILD == PALINDROME:
        Sequence2 = 'SITE'
    else:
        Sequence2 = 'NoSITE'
    ##-----------------------
            
    if ML!='3':
        if AASubs == 'SSS':
            if TotalMutations == 0:
                SType = 'NoMutation'
            else:
                SType = 'Synonym'
        elif len(re.findall("[Ss][Ss][Ss]", AASubs))==1:
            if TotalMutations == 0:
                SType = 'ConservativeNoSiteMut'
            else:
                SType = 'Conservative'
        elif (len(re.findall("[DSs][DSs][DSs]", AASubs))==1) or (len(re.findall("[DN][DN][DN]", AASubs))==1):
            SType = 'Deletion'
        elif (len(re.findall("D[NSs][NSs]", AASubs))==1) or (len(re.findall("[NSs]D[NSs]", AASubs))==1) or (len(re.findall("[NSs][NSs]D", AASubs))==1):
            SType = 'Deletion'
        elif (len(re.findall("I[DNSs][DNSs]", AASubs))==1) or (len(re.findall("[DNSs]I[DNSs]", AASubs))==1) or (len(re.findall("[DNSs][DNSs]I", AASubs))==1):
            SType = 'Deletion'
        elif len(re.findall("[NSs][NSs][NSs]", AASubs))==1:
            if TotalMutations == 0:
                SType = 'NoSynonymNoSiteMut'
            else:
                SType = 'NoSynonym'
    elif ML=='3':
        if AASubs == 'SSSS':
            if TotalMutations == 0:
                SType = 'NoMutation'
            else:
                SType = 'Synonym'
        elif len(re.findall("[Ss][Ss][Ss][Ss]", AASubs))==1:
            if TotalMutations == 0:
                SType = 'ConservativeNoSiteMut'
            else:
                SType = 'Conservative'
        elif len(re.findall("[DSs][DSs][DSs][DSs]", AASubs))==1:
            SType = 'Deletion'         
        elif len(re.findall("[NSs][NSs][NSs][NSs]", AASubs))==1:
            if TotalMutations == 0:
                SType = 'NoSynonymNoSiteMut'
            else:
                SType = 'NoSynonym'
    
    if ML == '3':
            #print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(line,Mutations1,Mutations2,Mutations3,Mutations4,TotalMutations,BL1,BL2,BL3,BL4,C1S,C2S,C3S,AASubs,SType,PalSite,CHILD,Sequence,HIPMut1,HIPMut2,HIPMut3,HIPMut4,HIPMut5,HIPMut6,HIPMut7,HIPMut8,SMT1,SMT2,SMT3,SMT4,SMT5,SMT6,SMT7,SMT8,TransitionCount,TransversionCount))
            output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(line,Mutations1,Mutations2,Mutations3,Mutations4,TotalMutations,BL1,BL2,BL3,BL4,AASubs,SType,PalSite,CHILD,Sequence,Sequence2,HIPMut1,HIPMut2,HIPMut3,HIPMut4,HIPMut5,HIPMut6,HIPMut7,HIPMut8,SMT1,SMT2,SMT3,SMT4,SMT5,SMT6,SMT7,SMT8,TransitionCount,TransversionCount,NoMutationCount,InsertionCount,DeletionCount))
    else:
        #print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(line,Mutations1,Mutations2,Mutations3,TotalMutations,BL1,BL2,BL3,C1S,C2S,C3S,AASubs,SType,PalSite,CHILD,Sequence,HIPMut1,HIPMut2,HIPMut3,HIPMut4,HIPMut5,HIPMut6,HIPMut7,HIPMut8,SMT1,SMT2,SMT3,SMT4,SMT5,SMT6,SMT7,SMT8,TransitionCount,TransversionCount))
        output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(line,Mutations1,Mutations2,Mutations3,TotalMutations,BL1,BL2,BL3,AASubs,SType,PalSite,CHILD,Sequence,Sequence2,HIPMut1,HIPMut2,HIPMut3,HIPMut4,HIPMut5,HIPMut6,HIPMut7,HIPMut8,SMT1,SMT2,SMT3,SMT4,SMT5,SMT6,SMT7,SMT8,TransitionCount,TransversionCount,NoMutationCount,InsertionCount,DeletionCount))

output.close()
