## Written by Claire Goul -- Oct 2023
## Based on algorithm to calculate FFAT Score described in Di Mattia et al, EMBO 2020
## accounts for phospho- FFATs
import pandas as pd
from Bio import SeqIO
import Bio.SeqIO.FastaIO
import sys

def getFFATs(residuesdict,fasta_sequences):
    top_FFATs_genome=[]
    with open("HumanProteome.fasta") as handle:
        for value in Bio.SeqIO.FastaIO.SimpleFastaParser(handle):
            seq=value[1]
            ID=value[0]
            if len(seq)>=19:
                for i in range(0,len(seq)-19):
                    subseq=seq[i:i+19]
                    if 'U' not in subseq:
                        tractn=residuedict[0][subseq[0]]+residuedict[1][subseq[1]]+residuedict[2][subseq[2]]+residuedict[3][subseq[3]]+residuedict[4][subseq[4]]+residuedict[5][subseq[5]]
                        score1=residuedict[6][subseq[6]]+residuedict[7][subseq[7]]+residuedict[8][subseq[8]]+residuedict[9][subseq[9]]+residuedict[10][subseq[10]]+residuedict[11][subseq[11]]+residuedict[12][subseq[12]]
                        tractc=residuedict[13][subseq[13]]+residuedict[14][subseq[14]]+residuedict[15][subseq[15]]+residuedict[16][subseq[16]]+residuedict[17][subseq[17]]+residuedict[18][subseq[18]]
                        tractcn=tractc+tractn
                        #score2
                        test=1.5
                        if tractn >=2:
                            test=1
                        if tractn >=3:
                            test1=0.5
                        elif tractn<3:
                            test1=test
                        if tractn>=4:
                            score2=0
                        else:
                            score2=test1
                        #score3
                        if tractc >=2:
                            test=1
                        elif tractc <2:
                            test = 1.5
                        if tractc >=3:
                            test1=0.5
                        elif tractc<3:
                            test1=test
                        if tractc>=4:
                            score3=0
                        else:
                            score3=test1
                        #score4  
                        if tractcn >=3:
                            test=1
                        elif tractcn <3:
                            test = 1.5
                        if tractcn >=5:
                            test1=0.5
                        elif tractcn<5:
                            test1=test
                        if tractcn>=6:
                            score4=0
                        else:
                            score4=test1
                        min234=min(score2,score3,score4)
                        #score5-- if pos 13 s D/E/S, test=1, else test=0;
                        #if pos 14 is D/E/S, test1=1, else test1=0
                        #same for E, same for S
                        #sum test and test1. ifsum >0, score5=0, else score5=4
                        #
                        test=0
                        test1=0
                        if subseq[12] in ['D','E','S']:
                            test=1
                        if subseq[13] in ['D','E','S']:
                            test1=1
                        sumtt1=test+test1
                        if sumtt1>0:
                            score5=0
                        else:
                            score5=4
                        finalscore=score1+min234+score5
                        if finalscore<4:
                            top_FFATs_genome.append([ID,subseq,finalscore])
        topFFsdf=pd.DataFrame(top_FFATs_genome,columns=['ID','FFAT','score'])
        topFFsdf.to_excel('topFFATsgenome.xlsx')

records = list(SeqIO.parse(r"HumanProteome.fasta", "fasta"))
fasta_sequences = SeqIO.parse(open('HumanProteome.fasta'),'fasta')
residuescores=pd.read_excel('lookupscore.xlsx')
residuescores.set_index('residue',inplace=True)
residuedict=residuescores.to_dict()
getFFATs(residuedict,fasta_sequences)
        
            
