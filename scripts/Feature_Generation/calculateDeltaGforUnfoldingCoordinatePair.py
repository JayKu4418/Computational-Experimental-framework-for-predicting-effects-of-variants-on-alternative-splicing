# Script to calculate the unfolding energy of a given coordinate region 
# Used to calculate unfolding energy for each structure in an ensemble for every mutation 
# in the mutation file 
import argparse
import logging
import subprocess
import Functions_Features.functionsToDetermineMotifStructure as fds
import Functions_Features.functionsToGetMutationSeqAndDMSfile as fms
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-f","--folderToWrite",type=str,help="Input a folder in which you want to access structural ensemble files into, make sure to include / at the end")
parser.add_argument("-t","--mutType",type=str,help="Input whether Muts or SNPs")
parser.add_argument("-b","--foldType",type=str,help="Input whether StructuralEnsemble or MFEstructures")
parser.add_argument("-r","--regionToUnfoldTitle",type=str,help="Input the region to unfold title for folder")
parser.add_argument("-c","--intronexonboundary", type=int,help="Input the location of the last position of the exon in the wildtype fasta file")
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")
parser.add_argument("-1","--start",type=int,help="1based-coordinate Start site with respect to original fasta file for unfolding")
parser.add_argument("-2","--end",type=int,help="1based-coordinate End site with respect to original fasta file for unfolding")
parser.add_argument("-p","--dmsornot",nargs='?',type=str,help="Indicate if Rsample or Partition or _DMSdata, this is optional")
parser.add_argument("-s","--seed", nargs='?',default=50,type=int,help="Input the seed for generating sample structures using stochastic, this is optional, default is 50")
parser.add_argument("-g","--genename",type=str,help="Input genename for CTFile")
parser.add_argument("-a","--statistic",type=str,help="Input statistic, Mean or Median to calculate")

args = parser.parse_args()

TMP_FOLDER=args.folderToWrite
MUTTYPE=args.mutType
FOLDTYPE=args.foldType
regionToUnfoldFolderTitle=args.regionToUnfoldTitle
MUTATION_FILE=args.mutationfile
START=args.start
END=args.end
DMSORNOT=args.dmsornot
SEED=str(args.seed)
GENE=args.genename
STATISTIC=args.statistic
INTRON_EXON_BOUNDARY=args.intronexonboundary

import os
if not os.path.exists(TMP_FOLDER+regionToUnfoldFolderTitle+"_"+MUTTYPE+"/"):
    os.makedirs(TMP_FOLDER+regionToUnfoldFolderTitle+"_"+MUTTYPE+"/")

print(MUTATION_FILE.split("/")[2].split(".")[0])

# Open up the mutation file which contains the position of mutations 
with open(MUTATION_FILE) as f:
    mutations = [line.strip().split("\t") for line in f]

to_write=[]
# Go through every mutation, get new coordinates with respect to the seq string since 
# mutation coordinates are given with respect to exon-intron junction
for mut in mutations:
    mutID = mut[0]
    # If wuldtype keep coordinates same
    if mutID == "WT":
        newSTART=START
        newEND=END
    # If double mutation
    elif len(mut)>5:
        if int(mut[1]) < 0:
            pos1 = INTRON_EXON_BOUNDARY+int(mut[1])
        else:
            pos1 = INTRON_EXON_BOUNDARY-1+int(mut[1])
        if int(mut[4]) < 0:
            pos2 = INTRON_EXON_BOUNDARY+int(mut[4])
        else:
            pos2 = INTRON_EXON_BOUNDARY-1+int(mut[4])
        
        WTbase1 = mut[2]
        MUTbase1 = mut[3]
        WTbase2 = mut[5]
        MUTbase2 = mut[6]
                
        [newSTART,newEND]=fms.getNewStartAndStopCoordinatesForDoubleMutations(WTbase1,MUTbase1,WTbase2,MUTbase2,pos1,pos2,START,END)
    # If single mutation
    else:
        if int(mut[1]) < 0:
            position = INTRON_EXON_BOUNDARY+int(mut[1])
        else:
            position = INTRON_EXON_BOUNDARY-1+int(mut[1])
        
        newSTART=START
        newEND=END
        # Only applicable to insertions and deletions in exon
        # Insertion because WTbase is NA 
        # Need to increase start site and end site by number of bases inserted 
        if position <= END and mut[2]=="NA":
            if position <= START:
                newSTART=newSTART+len(mut[3])
            newEND=newEND+len(mut[3])
        # Deletion because Mutbase is NA
        # Need to decrease the start site and end site by number of bases deleted
        if position <= END and mut[3]=="NA":
            if position <= START:
                newSTART=newSTART-len(mut[2])
            newEND=newEND-len(mut[2])
    
    average_energy = fds.getAverageDeltaGUnfoldingForCoords(TMP_FOLDER,FOLDTYPE,MUTTYPE,regionToUnfoldFolderTitle+"_"+MUTTYPE+"/",mutID,newSTART,newEND,GENE,SEED,STATISTIC,DMSORNOT)
    
    to_write.append([mutID,str(average_energy)])

with open(TMP_FOLDER+MUTATION_FILE.split("/")[2].split(".")[0]+"_"+STATISTIC+"DeltaGfor"+regionToUnfoldFolderTitle+"_"+FOLDTYPE+".tsv","w") as fw:
    for i in to_write:
        fw.write("\t".join(i))
        fw.write("\n")  