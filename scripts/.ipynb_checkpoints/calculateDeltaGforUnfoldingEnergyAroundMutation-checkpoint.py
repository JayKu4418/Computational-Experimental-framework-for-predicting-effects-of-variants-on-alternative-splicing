#import pandas as pd
import argparse
import logging
import subprocess
#import Functions_Features.functionsToUnfoldAndCalculateEnergyOfRegionOfInterest as fu
import Functions_Features.functionsToDetermineMotifStructure as fds
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-f","--folderToWrite",type=str,help="Input a folder in which you want to access structural ensemble files into")
parser.add_argument("-t","--mutType",type=str,help="Input whether Muts or SNPs")
parser.add_argument("-c","--foldType",type=str,help="Input whether StructuralEnsemble or MFEstructures")
parser.add_argument("-r","--regionToUnfoldTitle",type=str,help="Input the region to unfold title for folder")
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")
parser.add_argument("-w","--window",type=int,help="Specify window to unfold around mutation site")
parser.add_argument("-b","--intronexonboundary", type=int,help="Input the location of the last position of the exon in the wildtype fasta file")
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
WINDOW=args.window
INTRON_EXON_BOUNDARY=args.intronexonboundary
DMSORNOT=args.dmsornot
SEED=str(args.seed)
GENE=args.genename
STATISTIC=args.statistic

import os
if not os.path.exists(TMP_FOLDER+regionToUnfoldFolderTitle+"_"+MUTTYPE+"/"):
    os.makedirs(TMP_FOLDER+regionToUnfoldFolderTitle+"_"+MUTTYPE+"/")

print(MUTATION_FILE.split("/")[2].split(".")[0])

# Open up the mutation file which contains the position of mutations 
with open(MUTATION_FILE) as f:
    mutations = [line.strip().split("\t") for line in f]

to_write=[]
# Go through every mutation, run parition, followed by stochastic to generate ensemble of 1000, next convert the dot bracket structures to the element representation
for mut in mutations:
    mutID = mut[0]
    if mutID=="WT":
        continue
    print(mutID)
    # If a double mutation,
    if len(mut)>5:
        if int(mut[1]) < 0:
            pos1 = INTRON_EXON_BOUNDARY+int(mut[1])
        else:
            pos1 = INTRON_EXON_BOUNDARY-1+int(mut[1])
        if int(mut[4]) < 0:
            pos2 = INTRON_EXON_BOUNDARY+int(mut[4])
        else:
            pos2 = INTRON_EXON_BOUNDARY-1+int(mut[4])
        
        # Calculate the WT energy difference in unfolding
        average_energy_WT = fds.getAverageDeltaGUnfoldingForTwoSetsOfCoords(TMP_FOLDER,FOLDTYPE,MUTTYPE,regionToUnfoldFolderTitle+"_"+MUTTYPE+"/","WT",pos1 - WINDOW,pos1 + WINDOW,pos2 - WINDOW,pos2 + WINDOW,GENE,SEED,STATISTIC,DMSORNOT)
        
        # Make modifications to start and ends based on whether insertions or deletions
        WTbase1 = mut[2]
        MUTbase1 = mut[3]
        WTbase2 = mut[5]
        MUTbase2 = mut[6]
        
        # Determine what type of mutation for first mutations
    
        firstBaseSubstitute=False
        firstBaseInsertion=False
        firstBaseDeletion=False
    
        if WTbase1!="NA" and MUTbase1!="NA":
            firstBaseSubstitute=True
        elif WTbase1!="NA":
            firstBaseDeletion=True
        else:
            firstBaseInsertion=True
        
        # Determine what type of mutation for second mutations
    
        secondBaseSubstitute=False
        secondBaseInsertion=False
        secondBaseDeletion=False
    
        if WTbase2!="NA" and MUTbase2!="NA":
            secondBaseSubstitute=True
        elif WTbase2!="NA":
            secondBaseDeletion=True
        else:
            secondBaseInsertion=True
        
        # get end coord based on what type of mutation first mut is 
        if firstBaseSubstitute:
            val1ToAdd = 0
        elif firstBaseInsertion:
            val1ToAdd = len(MUTbase1)
        elif firstBaseDeletion:
            val1ToAdd = -len(WTbase1)
        
        # get amount of bases to add based on what type of mutation second mut is 
        if secondBaseSubstitute:
            val2ToAdd=0
        elif secondBaseInsertion:
            val2ToAdd=len(MUTbase2)
        elif secondBaseDeletion:
            val2ToAdd= -len(WTbase2)
        
        # Always the same for all
        START1 = pos1 - WINDOW
        END1 = pos1 + val1ToAdd + WINDOW
        # Start2 position is always whatever adjustments due to mut1 
        START2 = pos2 + val1ToAdd - WINDOW
        END2 = pos2 + val1ToAdd + val2ToAdd + WINDOW
        
        average_energy_MUT = fds.getAverageDeltaGUnfoldingForTwoSetsOfCoords(TMP_FOLDER,FOLDTYPE,MUTTYPE,regionToUnfoldFolderTitle+"_"+MUTTYPE+"/",mutID,START1,END1,START2,END2,GENE,SEED,STATISTIC,DMSORNOT)
        
    # If a single mutation
    else:
        # Get position with respect to the sequence provided
        if int(mut[1]) < 0:
            position = INTRON_EXON_BOUNDARY+int(mut[1])
        else:
            position = INTRON_EXON_BOUNDARY-1+int(mut[1])
        
        average_energy_WT = fds.getAverageDeltaGUnfoldingForCoords(TMP_FOLDER,FOLDTYPE,MUTTYPE,regionToUnfoldFolderTitle+"_"+MUTTYPE+"/","WT",position - WINDOW,position + WINDOW,GENE,SEED,STATISTIC,DMSORNOT)
        
        # Get the WT and MUT bases to determine if insertion, deletion or substituition 
        WTbase = mut[2]
        MUTbase = mut[3]
        
        # adjust end of mutation position so that same window is compared between WT and mutation
        # If deletion
        if MUTbase=="NA":
            END_MUT = position + WINDOW - len(WTbase)
        # If insertion
        elif WTbase=="NA":
            END_MUT = position + WINDOW + len(MUTbase)
        # If substitution
        else:
            END_MUT = position + WINDOW
        
        average_energy_MUT = fds.getAverageDeltaGUnfoldingForCoords(TMP_FOLDER,FOLDTYPE,MUTTYPE,regionToUnfoldFolderTitle+"_"+MUTTYPE+"/",mutID,position - WINDOW,END_MUT,GENE,SEED,STATISTIC,DMSORNOT)
    
    to_write.append([mutID,str(average_energy_MUT),str(average_energy_WT)])
    
with open(TMP_FOLDER+MUTATION_FILE.split("/")[2].split(".")[0]+"_"+STATISTIC+"DeltaGfor"+regionToUnfoldFolderTitle+"_"+FOLDTYPE+".tsv","w") as fw:
    for i in to_write:
        fw.write("\t".join(i))
        fw.write("\n")  