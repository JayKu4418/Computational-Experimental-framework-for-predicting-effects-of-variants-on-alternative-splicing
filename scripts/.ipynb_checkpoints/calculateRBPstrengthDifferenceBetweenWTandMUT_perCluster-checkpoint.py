import argparse
import Functions_Features.functionsToDetermineMotifStrength as fdm
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("-t","--folderToWrite",type=str,help="Input a folder in which you want to write files into")
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")
parser.add_argument("-q","--quantile",nargs='?',default=0.95,type=float,help="Input a quantile value to set a threshold strength score for each motif cluster, default is 0.95")
parser.add_argument("-r","--region",type=str,help="Input either exon or intron for region")
args = parser.parse_args()

TMPfolder=args.folderToWrite
MUTATION_FILE=args.mutationfile
QUANTILE=args.quantile
REGION=args.region

strength_threshold_dict=fdm.createRBPclusterThresholdDictionary(TMPfolder,QUANTILE)

with open(MUTATION_FILE) as f:
#with open("../data/MAPT_MUTs_ToTest.tsv") as f:
    mutations=[line.strip().split("\t") for line in f]

#mutsToIgnore=["Mut3","Mut10","Mut33"] 
to_write = []

# Get the IDs for the RBP motifs 
RBPmotif_ids = []
for filename in os.listdir(TMPfolder+"PWMs_RBPs_Compendium/"):
    if "_pwm.txt" in filename:
        RBPmotif_ids.append(filename.split("_")[0])

# Go through each mutation
for mut in mutations:
    mutID=mut[0]
    motifStrengths = fdm.getSumOfMotifScoreDiffsPerRBPmotif(TMPfolder,mutID,REGION,RBPmotif_ids,strength_threshold_dict)
    
    motifStrengths_forMut = [mutID]+motifStrengths
    
    to_write.append(motifStrengths_forMut)
    
with open(TMPfolder+MUTATION_FILE.split("/")[2].split(".")[0]+"_RBPs_Compendium_strengthsDifferences_perRBPmotif_in" + REGION + ".tsv","w") as fw:
#with open(TMPfolder+motifType+"_MUTsToTest_ScoreDifferences.tsv","w") as fw:
    fw.write("MutID")
    for motif in RBPmotif_ids:
        fw.write("\t")
        fw.write(motif)
    fw.write("\n")
    for i in to_write:
        fw.write("\t".join(i))
        fw.write("\n")