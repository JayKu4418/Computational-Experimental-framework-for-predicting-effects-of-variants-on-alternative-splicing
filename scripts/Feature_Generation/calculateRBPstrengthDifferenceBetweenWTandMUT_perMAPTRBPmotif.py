import argparse
import Functions_Features.functionsToDetermineMotifStrength as fdm
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")
parser.add_argument("-q","--quantile",nargs='?',default=0.95,type=float,help="Input a quantile value to set a threshold strength score for each motif cluster, default is 0.95")
args = parser.parse_args()

MUTATION_FILE=args.mutationfile
QUANTILE=args.quantile

strength_threshold_dict=fdm.createRBPclusterThresholdDictionary(QUANTILE)

with open(MUTATION_FILE) as f:
#with open("../data/MAPT_MUTs_ToTest.tsv") as f:
    mutations=[line.strip().split("\t") for line in f]

#mutsToIgnore=["Mut3","Mut10","Mut33"] 
to_write = []

# Get the IDs for the RBP motifs 
RBPmotif_ids = []
for filename in os.listdir("../data/MAPT_RBPs/"):
    if "_pwm.txt" in filename:
        RBPmotif_ids.append(filename.split("_")[0])

# Go through each mutation
for mut in mutations:
    mutID=mut[0]
    motifStrengths = fdm.getSumOfMotifScoreDiffsPerRBPmotif(mutID,RBPmotif_ids,strength_threshold_dict)
    
    motifStrengths_forMut = [mutID]+motifStrengths
    
    to_write.append(motifStrengths_forMut)
    
with open("../tmp/"+MUTATION_FILE.split("/")[2].split(".")[0]+"_MAPT_RBPs_strengthsDifferences_perRBPmotif.tsv","w") as fw:
#with open(TMPfolder+motifType+"_MUTsToTest_ScoreDifferences.tsv","w") as fw:
    fw.write("MutID")
    for motif in RBPmotif_ids:
        fw.write("\t")
        fw.write(motif)
    fw.write("\n")
    for i in to_write:
        fw.write("\t".join(i))
        fw.write("\n")