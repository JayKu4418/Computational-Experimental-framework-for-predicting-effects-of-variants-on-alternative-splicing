import Functions_Features.functionsToDetermineMotifStrength as fdm
import Functions_Features.functionsToDetermineMotifStructure as fds
import argparse
import os

# This script will calculate the difference between WT and MUT SRE scores 
parser = argparse.ArgumentParser()
parser.add_argument("-t","--structurefolder",type=str,help="Input the structure folder to get structural data")
parser.add_argument("-f","--rbpfolder",type=str,help="Input the RBP folder to get RBP motif data")
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")
parser.add_argument("-q","--quantile",nargs='?',default=0.95,type=float,help="Input a quantile value to set a threshold strength score for each motif cluster, default is 0.95")
parser.add_argument("-i","--boundary",type=int,help="Get the coordinate of the intron exon boundary with respect to the sequence being folded")
args = parser.parse_args()

TMPfolder=args.structurefolder
TMPfolder_RBP=args.rbpfolder
MUTATION_FILE=args.mutationfile
QUANTILE=args.quantile
INTRON_EXON_BOUNDARY=args.boundary

# This is a dictionary that holds the threshold strengths for each RBP for each motif type
strength_threshold_dict=fdm.createRBPclusterThresholdDictionary(TMPfolder_RBP,QUANTILE)

# Get the RBPmotifs
# Get the IDs for the RBP motifs 
RBPmotif_ids = []
for filename in os.listdir(TMPfolder_RBP+"PWMs_RBPs_Compendium/"):
    if "_pwm.txt" in filename:
        RBPmotif_ids.append(filename.split("_")[0])


# Mutation file
with open(MUTATION_FILE) as f:
    mutations=[line.strip().split("\t") for line in f]

to_write_toFile=[]    

# Go through every mut
for mut in mutations:
    
    mutID = mut[0]
    # Open up the file that contains the delta Gs of unfolding calculated for all possible SREs for a mutation
    with open(TMPfolder+"MedianDeltaGforUnfolding_RBPs_WTandMUTs/"+mutID+"_UnfoldedRBP_AllCoords_MedianDeltaGUnfolding.tsv") as f:
        deltaGs_for_muts = [line.strip().split("\t") for line in f]
    
    newboundary=INTRON_EXON_BOUNDARY
    # Only applicable to insertions and deletions in exon
    # Insertion because WTbase is NA 
    # Need to increase start site and end site by number of bases inserted 
    if mutID!="WT" and int(mut[1]) < 0 and mut[2]=="NA":
        newboundary=INTRON_EXON_BOUNDARY+len(mut[3])
    # Deletion because Mutbase is NA
    # Need to decrease the start site and end site by number of bases deleted
    if mutID!="WT" and int(mut[1]) < 0 and mut[3]=="NA":
        newboundary=INTRON_EXON_BOUNDARY-len(mut[2])
    
    medianDeltaGs_Exon=fds.getMedianDeltaGUnfoldingPerPerRBP(TMPfolder_RBP,mutID,RBPmotif_ids,"Exon",strength_threshold_dict,deltaGs_for_muts,newboundary)
    medianDeltaGs_Intron=fds.getMedianDeltaGUnfoldingPerPerRBP(TMPfolder_RBP,mutID,RBPmotif_ids,"Intron",strength_threshold_dict,deltaGs_for_muts,newboundary)
    
    medianDeltaGs_all = medianDeltaGs_Exon + medianDeltaGs_Intron
    
    to_write_toFile.append([mutID]+medianDeltaGs_all)
    
# Write file
with open(TMPfolder+MUTATION_FILE.split("/")[2].split(".")[0]+"_MedianDeltaGforUnfoldingRBPs_ExonAndIntron.tsv","w") as fw:
    fw.write("MutID")
    for i in RBPmotif_ids:
        fw.write("\t")
        fw.write(i+"_Exon")
    for i in RBPmotif_ids:
        fw.write("\t")
        fw.write(i+"_Intron")
    fw.write("\n")
    for i in to_write_toFile:
        fw.write("\t".join(i))
        fw.write("\n")