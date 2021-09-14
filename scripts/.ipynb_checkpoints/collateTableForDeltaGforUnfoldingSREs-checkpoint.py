import Functions_Features.functionsToDetermineMotifStrength as fdm
import Functions_Features.functionsToDetermineMotifStructure as fds
import argparse

# This script will calculate the difference between WT and MUT SRE scores 
parser = argparse.ArgumentParser()
parser.add_argument("-w","--structurefolder",type=str,help="Input the structure folder to get structural data")
parser.add_argument("-f","--srefolder",type=str,help="Input the SRE folder to get SRE motif data")
parser.add_argument("-t","--folderTitle",type=str,help="Input the title of the mutation file")
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")
parser.add_argument("-q","--quantile",nargs='?',default=0.95,type=float,help="Input a quantile value to set a threshold strength score for each motif cluster, default is 0.95")
parser.add_argument("-i","--boundary",type=int,help="Get the coordinate of the intron exon boundary with respect to the sequence being folded")
args = parser.parse_args()

TMPfolder=args.structurefolder
TMPfolder_SRE=args.srefolder
folderTitle=args.folderTitle
MUTATION_FILE=args.mutationfile
QUANTILE=args.quantile
INTRON_EXON_BOUNDARY=args.boundary

# This is dictionary for number of clusters per motif type
dict_NumCluster={"ESE":8,"ESS":7,"ISE":7,"ISS":8}
# This is a dictionary that holds the threshold strengths for each cluster for each motif type
strength_threshold_dict=fdm.createSREclusterThresholdDictionary(TMPfolder_SRE,dict_NumCluster,QUANTILE)
            
# Mutation file
with open(MUTATION_FILE) as f:
    mutations=[line.strip().split("\t") for line in f]

to_write_toFile=[]    

# Go through every mut
for mut in mutations:
    
    mutID = mut[0]
    # Open up the file that contains the delta Gs of unfolding calculated for all possible SREs for a mutation
    with open(TMPfolder+"MedianDeltaGforUnfolding_SREs_"+folderTitle+"/"+mutID+"_UnfoldedSRE_AllCoords_MedianDeltaGUnfolding.tsv") as f:
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
    
    
    medianDeltaGs_ESE=fds.getMedianDeltaGUnfoldingPerSRECluster(TMPfolder_SRE,folderTitle,mutID,"ESE",dict_NumCluster["ESE"],strength_threshold_dict,deltaGs_for_muts,newboundary)
    medianDeltaGs_ESS=fds.getMedianDeltaGUnfoldingPerSRECluster(TMPfolder_SRE,folderTitle,mutID,"ESS",dict_NumCluster["ESS"],strength_threshold_dict,deltaGs_for_muts,newboundary)
    medianDeltaGs_ISE=fds.getMedianDeltaGUnfoldingPerSRECluster(TMPfolder_SRE,folderTitle,mutID,"ISE",dict_NumCluster["ISE"],strength_threshold_dict,deltaGs_for_muts,newboundary)
    medianDeltaGs_ISS=fds.getMedianDeltaGUnfoldingPerSRECluster(TMPfolder_SRE,folderTitle,mutID,"ISS",dict_NumCluster["ISS"],strength_threshold_dict,deltaGs_for_muts,newboundary)
    
    medianDeltaGs_all = medianDeltaGs_ESE + medianDeltaGs_ESS + medianDeltaGs_ISE + medianDeltaGs_ISS
    
    to_write_toFile.append([mutID]+medianDeltaGs_all)
    
# Write file
with open(TMPfolder+MUTATION_FILE.split("/")[2].split(".")[0]+"_MedianDeltaGforUnfoldingSREs_PerCluster.tsv","w") as fw:
    for i in to_write_toFile:
        fw.write("\t".join(i))
        fw.write("\n")