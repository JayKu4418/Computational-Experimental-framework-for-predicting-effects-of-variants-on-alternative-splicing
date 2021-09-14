import Functions_Features.functionsToDetermineMotifStrength as fdm
import Functions_Features.functionsToDetermineMotifStructure as fds
from joblib import Parallel, delayed
import multiprocessing
import argparse

# This script will calculate the difference between WT and MUT SRE scores 
parser = argparse.ArgumentParser()
parser.add_argument("-w","--structurefolder",type=str,help="Input the structure folder to get structural data")
parser.add_argument("-f","--srefolder",type=str,help="Input the SRE folder to get SRE motif data")
parser.add_argument("-t","--folderTitle",type=str,help="Input the title of the mutation file")
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")
parser.add_argument("-q","--quantile",nargs='?',default=0.95,type=float,help="Input a quantile value to set a threshold strength score for each motif cluster, default is 0.95")
parser.add_argument("-p","--dmsornot",type=str,help="Indicate if Rsample or Partition")
parser.add_argument("-s","--seed", nargs='?',default=50,type=int,help="Input the seed for generating sample structures using stochastic, this is optional, default is 50")
parser.add_argument("-g","--genename",type=str,help="Input genename for CTFile")
parser.add_argument("-i","--boundary",type=int,help="Get the coordinate of the intron exon boundary with respect to the sequence being folded")

args = parser.parse_args()

TMPfolder=args.structurefolder
TMPfolder_SRE=args.srefolder
folderTitle=args.folderTitle
GENE=args.genename
SEED=args.seed
DMSORNOT=args.dmsornot
MUTATION_FILE=args.mutationfile
QUANTILE=args.quantile
INTRON_EXON_BOUNDARY=args.boundary

# This is dictionary for number of clusters per motif type
dict_NumCluster={"ESE":8,"ESS":7,"ISE":7,"ISS":8}
# This is a dictionary that holds the threshold strengths for each cluster for each motif type
strength_threshold_dict=fdm.createSREclusterThresholdDictionary(TMPfolder_SRE,dict_NumCluster,QUANTILE)

# This function will perform the unfolding of motif coordinates and calculate the deltaG of unfolding of each motif for given mutation
def unfoldSREmotifsForMutation(mutation):
    
    mutID=mutation[0]
    
    newboundary=INTRON_EXON_BOUNDARY
    # Only applicable to insertions and deletions in exon
    # Insertion because WTbase is NA 
    # Need to increase start site and end site by number of bases inserted 
    if mutID!="WT" and int(mutation[1]) < 0 and mutation[2]=="NA":
        newboundary=INTRON_EXON_BOUNDARY+len(mutation[3])
    # Deletion because Mutbase is NA
    # Need to decrease the start site and end site by number of bases deleted
    if mutID!="WT" and int(mutation[1]) < 0 and mutation[3]=="NA":
        newboundary=INTRON_EXON_BOUNDARY-len(mutation[2])
    
    
    coordsToUnfold_ESE=fds.getCoordsForSREMotifs(TMPfolder_SRE,folderTitle,mutID,"ESE",dict_NumCluster["ESE"],strength_threshold_dict,newboundary)
    coordsToUnfold_ESS=fds.getCoordsForSREMotifs(TMPfolder_SRE,folderTitle,mutID,"ESS",dict_NumCluster["ESS"],strength_threshold_dict,newboundary)
    coordsToUnfold_ISE=fds.getCoordsForSREMotifs(TMPfolder_SRE,folderTitle,mutID,"ISE",dict_NumCluster["ISE"],strength_threshold_dict,newboundary)
    coordsToUnfold_ISS=fds.getCoordsForSREMotifs(TMPfolder_SRE,folderTitle,mutID,"ISS",dict_NumCluster["ISS"],strength_threshold_dict,newboundary)
    # Combine all the coordinates
    combinedCoords_Mut = coordsToUnfold_ESE+coordsToUnfold_ESS+coordsToUnfold_ISE+coordsToUnfold_ISS
    # Get the unique coordinates
    unique_combinedCoords_Mut=[list(x) for x in set(tuple(x) for x in combinedCoords_Mut)]
    
    to_write=[]
    
    for coord in unique_combinedCoords_Mut:
        # Make the 1 based coordinate
        startCoord=int(coord[0])+1
        endCoord=int(coord[1])
        deltaGDiff_forCoord = fds.getMedianDeltaGUnfoldingForCoords(TMPfolder,folderTitle,"Unfold_SREs_"+folderTitle+"/",mutID,startCoord,endCoord,GENE,DMSORNOT,SEED)
        to_write.append([mutID,str(coord[0]),str(coord[1]),str(deltaGDiff_forCoord)])
        
    with open(TMPfolder+"MedianDeltaGforUnfolding_SREs_"+folderTitle+"/"+mutID+"_UnfoldedSRE_AllCoords_MedianDeltaGUnfolding.tsv","w") as fw:
        for i in to_write:
            fw.write("\t".join(i))
            fw.write("\n")
            
# Mutation file
with open(MUTATION_FILE) as f:
    mutations=[line.strip().split("\t") for line in f]

# Get number of cores
num_cores = multiprocessing.cpu_count()

# Run all the mutations in parallel
Parallel(n_jobs=num_cores)(delayed(unfoldSREmotifsForMutation)(mut) for mut in mutations)
