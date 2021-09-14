import argparse
import Functions_Features.functionsToDetermineMotifStrength as fdm
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-w","--tmpfolder",type=str,help="Input the upperlevel folder containing folder to Write to")
parser.add_argument("-t","--foldertitle",type=str,help="Input the title of the mutation file")
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")
parser.add_argument("-q","--quantile",nargs='?',default=0.95,type=float,help="Input a quantile value to set a threshold strength score for each motif cluster, default is 0.95")

args = parser.parse_args()

TMPfolder=args.tmpfolder
folderTitle=args.foldertitle
MUTATION_FILE=args.mutationfile
QUANTILE=args.quantile


dict_NumCluster={"ESE":8,"ESS":7,"ISE":7,"ISS":8}
strength_threshold_dict=fdm.createSREclusterThresholdDictionary(TMPfolder,dict_NumCluster,QUANTILE)

with open(MUTATION_FILE) as f:
#with open("../data/MAPT_MUTs_ToTest.tsv") as f:
    mutations=[line.strip().split("\t") for line in f]

#mutsToIgnore=["Mut3","Mut10","Mut33"] 
to_write = []
    
# Go through each mutation
for mut in mutations:
    mutID=mut[0]
    ESE_motifStrengths = fdm.getSumOfMotifScoreDiffsPerSRECluster(TMPfolder,folderTitle,mutID,"ESE",dict_NumCluster["ESE"],strength_threshold_dict)
    ESS_motifStrengths = fdm.getSumOfMotifScoreDiffsPerSRECluster(TMPfolder,folderTitle,mutID,"ESS",dict_NumCluster["ESS"],strength_threshold_dict)
    ISE_motifStrengths = fdm.getSumOfMotifScoreDiffsPerSRECluster(TMPfolder,folderTitle,mutID,"ISE",dict_NumCluster["ISE"],strength_threshold_dict)
    ISS_motifStrengths = fdm.getSumOfMotifScoreDiffsPerSRECluster(TMPfolder,folderTitle,mutID,"ISS",dict_NumCluster["ISS"],strength_threshold_dict)
    
    motifStrengths_forMut = [mutID]+ESE_motifStrengths+ESS_motifStrengths+ISE_motifStrengths+ISS_motifStrengths
    
    to_write.append(motifStrengths_forMut)
    
with open(TMPfolder+MUTATION_FILE.split("/")[2].split(".")[0]+"_SREstrengthsDifferences_perCluster.tsv","w") as fw:
#with open(TMPfolder+motifType+"_MUTsToTest_ScoreDifferences.tsv","w") as fw:
    fw.write("MutID")
    for motifType in ["ESE","ESS","ISE","ISS"]:
        for cluster in range(1,dict_NumCluster[motifType]+1):
            fw.write("\t")
            fw.write(motifType+"_Cluster"+str(cluster))
    fw.write("\n")
    for i in to_write:
        fw.write("\t".join(i))
        fw.write("\n")