import pandas as pd
import numpy as  np
import sys 

TypeLabel=sys.argv[1]

# This is finding diffs for splice sites and regions around splice site, where in the table written the first row is WT and the remaining rows are Mutations 
def findDiffBetweenWTandMUT_SingleCols(datafile,writefile):
    with open(datafile) as f:
        data = [line.strip().split("\t") for line in f][1:]
    
    WTval = [float(i[1]) for i in data if i[0]=="WT"][0]
    to_write=[]
    for mut in data:
        diff=WTval-float(mut[1])
        to_write.append([mut[0],str(diff)])
    with open(writefile,"w") as fw:
        fw.write("MutID"+"\t"+"Strength"+"\n")
        for i in to_write:
            fw.write("\t".join(i))
            fw.write("\n")

def getMeanOfNonZeroValues(x):
    y = list(x.values)
    nz_y = [i for i in y if i!=0]
    if len(nz_y)!=0:
        return(sum(nz_y)/len(nz_y))
    else:
        return(0.0)

# 5p Splice Site
findDiffBetweenWTandMUT_SingleCols("../tmp/MAPT_"+TypeLabel+"_5pSpliceSitesScores.tsv","../results/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffIn5pSpliceSitesScores.tsv")



# These lines of code will calculate the average SRE motif strength diff across all the clusters per category 
sre_scores = pd.read_csv("../tmp/MAPT_"+TypeLabel+"_SREstrengthsDifferences_perCluster.tsv",header=0,sep="\t")

sre_scores.index = sre_scores["MutID"]

strengths_Muts=[]
for motifType in ["ESE","ESS","ISE","ISS"]:

    motif_columns = [i for i in sre_scores.columns.values if motifType in i]

    sre_scores_motif = sre_scores[motif_columns]
    
    strengths_Muts.append(sre_scores_motif.apply(axis=1,func=getMeanOfNonZeroValues))
    #strengths_Muts.append(sre_scores_motif.apply(axis=1,func=np.mean))

strengths_Muts_DF = pd.DataFrame(strengths_Muts)
strengths_Muts_DF = strengths_Muts_DF.transpose()
strengths_Muts_DF.columns=["ESE","ESS","ISE","ISS"]
strengths_Muts_DF = strengths_Muts_DF.assign(MutID=strengths_Muts_DF.index.values)

strengths_Muts_DF.iloc[:,[4,0,1,2,3]].to_csv("../results/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffStrengthOfSREs_MeanOfClustersPerMotifType.tsv",sep="\t",header=True,index=False)

"""
# This is a function for finding diffs for local structural changes around region and the table is written with Mutation in one column and WT data in another column since region analyzed is different for every mutation
def findDiffBetweenWTandMUT_DoubleCols(datafile,writefile):
    with open(datafile) as f:
        data = [line.strip().split("\t") for line in f]
    
    to_write=[]
    to_write.append(["WT","0.0"])
    for mut in data:
        # WT val - MUT val
        diff=float(mut[2])-float(mut[1])
        to_write.append([mut[0],str(diff)])
    with open(writefile,"w") as fw:
        fw.write("MutID"+"\t"+"Strength"+"\n")
        for i in to_write:
            fw.write("\t".join(i))
            fw.write("\n")
"""