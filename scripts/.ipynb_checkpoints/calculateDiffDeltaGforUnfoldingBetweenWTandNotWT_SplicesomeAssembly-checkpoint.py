import pandas as pd
import numpy as  np
import sys 

foldingLabel=sys.argv[1]
TypeLabel=sys.argv[2]

# This is finding diffs for splice sites and regions around splice site, where in the table written the first row is WT and the remaining rows are Mutations 
def findDiffBetweenWTandMUT_SingleCols(datafile,writefile):
    with open(datafile) as f:
        data = [line.strip().split("\t") for line in f]
    
    WTval = [float(i[1]) for i in data if i[0]=="WT"][0]
    to_write=[]
    for mut in data:
        diff=WTval-float(mut[1])
        to_write.append([mut[0],str(diff)])
    with open(writefile,"w") as fw:
        fw.write("MutID"+"\t"+"Struc"+"\n")
        for i in to_write:
            fw.write("\t".join(i))
            fw.write("\n")


for stage in ["PreBComplex_6qx9","BComplex_5o9z","BComplex_6ahd","PreBactComplex_7abf","BactComplex_6ff4","BactComplex_5z56"]:

# Invivo data
    findDiffBetweenWTandMUT_SingleCols("../tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_LengthRNAinSpliceosome_"+stage+"_"+foldingLabel+".tsv","../results/InvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforUnfoldingLengthRNAinSpliceosome_"+stage+"_"+foldingLabel+".tsv")

# Exvivo data
    findDiffBetweenWTandMUT_SingleCols("../tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_LengthRNAinSpliceosome_"+stage+"_"+foldingLabel+".tsv","../results/ExvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforUnfoldingLengthRNAinSpliceosome_"+stage+"_"+foldingLabel+".tsv")

# No data
    findDiffBetweenWTandMUT_SingleCols("../tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_LengthRNAinSpliceosome_"+stage+"_"+foldingLabel+".tsv","../results/NoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforUnfoldingLengthRNAinSpliceosome_"+stage+"_"+foldingLabel+".tsv")