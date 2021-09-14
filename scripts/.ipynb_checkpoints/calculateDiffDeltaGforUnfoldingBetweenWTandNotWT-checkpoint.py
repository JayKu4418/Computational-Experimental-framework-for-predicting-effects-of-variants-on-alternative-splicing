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

# Invivo data
findDiffBetweenWTandMUT_SingleCols("../tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_SpliceSite_"+foldingLabel+".tsv","../results/InvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforUnfolding5pSpliceSite_"+foldingLabel+".tsv")
findDiffBetweenWTandMUT_SingleCols("../tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_MaxLengthRNAinSpliceosome_"+foldingLabel+".tsv","../results/InvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforUnfoldingMaxLengthRNAinSpliceosome_"+foldingLabel+".tsv")
findDiffBetweenWTandMUT_SingleCols("../tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_20basesAroundSpliceSite_"+foldingLabel+".tsv","../results/InvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforUnfolding20BasesAround5pSpliceSite_"+foldingLabel+".tsv")
findDiffBetweenWTandMUT_SingleCols("../tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforEnsemble.tsv","../results/InvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforEnsemble.tsv")
findDiffBetweenWTandMUT_SingleCols("../tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_DeltaGforMFEstructure.tsv","../results/InvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInDeltaGforMFEstructure.tsv")

# Exvivo data
findDiffBetweenWTandMUT_SingleCols("../tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_SpliceSite_"+foldingLabel+".tsv","../results/ExvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforUnfolding5pSpliceSite_"+foldingLabel+".tsv")
findDiffBetweenWTandMUT_SingleCols("../tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_MaxLengthRNAinSpliceosome_"+foldingLabel+".tsv","../results/ExvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforUnfoldingMaxLengthRNAinSpliceosome_"+foldingLabel+".tsv")
findDiffBetweenWTandMUT_SingleCols("../tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_20basesAroundSpliceSite_"+foldingLabel+".tsv","../results/ExvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforUnfolding20BasesAround5pSpliceSite_"+foldingLabel+".tsv")
findDiffBetweenWTandMUT_SingleCols("../tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforEnsemble.tsv","../results/ExvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforEnsemble.tsv")
findDiffBetweenWTandMUT_SingleCols("../tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_DeltaGforMFEstructure.tsv","../results/ExvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInDeltaGforMFEstructure.tsv")

# No data
findDiffBetweenWTandMUT_SingleCols("../tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_SpliceSite_"+foldingLabel+".tsv","../results/NoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforUnfolding5pSpliceSite_"+foldingLabel+".tsv")
findDiffBetweenWTandMUT_SingleCols("../tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_MaxLengthRNAinSpliceosome_"+foldingLabel+".tsv","../results/NoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforUnfoldingMaxLengthRNAinSpliceosome_"+foldingLabel+".tsv")
findDiffBetweenWTandMUT_SingleCols("../tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_20basesAroundSpliceSite_"+foldingLabel+".tsv","../results/NoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforUnfolding20BasesAround5pSpliceSite_"+foldingLabel+".tsv")
findDiffBetweenWTandMUT_SingleCols("../tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforEnsemble.tsv","../results/NoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforEnsemble.tsv")
findDiffBetweenWTandMUT_SingleCols("../tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_DeltaGforMFEstructure.tsv","../results/NoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInDeltaGforMFEstructure.tsv")



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
        fw.write("MutID"+"\t"+"Struc"+"\n")
        for i in to_write:
            fw.write("\t".join(i))
            fw.write("\n")

# Invivo data
findDiffBetweenWTandMUT_DoubleCols("../tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_LocalRegion_"+foldingLabel+".tsv","../results/InvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforLocalRegion_"+foldingLabel+".tsv")
# Exvivo data
findDiffBetweenWTandMUT_DoubleCols("../tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_LocalRegion_"+foldingLabel+".tsv","../results/ExvivoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforLocalRegion_"+foldingLabel+".tsv")
# No data
findDiffBetweenWTandMUT_DoubleCols("../tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/MAPT_"+TypeLabel+"_MedianDeltaGforUnfold_LocalRegion_"+foldingLabel+".tsv","../results/NoDMSdata_PooledRep_Exon10Intron10/MAPT_ModelFeature_"+TypeLabel+"_WTvsMUT_DiffInMedianDeltaGforLocalRegion_"+foldingLabel+".tsv")
      
