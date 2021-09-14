import pandas as pd
import numpy as  np

# This is finding diffs for 5' SS
with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_NoDMSdata_PooledRep/MAPT_SNPs_ToTest_MedianDeltaGforUnfold_20bpWindowExonIntronJunction.tsv") as f:
    singleMuts=[line.strip().split("\t") for line in f]
    
WTval = [float(i[1]) for i in singleMuts if i[0]=="WT"][0]

to_write=[]
for mut in singleMuts:
    diff=WTval-float(mut[1])
    to_write.append([mut[0],str(diff)])
    
with open("../results/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_NoDMSdata_PooledRep/MAPT_ModelFeature_SNPsToTest_DiffInMedianDeltaGforUnfolding20bpWindowExonIntronJunction.tsv","w") as fw:
    fw.write("MutID"+"\t"+"SS"+"\n")
    for i in to_write:
        fw.write("\t".join(i))
        fw.write("\n")

"""
# This is finding diffs for SREs        
with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_InvivoDMSdata_PooledRep/MAPT_SNPs_ToTest_MedianDeltaGforUnfoldingSREs_PerCluster.tsv") as f:
    singleMuts=[line.strip().split("\t") for line in f]  

WTval = [i for i in singleMuts if i[0]=="WT"][0][1:]
WTval = np.array([float(i) for i in WTval])

to_write=[]
for mut in singleMuts:
    mutID=mut[0]
    MUTval = np.array([float(i) for i in mut[1:]])
    diff = list(WTval-MUTval)
    diff = [str(i) for i in diff]
    to_write.append([mutID]+diff)

dict_NumCluster={"ESE":8,"ESS":7,"ISE":7,"ISS":8}
    
with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_InvivoDMSdata_PooledRep/MAPT_ModelFeature_SNPsToTest_DiffInMedianDeltaGforUnfoldingSREs_PerCluster.tsv","w") as fw:
    fw.write("MutID")
    for motifType in ["ESE","ESS","ISE","ISS"]:
        for cluster in range(1,dict_NumCluster[motifType]+1):
            fw.write("\t")
            fw.write(motifType+"_Cluster"+str(cluster))
    fw.write("\n")
    for i in to_write:
        fw.write("\t".join(i))
        fw.write("\n")



struc_scores = pd.read_csv("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_InvivoDMSdata_PooledRep/MAPT_ModelFeature_SNPsToTest_DiffInMedianDeltaGforUnfoldingSREs_PerCluster.tsv",sep="\t",header=0)
struc_scores.index = struc_scores["MutID"]

struc_Muts=[]
for motifType in ["ESE","ESS","ISE","ISS"]:

    motif_columns = [i for i in struc_scores.columns.values if motifType in i]

    struc_scores_motif = struc_scores[motif_columns]

    struc_Muts.append(struc_scores_motif.apply(axis=1,func=np.mean))

struc_Muts_DF = pd.DataFrame(struc_Muts)
struc_Muts_DF = struc_Muts_DF.transpose()
struc_Muts_DF.columns=["ESE","ESS","ISE","ISS"]
struc_Muts_DF = struc_Muts_DF.assign(MutID=struc_Muts_DF.index.values)

struc_Muts_DF.iloc[:,[4,0,1,2,3]].to_csv("../results/MAPT_ModelFeature_SNPsToTest_DiffInMedianDeltaGforUnfoldingSREs_MeanOfClustersPerMotifType.tsv",sep="\t",header=True,index=False)

# This is finding diffs for RBPs

with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_InvivoDMSdata_PooledRep/MAPT_SingleMutations_TrainingSet_MedianDeltaGforUnfoldingRBPs_ExonAndIntron.tsv") as f:
    struc_rbp_data_single = [line.strip().split("\t") for line in f]
with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_InvivoDMSdata_PooledRep/MAPT_DoubleMutations_TrainingSet_MedianDeltaGforUnfoldingRBPs_ExonAndIntron.tsv") as f:
    struc_rbp_data_double = [line.strip().split("\t") for line in f]
    
struc_rbp_data_all = struc_rbp_data_single[1:]+struc_rbp_data_double[1:]

WT_struc_rbp_data = np.array([float(i) for i in struc_rbp_data_single[1][1:]])

struc_diff_data_to_write=[]
for mutdata in struc_rbp_data_all:
    MUT_struc_rbp_data = np.array([float(i) for i in mutdata[1:]])
    diff_data_for_mut = WT_struc_rbp_data-MUT_struc_rbp_data
    diff_data_for_mut_str = [str(i) for i in diff_data_for_mut]
    struc_diff_data_to_write.append([mutdata[0]] + diff_data_for_mut_str)
    
with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_InvivoDMSdata_PooledRep/MAPT_ModelFeature_AllMutations_WTvsMUT_DiffInMedianDeltaGforUnfoldingRBPs_ExonAndIntron.tsv","w") as fw:
    fw.write("\t".join(struc_rbp_data_single[0]))
    fw.write("\n")
    for i in struc_diff_data_to_write:
        fw.write("\t".join(i))
        fw.write("\n")
        
struc_scores = pd.read_csv("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_InvivoDMSdata_PooledRep/MAPT_ModelFeature_AllMutations_WTvsMUT_DiffInMedianDeltaGforUnfoldingRBPs_ExonAndIntron.tsv",sep="\t",header=0)
struc_scores.index = struc_scores["MutID"]

struc_Muts=[]
for motifType in ["Exon","Intron"]:

    motif_columns = [i for i in struc_scores.columns.values if motifType in i]

    struc_scores_motif = struc_scores[motif_columns]

    struc_Muts.append(struc_scores_motif.apply(axis=1,func=np.mean))

struc_Muts_DF = pd.DataFrame(struc_Muts)
struc_Muts_DF = struc_Muts_DF.transpose()
struc_Muts_DF.columns=["Exon","Intron"]
struc_Muts_DF = struc_Muts_DF.assign(MutID=struc_Muts_DF.index.values)

struc_Muts_DF.iloc[:,[2,0,1]].to_csv("../results/MAPT_ModelFeature_AllMutations_WTvsMUT_DiffInMedianDeltaGforUnfoldingRBPs_MeanOfExonOrIntron.tsv",sep="\t",header=True,index=False)
"""