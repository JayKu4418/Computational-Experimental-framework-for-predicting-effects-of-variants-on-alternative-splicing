import numpy as np
import pandas as pd

# This is for Splice sites

with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_SREmotifs/MAPT_SNPs_ToTest_5pSpliceSitesScores.tsv") as f:
    singleMuts=[line.strip().split("\t") for line in f][1:]

WTval = [float(i[1]) for i in singleMuts if i[0]=="WT"][0]


to_write=[]
for mut in singleMuts:
    diff=float(mut[1])-WTval
    to_write.append([mut[0],str(diff)])
    
with open("../results/MAPT_ModelFeature_SNPsToTest_DiffIn5pSpliceSiteStrength.tsv","w") as fw:
    fw.write("MutID"+"\t"+"SS"+"\n")
    for i in to_write:
        fw.write("\t".join(i))
        fw.write("\n")


# This is for SREs
sre_scores = pd.read_csv("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_SREmotifs/MAPT_SNPs_ToTest_SREstrengthsDifferences_perCluster.tsv",header=0,sep="\t")

sre_scores.index = sre_scores["MutID"]

strengths_Muts=[]
for motifType in ["ESE","ESS","ISE","ISS"]:

    motif_columns = [i for i in sre_scores.columns.values if motifType in i]

    sre_scores_motif = sre_scores[motif_columns]

    strengths_Muts.append(sre_scores_motif.apply(axis=1,func=np.mean))

strengths_Muts_DF = pd.DataFrame(strengths_Muts)
strengths_Muts_DF = strengths_Muts_DF.transpose()
strengths_Muts_DF.columns=["ESE","ESS","ISE","ISS"]
strengths_Muts_DF = strengths_Muts_DF.assign(MutID=strengths_Muts_DF.index.values)

strengths_Muts_DF.iloc[:,[4,0,1,2,3]].to_csv("../results/MAPT_ModelFeature_SNPsToTest_DiffStrengthOfSREs_MeanOfClustersPerMotifType.tsv",sep="\t",header=True,index=False)

# This is for RBPs
"""
with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_RBPmotifs/MAPT_SingleMutations_TrainingSet_RBPs_Compendium_strengthsDifferences_perRBPmotif_inExon.tsv") as f:
    singleMuts=[line.strip().split("\t") for line in f]

with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_RBPmotifs/MAPT_DoubleMutations_TrainingSet_RBPs_Compendium_strengthsDifferences_perRBPmotif_inExon.tsv") as f:
    doubleMuts=[line.strip().split("\t") for line in f]
    
to_Write_Exon = singleMuts[1:] + doubleMuts[1:]

with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_RBPmotifs/MAPT_ModelFeature_AllMutations_WTvsMUT_DiffStrengthOfRBPs_RBPs_Compendium_perRBPmotif_inExon.tsv","w") as fw:
    fw.write("\t".join(singleMuts[0]))
    fw.write("\n")
    for i in to_Write_Exon:
        fw.write("\t".join(i))
        fw.write("\n")

with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_RBPmotifs/MAPT_SingleMutations_TrainingSet_RBPs_Compendium_strengthsDifferences_perRBPmotif_inIntron.tsv") as f:
    singleMuts=[line.strip().split("\t") for line in f]

with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_RBPmotifs/MAPT_DoubleMutations_TrainingSet_RBPs_Compendium_strengthsDifferences_perRBPmotif_inIntron.tsv") as f:
    doubleMuts=[line.strip().split("\t") for line in f]
    
to_Write_Intron = singleMuts[1:] + doubleMuts[1:]

with open("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_RBPmotifs/MAPT_ModelFeature_AllMutations_WTvsMUT_DiffStrengthOfRBPs_RBPs_Compendium_perRBPmotif_inIntron.tsv","w") as fw:
    fw.write("\t".join(singleMuts[0]))
    fw.write("\n")
    for i in to_Write_Intron:
        fw.write("\t".join(i))
        fw.write("\n")
        
strengths_Muts=[]
        
rbp_scores_exon = pd.read_csv("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_RBPmotifs/MAPT_ModelFeature_AllMutations_WTvsMUT_DiffStrengthOfRBPs_RBPs_Compendium_perRBPmotif_inExon.tsv",header=0,sep="\t")

#rbp_scores_exon.index = rbp_scores_exon["MutID"]

motif_columns = [i for i in rbp_scores_exon.columns.values if "MutID" not in i]

strengths_Muts.append(rbp_scores_exon[motif_columns].apply(axis=1,func=np.mean))

rbp_scores_intron = pd.read_csv("../tmp/Model_MAPTsplicing_MAPTExon10Intron10withDMSdata_RBPmotifs/MAPT_ModelFeature_AllMutations_WTvsMUT_DiffStrengthOfRBPs_RBPs_Compendium_perRBPmotif_inIntron.tsv",header=0,sep="\t")

motif_columns = [i for i in rbp_scores_intron.columns.values if "MutID" not in i]

strengths_Muts.append(rbp_scores_intron[motif_columns].apply(axis=1,func=np.mean))

strengths_Muts_DF = pd.DataFrame(strengths_Muts)
strengths_Muts_DF = strengths_Muts_DF.transpose()
strengths_Muts_DF.columns=["Exon","Intron"]
strengths_Muts_DF = strengths_Muts_DF.assign(MutID=rbp_scores_exon["MutID"].values)

strengths_Muts_DF.iloc[:,[2,0,1]].to_csv("../results/MAPT_ModelFeature_AllMutations_WTvsMUT_DiffStrengthOfRBPs_MeanOfExonOrIntron.tsv",sep="\t",header=True,index=False)
"""

