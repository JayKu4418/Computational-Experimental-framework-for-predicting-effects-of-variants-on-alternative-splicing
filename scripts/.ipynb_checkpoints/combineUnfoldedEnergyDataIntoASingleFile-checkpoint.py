import pandas as pd

DATAFOLDER="/data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/"
"""
for stage in ["PreBComplex_6qx9","BComplex_5o9z","PreBactComplex_7abf","BactComplex_5z56"]:
    all_data = []
    colNames_data = []
    data_WT = pd.read_csv(DATAFOLDER+"Unfold_LengthRNAinSpliceosome_"+stage+"_Muts/WT_RsampleOutput_NoMaxDist_Seed50_UnfoldedMotif_EnergyDifference.txt",header=None)
    all_data.append(list(data_WT[0].values))
    colNames_data.append("WT")
    for mut in range(1,60):
        data_MUT = pd.read_csv(DATAFOLDER+"Unfold_LengthRNAinSpliceosome_"+stage+"_Muts/Mut"+str(mut)+"_RsampleOutput_NoMaxDist_Seed50_UnfoldedMotif_EnergyDifference.txt",header=None)
        all_data.append(list(data_MUT[0].values))
        colNames_data.append("Mut"+str(mut))
    data_ToWrite = pd.DataFrame(all_data)
    data_ToWrite_Transposed = data_ToWrite.transpose()
    data_ToWrite_Transposed.columns = colNames_data
    
    data_ToWrite_Transposed.to_csv(DATAFOLDER+"AllMuts_Unfold_LengthRNAinSpliceosome_"+stage+".tsv",sep="\t",header=True,index=False)
    
for stage in ["PreBComplex_6qx9","BComplex_5o9z","PreBactComplex_7abf","BactComplex_5z56"]:
    all_data = []
    colNames_data = []
    data_WT = pd.read_csv(DATAFOLDER+"Unfold_LengthRNAinSpliceosome_"+stage+"_SNPs/WT_RsampleOutput_NoMaxDist_Seed50_UnfoldedMotif_EnergyDifference.txt",header=None)
    all_data.append(list(data_WT[0].values))
    colNames_data.append("WT")
    for mut in range(1,59):
        data_MUT = pd.read_csv(DATAFOLDER+"Unfold_LengthRNAinSpliceosome_"+stage+"_SNPs/SNP"+str(mut)+"_RsampleOutput_NoMaxDist_Seed50_UnfoldedMotif_EnergyDifference.txt",header=None)
        all_data.append(list(data_MUT[0].values))
        colNames_data.append("SNP"+str(mut))
    data_ToWrite = pd.DataFrame(all_data)
    data_ToWrite_Transposed = data_ToWrite.transpose()
    data_ToWrite_Transposed.columns = colNames_data
    
    data_ToWrite_Transposed.to_csv(DATAFOLDER+"AllSNPs_Unfold_LengthRNAinSpliceosome_"+stage+".tsv",sep="\t",header=True,index=False)
"""    
for stage in ["BactComplex_5z56"]:
    all_data = []
    colNames_data = []
    data_WT = pd.read_csv(DATAFOLDER+"Unfold_LengthRNAinSpliceosome_"+stage+"_CompleteMutagenesis/WT_RsampleOutput_NoMaxDist_Seed50_UnfoldedMotif_EnergyDifference.txt",header=None)
    all_data.append(list(data_WT[0].values))
    colNames_data.append("WT")
    for mut in range(1,703):
        data_MUT = pd.read_csv(DATAFOLDER+"Unfold_LengthRNAinSpliceosome_"+stage+"_CompleteMutagenesis/ID"+str(mut)+"_RsampleOutput_NoMaxDist_Seed50_UnfoldedMotif_EnergyDifference.txt",header=None)
        all_data.append(list(data_MUT[0].values))
        colNames_data.append("ID"+str(mut))
    data_ToWrite = pd.DataFrame(all_data)
    data_ToWrite_Transposed = data_ToWrite.transpose()
    data_ToWrite_Transposed.columns = colNames_data
    
    data_ToWrite_Transposed.to_csv(DATAFOLDER+"CompleteMutagenesis_Unfold_LengthRNAinSpliceosome_"+stage+".tsv",sep="\t",header=True,index=False)