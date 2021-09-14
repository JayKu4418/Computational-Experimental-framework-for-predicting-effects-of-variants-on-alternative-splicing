#!/bin/bash

TypeOfFold=$1
MutFile_Name=$2
IsMut=$3
# # Unfold different lengths of RNA in spliceosome in different stages of assembly

# ## Structural Ensemble calculations

# ### Unfold 2 exonic bases and 8 intronic bases for each sequence (WT and mutations) and calculate the difference in delta G -> from splicesome assembly model PreB complex 6QX9
# Invivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_PreBComplex_6qx9 -m ../data/MAPT_${MutFile_Name}.tsv -1 92 -2 101 -p DMSdata -g MAPT -a Median -c 93
# Exvivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_PreBComplex_6qx9 -m ../data/MAPT_${MutFile_Name}.tsv -1 92 -2 101 -p DMSdata -g MAPT -a Median -c 93
# No data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_PreBComplex_6qx9 -m ../data/MAPT_${MutFile_Name}.tsv -1 92 -2 101 -p NoDMSdata -g MAPT -a Median -c 93
# ### Unfold 10 exonic bases and 17 intronic bases for each sequence (WT and mutations) and calculate the difference in delta G -> from splicesome assembly model B complex 5o9z
# Invivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_BComplex_5o9z -m ../data/MAPT_${MutFile_Name}.tsv -1 84 -2 110 -p DMSdata -g MAPT -a Median -c 93
# Exvivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_BComplex_5o9z -m ../data/MAPT_${MutFile_Name}.tsv -1 84 -2 110 -p DMSdata -g MAPT -a Median -c 93
# No data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_BComplex_5o9z -m ../data/MAPT_${MutFile_Name}.tsv -1 84 -2 110 -p NoDMSdata -g MAPT -a Median -c 93

# ### Unfold 12 exonic bases and 23 intronic bases for each sequence (WT and mutations) and calculate the difference in delta G -> from splicesome assembly model B complex 6ahd
# Invivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_BComplex_6ahd -m ../data/MAPT_${MutFile_Name}.tsv -1 82 -2 116 -p DMSdata -g MAPT -a Median -c 93
# Exvivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_BComplex_6ahd -m ../data/MAPT_${MutFile_Name}.tsv -1 82 -2 116 -p DMSdata -g MAPT -a Median -c 93
# No data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_BComplex_6ahd -m ../data/MAPT_${MutFile_Name}.tsv -1 82 -2 116 -p NoDMSdata -g MAPT -a Median -c 93

# ### Unfold 9 exonic bases and 20 intronic bases for each sequence (WT and mutations) and calculate the difference in delta G -> from splicesome assembly model Pre-Bact complex 7abf
# Invivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_PreBactComplex_7abf -m ../data/MAPT_${MutFile_Name}.tsv -1 85 -2 113 -p DMSdata -g MAPT -a Median -c 93
# Exvivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_PreBactComplex_7abf -m ../data/MAPT_${MutFile_Name}.tsv -1 85 -2 113 -p DMSdata -g MAPT -a Median -c 93
# No data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_PreBactComplex_7abf -m ../data/MAPT_${MutFile_Name}.tsv -1 85 -2 113 -p NoDMSdata -g MAPT -a Median -c 93

# ### Unfold 10 exonic bases and 20 intronic bases for each sequence (WT and mutations) and calculate the difference in delta G -> from splicesome assembly model Bact complex 6ff4
# Invivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_BactComplex_6ff4 -m ../data/MAPT_${MutFile_Name}.tsv -1 84 -2 113 -p DMSdata -g MAPT -a Median -c 93
# Exvivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_BactComplex_6ff4 -m ../data/MAPT_${MutFile_Name}.tsv -1 84 -2 113 -p DMSdata -g MAPT -a Median -c 93
# No data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_BactComplex_6ff4 -m ../data/MAPT_${MutFile_Name}.tsv -1 84 -2 113 -p NoDMSdata -g MAPT -a Median -c 93

# ### Unfold 12 exonic bases and 31 intronic bases for each sequence (WT and mutations) and calculate the difference in delta G -> from splicesome assembly model Bact complex 5Z56
# Invivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_BactComplex_5z56 -m ../data/MAPT_${MutFile_Name}.tsv -1 82 -2 124 -p DMSdata -g MAPT -a Median -c 93
# Exvivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_BactComplex_5z56 -m ../data/MAPT_${MutFile_Name}.tsv -1 82 -2 124 -p DMSdata -g MAPT -a Median -c 93
# No data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t ${IsMut} -b ${TypeOfFold} -r Unfold_LengthRNAinSpliceosome_BactComplex_5z56 -m ../data/MAPT_${MutFile_Name}.tsv -1 82 -2 124 -p NoDMSdata -g MAPT -a Median -c 93

#python calculateDiffDeltaGforUnfoldingBetweenWTandNotWT_SplicesomeAssembly.py ${TypeOfFold} ${MutFile_Name}