#!/bin/bash
# coding: utf-8

# # Steps for generating features 

# ## Feature 1 : Unfolding energy of regions: Splice Site, region around mutation, larger region around splice site
# * Splice site is defined as 3 bases upstream and 6 bases downstream of exon-intron boundary

# ## Ensemble based calculations

# ### Generate structural ensemble for each sequence (WT and mutations) incorporating DMS data collected for WT sequence

# Invivo data
#python generateStructuralEnsemble_NoMaxPairingDist_RNAStructure.py -f ../data/MAPT_exon10intron10withDMSdata.fa -m  ../data/MAPT_SingleMutations_TrainingSet.tsv -b 93 -d ../data/Invivo_Exon10Intron10_PooledRep_RenormalizedPerNucleotide_ScaledGUs_GsAndUsSetToLargeNegativeValue_RewrittenCoordinatesIncludeAllExon10.shape -t /home/jkumar/Projects/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/StructuralEnsemble_Muts/  

# Exvivo data
#python generateStructuralEnsemble_NoMaxPairingDist_RNAStructure.py -f ../data/MAPT_exon10intron10withDMSdata.fa -m  ../data/MAPT_SingleMutations_TrainingSet.tsv -b 93 -d ../data/Exvivo_Exon10Intron10_PooledRep_RenormalizedPerNucleotide_ScaledGUs_GsAndUsSetToLargeNegativeValue_RewrittenCoordinatesIncludeAllExon10.shape -t /home/jkumar/Projects/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/StructuralEnsemble_Muts/

# no data 
#python generateStructuralEnsemble_NoMaxPairingDist_RNAStructure.py -f ../data/MAPT_exon10intron10withDMSdata.fa -m  ../data/MAPT_SingleMutations_TrainingSet.tsv -b 93 -t /home/jkumar/Projects/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/StructuralEnsemble_Muts/


# ### Calculate the median energy of the ensemble for each mutation and WT

# Invivo data
#python calulateMedianDeltaGForEnsemblePerMutation.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -m ../data/MAPT_SingleMutations_TrainingSet.tsv -p Rsample -a Median

# Exvivo data
#python calulateMedianDeltaGForEnsemblePerMutation.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -m ../data/MAPT_SingleMutations_TrainingSet.tsv -p Rsample -a Median

# No data
#python calulateMedianDeltaGForEnsemblePerMutation.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -m ../data/MAPT_SingleMutations_TrainingSet.tsv -p Partition -a Median


# ### Unfold just the splice site for each sequence (WT and mutations) and calculate the difference in delta G 

# Invivo data
#python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b StructuralEnsemble -r Unfold_SpliceSite -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 91 -2 99 -p Rsample -g MAPT -a Median -c 93

# Exvivo data
#python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b StructuralEnsemble -r Unfold_SpliceSite -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 91 -2 99 -p Rsample -g MAPT -a Median -c 93

# No data
#python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b StructuralEnsemble -r Unfold_SpliceSite -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 91 -2 99 -p Partition -g MAPT -a Median -c 93


# ### Unfold the local region around each mutation and calculate the difference in delta G - use window of 5
# Invivo data
#python calculateDeltaGforUnfoldingEnergyAroundMutation.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -c StructuralEnsemble -r Unfold_LocalRegion -m ../data/MAPT_SingleMutations_TrainingSet.tsv -w 5 -b 93 -p Rsample -g MAPT -a Median

# Exvivo data
#python calculateDeltaGforUnfoldingEnergyAroundMutation.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -c StructuralEnsemble -r Unfold_LocalRegion -m ../data/MAPT_SingleMutations_TrainingSet.tsv -w 5 -b 93 -p Rsample -g MAPT -a Median

# No data
#python calculateDeltaGforUnfoldingEnergyAroundMutation.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -c StructuralEnsemble -r Unfold_LocalRegion -m ../data/MAPT_SingleMutations_TrainingSet.tsv -w 5 -b 93 -p Partition -g MAPT -a Median


# ### Unfold 10 bases around the splice site for each sequence (WT and mutations) and calculate the difference in delta G 
# Invivo data
#python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b StructuralEnsemble -r Unfold_20basesAroundSpliceSite -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 84 -2 103 -p Rsample -g MAPT -a Median -c 93

# Exvivo data
#python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b StructuralEnsemble -r Unfold_20basesAroundSpliceSite -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 84 -2 103 -p Rsample -g MAPT -a Median -c 93

# No data
#python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b StructuralEnsemble -r Unfold_20basesAroundSpliceSite -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 84 -2 103 -p Partition -g MAPT -a Median -c 93

# ### Unfold 12 exonic bases and 31 intronic bases for each sequence (WT and mutations) and calculate the difference in delta G -> from splicesome assembly model

# Invivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b StructuralEnsemble -r Unfold_MaxLengthRNAinSpliceosome -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 82 -2 124 -p Rsample -g MAPT -a Median -c 93

# Exvivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b StructuralEnsemble -r Unfold_MaxLengthRNAinSpliceosome -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 82 -2 124 -p Rsample -g MAPT -a Median -c 93

# No data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b StructuralEnsemble -r Unfold_MaxLengthRNAinSpliceosome -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 82 -2 124 -p Partition -g MAPT -a Median -c 93

# ## MFE based calculations
# ### Calculate the energy of the MFE of sequence 
# Invivo data
#python generateMFEstructure_NoMaxPairingDist_RNAStructure.py -f ../data/MAPT_exon10intron10withDMSdata.fa -m  ../data/MAPT_SingleMutations_TrainingSet.tsv -b 93 -d ../data/Invivo_Exon10Intron10_PooledRep_RenormalizedPerNucleotide_ScaledGUs_GsAndUsSetToLargeNegativeValue_RewrittenCoordinatesIncludeAllExon10.shape -t /home/jkumar/Projects/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -q Muts  

# Exvivo data
#python generateMFEstructure_NoMaxPairingDist_RNAStructure.py -f ../data/MAPT_exon10intron10withDMSdata.fa -m  ../data/MAPT_SingleMutations_TrainingSet.tsv -b 93 -d ../data/Exvivo_Exon10Intron10_PooledRep_RenormalizedPerNucleotide_ScaledGUs_GsAndUsSetToLargeNegativeValue_RewrittenCoordinatesIncludeAllExon10.shape -t /home/jkumar/Projects/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -q Muts  

# No data
#python generateMFEstructure_NoMaxPairingDist_RNAStructure.py -f ../data/MAPT_exon10intron10withDMSdata.fa -m  ../data/MAPT_SingleMutations_TrainingSet.tsv -b 93 -t /home/jkumar/Projects/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -q Muts  


# ### Unfold splice site
# Invivo data
#python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b MFEstructures -r Unfold_SpliceSite -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 91 -2 99 -p DMSdata -g MAPT -a Median -c 93

# Exvivo data
#python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b MFEstructures -r Unfold_SpliceSite -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 91 -2 99 -p DMSdata -g MAPT -a Median -c 93

# No data
#python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b MFEstructures -r Unfold_SpliceSite -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 91 -2 99 -p NoDMSdata -g MAPT -a Median -c 93


# ### Unfold the local region around each mutation and calculate the difference in delta G - use window of 5
# Invivo data
#python calculateDeltaGforUnfoldingEnergyAroundMutation.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -c MFEstructures -r Unfold_LocalRegion -m ../data/MAPT_SingleMutations_TrainingSet.tsv -w 5 -b 93 -p DMSdata -g MAPT -a Median

# Exvivo data
#python calculateDeltaGforUnfoldingEnergyAroundMutation.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -c MFEstructures -r Unfold_LocalRegion -m ../data/MAPT_SingleMutations_TrainingSet.tsv -w 5 -b 93 -p DMSdata -g MAPT -a Median

# No data
#python calculateDeltaGforUnfoldingEnergyAroundMutation.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -c MFEstructures -r Unfold_LocalRegion -m ../data/MAPT_SingleMutations_TrainingSet.tsv -w 5 -b 93 -p NoDMSdata -g MAPT -a Median


# ### Unfold 10 bases around the splice site for each sequence (WT and mutations) and calculate the difference in delta G 
# Invivo data
#python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b MFEstructures -r Unfold_20basesAroundSpliceSite -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 84 -2 103 -p DMSdata -g MAPT -a Median -c 93

# Exvivo data
#python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b MFEstructures -r Unfold_20basesAroundSpliceSite -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 84 -2 103 -p DMSdata -g MAPT -a Median -c 93

# No data
#python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b MFEstructures -r Unfold_20basesAroundSpliceSite -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 84 -2 103 -p NoDMSdata -g MAPT -a Median -c 93


# ### Unfold 12 exonic bases and 31 intronic bases for each sequence (WT and mutations) and calculate the difference in delta G -> from splicesome assembly model
# Invivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b MFEstructures -r Unfold_MaxLengthRNAinSpliceosome -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 82 -2 124 -p DMSdata -g MAPT -a Median -c 93

# Exvivo data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/ExvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b MFEstructures -r Unfold_MaxLengthRNAinSpliceosome -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 82 -2 124 -p DMSdata -g MAPT -a Median -c 93
# No data
python calculateDeltaGforUnfoldingCoordinatePair.py -f /home/jkumar/Projects/Model_MAPTsplicing/tmp/NoDMSdata_PooledRep_Exon10Intron10_Structures/ -t Muts -b MFEstructures -r Unfold_MaxLengthRNAinSpliceosome -m ../data/MAPT_SingleMutations_TrainingSet.tsv -1 82 -2 124 -p NoDMSdata -g MAPT -a Median -c 93

# ### Calculate the difference in delta G between WT and MUT for all regions interested

# Calculation for invivo, exvivo and nodata all done in a single script
python calculateDiffDeltaGforUnfoldingBetweenWTandNotWT.py StructuralEnsemble SingleMutations_TrainingSet
python calculateDiffDeltaGforUnfoldingBetweenWTandNotWT.py MFEstructures SingleMutations_TrainingSet

