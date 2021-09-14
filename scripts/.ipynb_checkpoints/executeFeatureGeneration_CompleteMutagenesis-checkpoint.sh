#!/bin/bash

#python calculate5pSpliceSiteStrength.py -t /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/ -w ../data/MAPT_exon10intron10withDMSdata.fa -b 93 -m ../data/MAPT_CompleteMutagenesis.tsv

#python getLocationAndStrengthOfSREbasedOnPWM.py -w /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/ -t CompleteMutagenesis -f ../data/MAPT_exon10intron10withDMSdata.fa -b 93 -m ../data/MAPT_CompleteMutagenesis.tsv

#python calculateSREstrengthDifferenceBetweenWTandMUT_perCluster.py -w /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/ -t CompleteMutagenesis -m ../data/MAPT_CompleteMutagenesis.tsv

#python calculateDiffInStrengthforMotifsBetweenWTandNotWT.py CompleteMutagenesis

#python generateStructuralEnsemble_NoMaxPairingDist_RNAStructure.py -f ../data/MAPT_exon10intron10withDMSdata.fa -m  ../data/MAPT_CompleteMutagenesis.tsv -b 93 -d ../data/Invivo_Exon10Intron10_PooledRep_RenormalizedPerNucleotide_ScaledGUs_GsAndUsSetToLargeNegativeValue_RewrittenCoordinatesIncludeAllExon10.shape -t /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/StructuralEnsemble_CompleteMutagenesis/

#python calulateMedianDeltaGForEnsemblePerMutation.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t CompleteMutagenesis -m ../data/MAPT_CompleteMutagenesis.tsv -p Rsample -a Median

#python calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t CompleteMutagenesis -b StructuralEnsemble -r Unfold_LengthRNAinSpliceosome_BactComplex_6ff4 -m ../data/MAPT_CompleteMutagenesis.tsv -1 84 -2 113 -p Rsample -g MAPT -a Median -c 93


python3 calculateDeltaGforUnfoldingCoordinatePair.py -f /data/data2/Projects/jkumar/Model_MAPTsplicing/tmp/InvivoDMSdata_PooledRep_Exon10Intron10_Structures/ -t CompleteMutagenesis -b StructuralEnsemble -r Unfold_LengthRNAinSpliceosome_BactComplex_5z56 -m /data/data2/Projects/jkumar/Model_MAPTsplicing/data/MAPT_CompleteMutagenesis.tsv -1 82 -2 124 -p Rsample -g MAPT -a Median -c 93