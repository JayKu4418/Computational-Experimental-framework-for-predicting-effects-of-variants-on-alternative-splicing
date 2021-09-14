import sys
import Functions_Features.functionsToDetermineMotifStructure as fds

typeOfData=sys.argv[1]
dmsornot=sys.argv[2]

TMP_FOLDER="/home/jkumar/Projects/Model_MAPTsplicing/tmp/"+typeOfData+"data_PooledRep_Exon10Intron10_Structures/"

average_energy_PreB = fds.getAverageDeltaGUnfoldingForCoords(TMP_FOLDER,"StructuralEnsemble","Muts","SplicesomeAssemblyStages_WT/PreB_","WT",92,101,"MAPT",50,"Median",dmsornot)

average_energy_B_5o9z = fds.getAverageDeltaGUnfoldingForCoords(TMP_FOLDER,"StructuralEnsemble","Muts","SplicesomeAssemblyStages_WT/B_5o9z_","WT",84,110,"MAPT",50,"Median",dmsornot)

average_energy_B_6ahd = fds.getAverageDeltaGUnfoldingForCoords(TMP_FOLDER,"StructuralEnsemble","Muts","SplicesomeAssemblyStages_WT/B_6ahd_","WT",82,116,"MAPT",50,"Median",dmsornot)

average_energy_PreBact = fds.getAverageDeltaGUnfoldingForCoords(TMP_FOLDER,"StructuralEnsemble","Muts","SplicesomeAssemblyStages_WT/PreBact_","WT",85,113,"MAPT",50,"Median",dmsornot)

average_energy_Bact_6ff4 = fds.getAverageDeltaGUnfoldingForCoords(TMP_FOLDER,"StructuralEnsemble","Muts","SplicesomeAssemblyStages_WT/Bact_6ff4_","WT",84,113,"MAPT",50,"Median",dmsornot)

average_energy_Bact_5z56 = fds.getAverageDeltaGUnfoldingForCoords(TMP_FOLDER,"StructuralEnsemble","Muts","SplicesomeAssemblyStages_WT/Bact_5z56","WT",82,124,"MAPT",50,"Median",dmsornot)


print([average_energy_PreB,average_energy_B_5o9z,average_energy_B_6ahd,average_energy_PreBact,average_energy_Bact_6ff4,average_energy_Bact_5z56])


