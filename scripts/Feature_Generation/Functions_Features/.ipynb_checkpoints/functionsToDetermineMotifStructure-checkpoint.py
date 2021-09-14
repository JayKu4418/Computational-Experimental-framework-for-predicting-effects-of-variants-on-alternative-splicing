import Functions_Features.functionsToUnfoldAndCalculateEnergyOfRegionOfInterest as fcu
import subprocess
import numpy as np

# This function will get the median delta G of unfolding for a given start and end site
def getAverageDeltaGUnfoldingForCoords(folder,foldingType,MutType,subMotiffolder,mutID,start,end,genename,seed,statistic,dmsornot=""):
    print(dmsornot)
    if foldingType=="StructuralEnsemble":
        middle_part =dmsornot+"Output_NoMaxDist_Seed"+str(seed)
    else:
        middle_part ="MFEstructure_NoMaxDist" + ("_" + dmsornot if dmsornot else "")
        
    fcu.convertCTFile(folder+foldingType+"_"+MutType+"/"+mutID+"_"+middle_part+".ct",[start,end],genename,folder+subMotiffolder+mutID+"_"+middle_part+"_UnfoldedMotif.ct")
    
    subprocess.check_output(["efn2-smp",folder+subMotiffolder+mutID+"_"+middle_part+"_UnfoldedMotif.ct", folder+subMotiffolder+mutID+"_"+middle_part+"_UnfoldedMotif_efn2.txt"])
    
    energy_diff = fcu.getEnergyOfUnfoldedRegion(folder+foldingType+"_"+MutType+"/"+mutID+"_"+middle_part+"_efn2.txt",folder+subMotiffolder+mutID+"_"+middle_part+"_UnfoldedMotif_efn2.txt")
    
    with open(folder+subMotiffolder+mutID+"_"+middle_part+"_UnfoldedMotif_EnergyDifference.txt","w") as fw:
        for j in energy_diff:
            fw.write(str(j))
            fw.write("\n")
    
    if statistic=="Mean":
        return(np.mean(energy_diff))
    else:
        return(np.median(energy_diff))

# This function will get the median delta G of unfolding if two regions need to be unfolded -> this is specifically for double mutations
def getAverageDeltaGUnfoldingForTwoSetsOfCoords(folder,foldingType,MutType,subMotiffolder,mutID,start1,end1,start2,end2,genename,seed,statistic,dmsornot=""):
    print(dmsornot)
    if foldingType=="StructuralEnsemble":
        middle_part =dmsornot+"Output_NoMaxDist_Seed"+str(seed)
    else:
        middle_part ="MFEstructure_NoMaxDist" + ("_" + dmsornot if dmsornot else "")
        
    fcu.convertCTFile(folder+foldingType+"_"+MutType+"/"+mutID+"_"+middle_part+".ct",[start1,end1],genename,folder+subMotiffolder+mutID+"_"+middle_part+"_UnfoldedMotif_Part1.ct")
    fcu.convertCTFile(folder+subMotiffolder+mutID+"_"+middle_part+"_UnfoldedMotif_Part1.ct",[start2,end2],genename,folder+subMotiffolder+mutID+"_"+middle_part+"_UnfoldedMotif_Part2.ct")
    
    subprocess.check_output(["efn2-smp",folder+subMotiffolder+mutID+"_"+middle_part+"_UnfoldedMotif_Part2.ct", folder+subMotiffolder+mutID+"_"+middle_part+"_UnfoldedMotif_efn2.txt"])
    
    energy_diff = fcu.getEnergyOfUnfoldedRegion(folder+foldingType+"_"+MutType+"/"+mutID+"_"+middle_part+"_efn2.txt",folder+subMotiffolder+mutID+"_"+middle_part+"_UnfoldedMotif_efn2.txt")
    
    with open(folder+subMotiffolder+mutID+"_"+middle_part+"_UnfoldedMotif_EnergyDifference.txt","w") as fw:
        for j in energy_diff:
            fw.write(str(j))
            fw.write("\n")
    
    if statistic=="Mean":
        return(np.mean(energy_diff))
    else:
        return(np.median(energy_diff))

    
# This function will get the coordinates that meet a specific threshold for a motif type and cluster
# The coordinates will be unique
def getCoordsForSREMotifs(folder,foldertitle,mutID,motifType,clust_num,strength_threshold_dict,boundary):
    # Get mutation motifs 
    with open(folder+"MotifScores_"+foldertitle+"/"+mutID+"_"+motifType+"_MotifScores.tsv") as f:
        MUT_motifs = [line.strip().split("\t") for line in f]
    
    # This is the empty vector that will contain all the unique coordinates to unfold for that motif type
    coordsToUnfold=[]
    
    for cluster in range(1,clust_num+1):
        
        #print strength_threshold_dict[motifName+"_Cluster"+str(cluster)]
        if strength_threshold_dict[motifType+"_Cluster"+str(cluster)] > 0:
            cluster_threshold=strength_threshold_dict[motifType+"_Cluster"+str(cluster)]
        else:
            cluster_threshold=0
        #print cluster_threshold
        # we only want to grab coords whoese 
        rel_MUT_motifs=[i for i in MUT_motifs if i[0]=="Cluster"+str(cluster) and float(i[4])>cluster_threshold]
        # Go through each of relative motif locations and if the coordinates have not already been added to the list, add them
        for i in rel_MUT_motifs:
            if [i[1],i[2]] not in coordsToUnfold:
                if motifType in ["ISE","ISS"]:
                    coordsToUnfold.append([str(int(i[1])+boundary),str(int(i[2])+boundary)])
                else:
                    coordsToUnfold.append([i[1],i[2]])
        
    return(coordsToUnfold)

def getMedianDeltaGUnfoldingPerSRECluster(folder,foldertitle,mutID,motifType,clust_num,strength_threshold_dict,deltaGforMut,boundary):
    
    # Only have boundary value be valid if intronic motif type
    if motifType in ["ESE","ESS"]:
        boundary=0
    
    with open(folder+"MotifScores_"+foldertitle+"/"+mutID+"_"+motifType+"_MotifScores.tsv") as f:
        MUT_motifs = [line.strip().split("\t") for line in f]
    
    # This is the empty vector that will contain all the median deltaGs for each cluster
    deltaG_Clusters=[]
    
    for cluster in range(1,clust_num+1):
        
        #print strength_threshold_dict[motifName+"_Cluster"+str(cluster)]
        if strength_threshold_dict[motifType+"_Cluster"+str(cluster)] > 0:
            cluster_threshold=strength_threshold_dict[motifType+"_Cluster"+str(cluster)]
        else:
            cluster_threshold=0
        #print cluster_threshold
        # we only want to grab coords whoese 
        rel_MUT_motif_Coords=[[str(int(i[1])+boundary),str(int(i[2])+boundary)] for i in MUT_motifs if i[0]=="Cluster"+str(cluster) and float(i[4])>cluster_threshold]
        
        coord_DeltaG = [float(i[3]) for i in deltaGforMut if [i[1],i[2]] in rel_MUT_motif_Coords]
        
        deltaG_Clusters.append(str(np.median(coord_DeltaG)))
        
    return(deltaG_Clusters)
        

# This function will get the coordinates that meet a specific threshold for a motif type and cluster
# The coordinates will be unique
def getCoordsForRBPMotifs(folder,foldertitle,mutID,RBPmotifs,region,strength_threshold_dict,boundary):
    
    with open(folder+"MotifScores_"+foldertitle+"/"+mutID+"_"+region+"_RBPs_Compendium_MotifScores.tsv") as f:
        MUT_motifs = [line.strip().split("\t") for line in f]
    
    # This is the empty vector that will contain all the unique coordinates to unfold for that motif type
    coordsToUnfold=[]
    
    for motif in RBPmotifs:
        if strength_threshold_dict[motif] > 0:
            cluster_threshold=strength_threshold_dict[motif]
        else:
            cluster_threshold=0
            #print cluster_threshold
                
        rel_MUT_motifs=[i for i in MUT_motifs if i[0]==motif and float(i[4])>cluster_threshold]
        
        # Go through each of relative motif locations and if the coordinates have not already been added to the list, add them
        for i in rel_MUT_motifs:
            if [i[1],i[2]] not in coordsToUnfold:
                if region=="Intron":
                    coordsToUnfold.append([str(int(i[1])+boundary),str(int(i[2])+boundary)])
                else:
                    coordsToUnfold.append([i[1],i[2]])
        
    return(coordsToUnfold)


def getMedianDeltaGUnfoldingPerPerRBP(folder,foldertitle,mutID,RBPmotifs,region,strength_threshold_dict,deltaGforMut,boundary):
    
    # Only have boundary value be valid if intronic motif type
    if region=="Exon":
        boundary=0
    
    with open(folder+"MotifScores_"+foldertitle+"/"+mutID+"_"+region+"_RBPs_Compendium_MotifScores.tsv") as f:
        MUT_motifs = [line.strip().split("\t") for line in f]
    
    # This is the empty vector that will contain all the median deltaGs for each RBP motif
    deltaG_RBPs=[]
    
    for motif in RBPmotifs:
        if strength_threshold_dict[motif] > 0:
            cluster_threshold=strength_threshold_dict[motif]
        else:
            cluster_threshold=0
        
        # we only want to grab coords whoese 
        rel_MUT_motif_Coords=[[str(int(i[1])+boundary),str(int(i[2])+boundary)] for i in MUT_motifs if i[0]==motif and float(i[4])>cluster_threshold]
        
        if len(rel_MUT_motif_Coords)>0:
            coord_DeltaG = [float(i[3]) for i in deltaGforMut if [i[1],i[2]] in rel_MUT_motif_Coords]
        
            deltaG_RBPs.append(str(np.median(coord_DeltaG)))
        else:
            deltaG_RBPs.append(str(0.00))
        
    return(deltaG_RBPs)