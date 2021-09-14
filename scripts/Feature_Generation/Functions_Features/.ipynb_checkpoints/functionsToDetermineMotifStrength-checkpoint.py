import pandas as pd
import os

def getMotifScoresForSequence(seq,motifPWMfile):
    #Open up motif PWM file, it should be arranged as rows positions and columns ACGTs
    with open(motifPWMfile) as f:
        pwm = [line.strip().split("\t") for line in f]
    # So pwm is now a list of lists where each list is for each position and within each list 
    # are the values for A C G and T
    
    # Store length of motif
    motif_length = len(pwm)
    
    # We create a dictionary where each key is a position of the motif and each value for key is 
    # another dictionary where the key is the nucleotide and the value is the score
    pwm_dict = {}
    for pos_num in range(motif_length):
        ACTG_pos = pwm[pos_num]
        pwm_dict[pos_num] = {}
        pwm_dict[pos_num]['A'] = float(ACTG_pos[0])
        pwm_dict[pos_num]['C'] = float(ACTG_pos[1])
        pwm_dict[pos_num]['G'] = float(ACTG_pos[2])
        pwm_dict[pos_num]['U'] = float(ACTG_pos[3])
    
    # List will contain coordinates and sequence and score 
    validmotifs = []
    
    # Now go through window of motif length through the given sequence and calculate a score 
    # If score above the threshold return that window
    for j in range(0,len(seq)-motif_length+1):
        seqToCheck = seq[j:j+motif_length]
        score_seq = 0
        for i in range(motif_length):
            score_seq += pwm_dict[i][seqToCheck[i]]
        #print score_seq
        #if score_seq > threshold:
        validmotifs.append([str(j),str(j+motif_length),seqToCheck,str(score_seq)])
    
    return(validmotifs)

# Write the motif scores 
# Exon10_ESE, 9
def writeSREMotifScoresForSequence(folder,foldertitle,mutID,motiftype,numclusters,seq):
    with open(folder+"MotifScores_"+foldertitle+"/"+mutID+"_"+motiftype+"_MotifScores.tsv","w") as fw:
        for clust in range(1,numclusters+1):
            motifs_clust = getMotifScoresForSequence(seq,folder+"PWMsForSREs/"+motiftype+"_Cluster"+str(clust)+"_PWM.txt")
            for m in motifs_clust:
                fw.write("Cluster"+str(clust)+"\t")
                fw.write("\t".join(m))
                fw.write("\n")

# Get the difference between WT motif scores and MUT motif scores for a given cluster
# This is specific to a motif type and cluster 
def getSumOfMotifScoreDiffsPerSRECluster(folder,foldertitle,mutID,motifType,clust_num,strength_threshold_dict):
    
    with open(folder+"MotifScores_"+foldertitle+"/"+"WT_"+motifType+"_MotifScores.tsv") as f:
        WT_motifs = [line.strip().split("\t") for line in f]
    # Get mutation motifs 
    with open(folder+"MotifScores_"+foldertitle+"/"+mutID+"_"+motifType+"_MotifScores.tsv") as f:
        MUT_motifs = [line.strip().split("\t") for line in f]
    
    diff_for_cluster=[]
    
    for cluster in range(1,clust_num+1):
        if strength_threshold_dict[motifType+"_Cluster"+str(cluster)] > 0:
            cluster_threshold=strength_threshold_dict[motifType+"_Cluster"+str(cluster)]
        else:
            cluster_threshold=0
        #print cluster_threshold
                
        rel_WT_motifs=[float(i[4]) for i in WT_motifs if i[0]=="Cluster"+str(cluster) and float(i[4])>cluster_threshold]
        rel_MUT_motifs=[float(i[4]) for i in MUT_motifs if i[0]=="Cluster"+str(cluster) and float(i[4])>cluster_threshold]
        
        diff=sum(rel_MUT_motifs)-sum(rel_WT_motifs)
        
        diff_for_cluster.append(str(diff))

    return diff_for_cluster

# Get the difference between WT motif scores and MUT motif scores for a given cluster
# This is specific to a motif type and cluster 
def getRatioOfMotifScoreDiffsPerSRECluster(folder,foldertitle,mutID,motifType,clust_num,strength_threshold_dict):
    
    with open(folder+"MotifScores_"+foldertitle+"/"+"WT_"+motifType+"_MotifScores.tsv") as f:
        WT_motifs = [line.strip().split("\t") for line in f]
    # Get mutation motifs 
    with open(folder+"MotifScores_"+foldertitle+"/"+mutID+"_"+motifType+"_MotifScores.tsv") as f:
        MUT_motifs = [line.strip().split("\t") for line in f]
    
    ratio_for_cluster=[]
    
    for cluster in range(1,clust_num+1):
        if strength_threshold_dict[motifType+"_Cluster"+str(cluster)] > 0:
            cluster_threshold=strength_threshold_dict[motifType+"_Cluster"+str(cluster)]
        else:
            cluster_threshold=0
        #print cluster_threshold
                
        rel_WT_motifs=[float(i[4]) for i in WT_motifs if i[0]=="Cluster"+str(cluster) and float(i[4])>cluster_threshold]
        rel_MUT_motifs=[float(i[4]) for i in MUT_motifs if i[0]=="Cluster"+str(cluster) and float(i[4])>cluster_threshold]
        
        ratio=float(sum(rel_MUT_motifs)+100)/(sum(rel_WT_motifs)+100)
        
        ratio_for_cluster.append(str(ratio))

    return ratio_for_cluster

def createSREclusterThresholdDictionary(folder,dict_NumCluster,whichquantile):
    # Build a dictionary that contains the threshold strength values for each motif and cluster
    strength_threshold_dict={}
    for motifType in ["ESE","ESS","ISE","ISS"]:
        for cluster in range(1,dict_NumCluster[motifType]+1):
            strength_dist = pd.read_csv(folder+"PWMsForSREs/"+motifType+"_Cluster"+str(cluster)+"_kmers_strengths.tsv",header=None)
            strength_dist.columns=["strengths"]
            #print(strength_dist.head())
            cluster_quantile=strength_dist.quantile(whichquantile)[0]
            #print(cluster_quantile)
            strength_threshold_dict[motifType+"_Cluster"+str(cluster)]=cluster_quantile
            
    return(strength_threshold_dict)

"""
def writeRBPMotifScoresForSequence(folder,foldertitle,mutID,seq,region):
    with open(folder+"MotifScores_"+foldertitle+"/"+mutID+"_"+region+"_RBPs_Compendium_MotifScores.tsv","w") as fw:
        for filename in os.listdir(folder+"PWMs_RBPs_Compendium/"):
            if ".txt" in filename:
                motifs_for_PWM = getMotifScoresForSequence(seq,folder+"PWMs_RBPs_Compendium/"+filename)
                for m in motifs_for_PWM:
                    fw.write(filename.split("_")[0]+"\t")
                    fw.write("\t".join(m))
                    fw.write("\n")
"""

def writeRBPMotifScoresForSequence(mutID,seq):
    with open("../tmp/MotifScores_Muts/"+mutID+"_MAPT_RBPs_MotifScores.tsv","w") as fw:
        for filename in os.listdir("../data/MAPT_RBPs/"):
            if ".txt" in filename:
                motifs_for_PWM = getMotifScoresForSequence(seq,"../data/MAPT_RBPs/"+filename)
                for m in motifs_for_PWM:
                    fw.write(filename.split("_")[0]+"\t")
                    fw.write("\t".join(m))
                    fw.write("\n")

# Get the difference between WT motif scores and MUT motif scores for a given RBP motif
"""
def getSumOfMotifScoreDiffsPerRBPmotif(folder,foldertitle,mutID,region,RBPmotifs,strength_threshold_dict):
    
    with open(folder+"MotifScores_"+foldertitle+"/"+"WT_"+region+"_RBPs_Compendium_MotifScores.tsv") as f:
        WT_motifs = [line.strip().split("\t") for line in f]
    # Get mutation motifs 
    with open(folder+"MotifScores_"+foldertitle+"/"+mutID+"_"+region+"_RBPs_Compendium_MotifScores.tsv") as f:
        MUT_motifs = [line.strip().split("\t") for line in f]
    
    diff_for_cluster=[]
    
    for motif in RBPmotifs:
        if strength_threshold_dict[motif] > 0:
            cluster_threshold=strength_threshold_dict[motif]
        else:
            cluster_threshold=0
            #print cluster_threshold
                
        rel_WT_motifs=[float(i[4]) for i in WT_motifs if i[0]==motif and float(i[4])>cluster_threshold]
        rel_MUT_motifs=[float(i[4]) for i in MUT_motifs if i[0]==motif and float(i[4])>cluster_threshold]
        
        diff=sum(rel_MUT_motifs)-sum(rel_WT_motifs)
        
        diff_for_cluster.append(str(diff))

    return diff_for_cluster
"""
def getSumOfMotifScoreDiffsPerRBPmotif(mutID,RBPmotifs,strength_threshold_dict):
    
    with open("../tmp/MotifScores_Muts/WT_MAPT_RBPs_MotifScores.tsv") as f:
        WT_motifs = [line.strip().split("\t") for line in f]
    # Get mutation motifs 
    with open("../tmp/MotifScores_Muts/"+mutID+"_MAPT_RBPs_MotifScores.tsv") as f:
        MUT_motifs = [line.strip().split("\t") for line in f]
    
    diff_for_cluster=[]
    
    for motif in RBPmotifs:
        if strength_threshold_dict[motif] > 0:
            cluster_threshold=strength_threshold_dict[motif]
        else:
            cluster_threshold=0
            #print cluster_threshold
                
        rel_WT_motifs=[float(i[4]) for i in WT_motifs if i[0]==motif and float(i[4])>cluster_threshold]
        rel_MUT_motifs=[float(i[4]) for i in MUT_motifs if i[0]==motif and float(i[4])>cluster_threshold]
        
        diff=sum(rel_MUT_motifs)-sum(rel_WT_motifs)
        
        diff_for_cluster.append(str(diff))

    return diff_for_cluster
# Get the difference between WT motif scores and MUT motif scores for a given RBP motif
def getRatioOfMotifScoreDiffsPerRBPmotif(folder,foldertitle,mutID,region,RBPmotifs,strength_threshold_dict):
    
    with open(folder+"MotifScores_"+foldertitle+"/"+"WT_"+region+"_RBPs_Compendium_MotifScores.tsv") as f:
        WT_motifs = [line.strip().split("\t") for line in f]
    # Get mutation motifs 
    with open(folder+"MotifScores_"+foldertitle+"/"+mutID+"_"+region+"_RBPs_Compendium_MotifScores.tsv") as f:
        MUT_motifs = [line.strip().split("\t") for line in f]
    
    ratio_for_cluster=[]
    
    for motif in RBPmotifs:
        if strength_threshold_dict[motif] > 0:
            cluster_threshold=strength_threshold_dict[motif]
        else:
            cluster_threshold=0
            #print cluster_threshold
                
        rel_WT_motifs=[float(i[4]) for i in WT_motifs if i[0]==motif and float(i[4])>cluster_threshold]
        rel_MUT_motifs=[float(i[4]) for i in MUT_motifs if i[0]==motif and float(i[4])>cluster_threshold]
        
        ratio=float(sum(rel_MUT_motifs)+100)/(sum(rel_WT_motifs)+100)
        
        ratio_for_cluster.append(str(ratio))

    return ratio_for_cluster

"""
def createRBPclusterThresholdDictionary(folder,whichquantile):
    # Build a dictionary that contains the threshold strength values for each motif and cluster
    strength_threshold_dict={}
    for filename in os.listdir(folder+"PWMs_RBPs_Compendium/"):
        if ".txt" in filename:
            strength_dist = pd.read_csv(folder+"PWMs_RBPs_Compendium/"+filename.split("_")[0]+"_kmers_strengths.tsv",header=None)
            strength_dist.columns=["strengths"]
            #print(strength_dist.head())
            cluster_quantile=strength_dist.quantile(whichquantile)[0]
            #print(cluster_quantile)
            strength_threshold_dict[filename.split("_")[0]]=cluster_quantile
            
    return(strength_threshold_dict)
"""
def createRBPclusterThresholdDictionary(whichquantile):
    # Build a dictionary that contains the threshold strength values for each motif and cluster
    strength_threshold_dict={}
    for filename in os.listdir("../data/MAPT_RBPs/"):
        if ".txt" in filename:
            strength_dist = pd.read_csv("../data/MAPT_RBPs/"+filename.split("_")[0]+"_kmers_strengths.tsv",header=None)
            strength_dist.columns=["strengths"]
            #print(strength_dist.head())
            cluster_quantile=strength_dist.quantile(whichquantile)[0]
            #print(cluster_quantile)
            strength_threshold_dict[filename.split("_")[0]]=cluster_quantile
            
    return(strength_threshold_dict)