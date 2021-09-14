# Remember for seqlogo you need to run it using python3
import seqlogo
import argparse
import random
import sys 
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-t","--folderToWrite",type=str,help="Input a folder in which you want to write files into")
parser.add_argument("-m","--motifType",type=str,help="Input one of the of the motif types: ESE, ESS, ISE, ISS")
parser.add_argument("-n","--numberClusters",type=int,help="Input the number of clusters for the motif type")
parser.add_argument("-p","--pseudocount",nargs='?',default=1.0,type=float,help="Input pseudocount for calculating PWM, the default value is 1")


args = parser.parse_args()

TMPfolder=args.folderToWrite
motifType=args.motifType
numClusters=args.numberClusters
pseudocount=args.pseudocount

# Format a float number to have 3 decimal places
def format(value):
    return "%.3f" % value

# Go through every cluster and grab the alignment ils 
for clust in range(1,numClusters+1):
    #print(clust)
    with open(TMPfolder+motifType+"Clusters/"+motifType+"_Cluster"+str(clust)+"_Alignment.fa") as f:
        lines = [line.strip() for line in f]
    motifs = [i for i in lines if ">" not in i]
    #print(motifs)
    motif_length = max([len(i) for i in motifs])
    num_motifs = len(motifs)
    # Get the number of counts of A,C,T and G 
    # This will create a matrix of 4 rows where each row is a nucleotide and columns of length of motif length
    counts_ACTG = np.zeros((4,motif_length))
    #print(motif_length)
    for i in range(motif_length):
        char_at_i = [j[i] for j in motifs]
        counts_ACTG[0,i] = char_at_i.count("A")
        counts_ACTG[1,i] = char_at_i.count("C")
        counts_ACTG[2,i] = char_at_i.count("G")
        counts_ACTG[3,i] = char_at_i.count("T")
    
    #Empty array for cumulative counts
    counts_ACTG_cum = np.zeros((4,motif_length))
    
    # Generate random counts for gap values a 100 times and add to counts
    # This is in order to randomly distribute the gaps between the 4 nucleotides 
    # Otherwise we are left with instances of gaps but it is unclear where to put these
    for z in range(100):
        rand_counts = np.zeros((4,motif_length))
        for i in range(motif_length):
            char_at_i = [j[i] for j in motifs]
            num_gaps = char_at_i.count("-")
            if num_gaps == 0:
                rand_counts_i = [0,0,0,0]
            else:
                a = random.randint(0, num_gaps)
                b = random.randint(0, num_gaps-a)
                c = random.randint(0, num_gaps-(a+b))
                rand_counts_i = [a, b, c, num_gaps - (a+b+c)]
            rand_counts[:,i] = rand_counts_i
        # The matrix rand_counts will account for the gaps randomly distributed across nucleotides
        # We are going to randomly do this a 100 times and keep adding random counts 
        counts_ACTG_toAdd = (counts_ACTG + rand_counts + (pseudocount/4))/(num_motifs+pseudocount)
        counts_ACTG_cum += counts_ACTG_toAdd
        #counts_ACTG_ToVisualize = counts_ACTG + rand_counts
    
    # Divide the cumulative counts by 100 since we did this a 100 times
    counts_ACTG_withpseudocount = counts_ACTG_cum/100
    #print(np.sum(counts_ACTG_withpseudocount,axis=0))
    # Assume a background nucleotide rate for each nucleoide as 0.25
    PWM_cluster = np.log2(counts_ACTG_withpseudocount/0.25)
    #print(PWM_cluster)
    with open(TMPfolder+motifType+"_Cluster"+str(clust)+"_PWM.txt","w") as fw:
        for line in PWM_cluster.transpose():
            for i in line:
                fw.write(str(format(i)))
                fw.write("\t")
            fw.write("\n")
    PPM_seqLogo = seqlogo.Ppm(counts_ACTG_withpseudocount.transpose())
    #print(counts_ACTG_ToVisualize.transpose())
    #print(seqlogo.pfm2pwm(counts_ACTG_ToVisualize.transpose(), background = 0.25, pseudocount = 0.8))
    
    #print(clust)
    seqlogo.seqlogo(PPM_seqLogo, ic_scale = False, filename=TMPfolder+motifType+"Clusters/"+motifType+"_Cluster"+str(clust)+"_fromPWM.png", format = 'png', size = 'medium')