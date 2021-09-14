from itertools import product
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t","--folderToWrite",type=str,help="Input a folder in which you want to write files into")
parser.add_argument("-m","--motifType",type=str,help="Input one of the of the motif types: ESE, ESS, ISE, ISS")
parser.add_argument("-n","--numberClusters",type=int,help="Input the number of clusters for the motif type")
args = parser.parse_args()

TMPfolder=args.folderToWrite
motif=args.motifType
numclusters=args.numberClusters

all_scores=[]

for clust in range(1,numclusters+1):
    with open(TMPfolder+motif+"_Cluster"+str(clust)+"_PWM.txt") as f:
        pwm = [line.strip().split("\t") for line in f]
    # So pwm is now a list of lists where each list is for each position and within each list 
    # are the values for A C G and U
    
    # Store length of motif
    motif_length = len(pwm)
    print(motif_length)
    
    all_kmers=[]
    for i in product('ACGU',repeat=motif_length):
        all_kmers.append("".join(i))
    print(len(all_kmers))
    # We create a dictionary where each key is a position of the motif and each value for key is 
    # another dictionary where the key is the nucleotide and the value is the score
    pwm_dict = {}
    for pos_num in range(motif_length):
        ACGU_pos = pwm[pos_num]
        pwm_dict[pos_num] = {}
        pwm_dict[pos_num]['A'] = float(ACGU_pos[0])
        pwm_dict[pos_num]['C'] = float(ACGU_pos[1])
        pwm_dict[pos_num]['G'] = float(ACGU_pos[2])
        pwm_dict[pos_num]['U'] = float(ACGU_pos[3])
       
    scores_p=[]
    for seq in all_kmers:
        score_seq=0
        for i in range(motif_length):
            score_seq += pwm_dict[i][seq[i]]
        scores_p.append(score_seq)
    
    print(len(scores_p))
    
    with open(TMPfolder+motif+"_Cluster"+str(clust)+"_kmers_strengths.tsv","w") as fw:
        for i in scores_p:
            fw.write(str(i)+"\n")