# This is a script to convert the dissimilarity table into a matrix 
#import pandas as pd
import argparse
import logging
import numpy as np
from scipy.cluster import hierarchy
import matplotlib
import getClustersFromDendogram as gcd
import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument("-t","--folderToWrite",type=str,help="Input a folder in which you want to write files into")
parser.add_argument("-m","--motifType",type=str,help="Input one of the of the motif types: ESE, ESS, ISE, ISS")
parser.add_argument("-l","--thresholdScore",type=float,help="Input a cut threshold to indicate where in the dendogram we want to make the cut")
args = parser.parse_args()

TMPFOLDER=args.folderToWrite
motiftype=args.motifType
threshold=args.thresholdScore

# Open up motif file
with open("../../data/"+motiftype+"_Overrepresented.tsv") as f:
    motifs = [line.strip() for line in f]

# Create dictionary for number to motif and vice-versa
motifToNum_dict = {motifs[i]:i for i in range(len(motifs))}
#numToMotif_dict = {i:motifs[i] for i in range(len(motifs))}

# Open up the dissimilarity score table file
with open(TMPFOLDER+motiftype+"_DissimilarityScore.tsv") as f:
    score_motifs = [line.strip().split("\t") for line in f]

# Convert dissimilarity table to a matrix, remember i,j entry same as j,i entry
dissimilarityScore_matrix = np.zeros((len(motifs),len(motifs)))
for pair in score_motifs:
    r_pair = motifToNum_dict[pair[0]]
    c_pair = motifToNum_dict[pair[1]]
    score = int(pair[2])
    dissimilarityScore_matrix[r_pair,c_pair]=score
    dissimilarityScore_matrix[c_pair,r_pair]=score

#Size of dissimilarity matrix
z = dissimilarityScore_matrix.shape[0]
# pdist list
pdist = []
# Get the flat array pdist 
for i in range(z):
    for j in range(i+1,z):
        pdist.append(dissimilarityScore_matrix[i,j])

print(len(pdist))

print(np.array(pdist))

matplotlib.use('Agg')

# This has to be called after matplotlib.use('Agg')
import matplotlib.pyplot as plt

Z = hierarchy.linkage(pdist, 'average')

#print Z
fig = plt.figure(figsize=(25, 10))

colors_to_use=["#ca6f66","#0eb690","#8dc57d","#edb953","#ee86be","#647fc5","#7a437d","#3a7cb4","#be2b84","#e86615","#df5cfc","#32b3ac","#9e18f3","#b0a499","#ec0e3b"]
hierarchy.set_link_color_palette(colors_to_use)
#dn = dendrogram(Z,color_threshold=3.47,orientation='left',labels=motifs)
dn = hierarchy.dendrogram(Z,color_threshold=threshold,orientation='left',labels=motifs)

#print(dn['color_list'])
#print(dn['ivl'])

plt.vlines(x=threshold,ymin=0,ymax=3000)

plt.savefig(TMPFOLDER+motiftype+'_dendogram.svg', format='svg', bbox_inches='tight',)

plt.show()

#from scipy.cluster.hierarchy import fcluster
#clusters = fcluster(Z, 3.47, criterion='distance')
#print(clusters)

os.chdir(TMPFOLDER)
print(os.getcwd())
cluster_classes = gcd.get_cluster_classes(dn)

k = 1 

print(cluster_classes)

# For every cluster, re-align the sequences using clustalw2 to get an alignment
for key in cluster_classes.keys():
    if len(cluster_classes[key]) > 4:
        
        with open(motiftype+"_Cluster"+str(k)+".fa","w") as fw:
            for item in list(set(cluster_classes[key])):
                #print item
                fw.write(">"+item+"\n")
                fw.write(item+"\n")
        
        par1 = "-INFILE="+motiftype+"_Cluster"+str(k)+".fa"
        print(par1)
        par2 = "-OUTFILE="+motiftype+"_Cluster"+str(k)+"_Alignment.fa"
        print(par2)
        subprocess.check_output(["clustalw2",par1,par2,"-OUTPUT=FASTA"])
        k += 1

os.chdir("/home/jkumar/Projects/Model_MAPTsplicing/scripts")
print(os.getcwd())