import argparse
from Bio import pairwise2
import numpy as np
from itertools import combinations

parser = argparse.ArgumentParser()
parser.add_argument("-t","--folderToWrite",type=str,help="Input a folder in which you want to write files into")
parser.add_argument("-m","--motifType",type=str,help="Input one of the of the motif types: ESE, ESS, ISE, ISS")

args = parser.parse_args()

TMPFOLDER=args.folderToWrite
motiftype=args.motifType

import os
if not os.path.exists(TMPFOLDER):
    os.makedirs(TMPFOLDER)
    
# Open up motif file
# Each line has a motif
with open("../../data/"+motiftype+"_Overrepresented.tsv") as f:
    motifs = [line.strip() for line in f]

# Get every combination of motif pairs
motifcombinations = [comb for comb in combinations(motifs, 2)]

# Go through every motif combo -> perform a local alignment to get dissimilarity score
motif_dissimilarity = []
for c in motifcombinations:
    # Align
    alignments = pairwise2.align.localms(c[0],c[1], 5, -4, -2, -0.5)
    
    # Initialize minimum score to Inf
    min_score = 20
    # Go through every alignment and dissimilarity distance defined as the number of shifts plus the number of mismatches in the best local alignment of the two hexamers
    for a in alignments:
        seq1 = a[0]
        seq2 = a[1]
        s = 0
        for i in range(len(seq1)):
            if seq1[i]=='-' or seq2[i]=='-' or seq1[i]!=seq2[i]:
                s += 1
        # If this is the lowest score thus far, make it the new minimum score
        if s < min_score:
            min_score = s
            
    # Append the two motifs and the min dissimilarity score to list     
    motif_dissimilarity.append([c[0],c[1],str(min_score)])

# Write to file    
with open(TMPFOLDER+motiftype+"_DissimilarityScore.tsv",'w') as fw:
    for m in motif_dissimilarity:
        fw.write("\t".join(m))
        fw.write("\n")