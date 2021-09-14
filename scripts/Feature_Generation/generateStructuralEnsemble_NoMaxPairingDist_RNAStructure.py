# Script generates a structural ensemble ->CT files for a given sequence and a set of mutations 

#import pandas as pd
import argparse
import logging
import Functions_Features.functionsToGetMutationSeqAndDMSfile as fmd
from joblib import Parallel, delayed
import multiprocessing

parser = argparse.ArgumentParser()
parser.add_argument("-f","--wtfastafile", help="Input the wildtype fasta file",type=str)
parser.add_argument("-b","--intronexonboundary", help="Input the location of the last position of the exon in the wildtype fasta file",type=int)
parser.add_argument("-s","--seed", nargs='?',default=50,type=int,help="Input the seed for generating sample structures using stochastic")
parser.add_argument("-d","--dmsfile",nargs="?",type=str,help="Input a DMS file if you want to use Rsample to generate structures otherwise partition is used")
parser.add_argument("-t","--folderToWrite",type=str,help="Input a folder in which you want to write structural ensemble files into")
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")
args = parser.parse_args()

WT_FASTA_FILE=args.wtfastafile
INTRON_EXON_BOUNDARY=args.intronexonboundary
SEED=args.seed
DMS_FILE = args.dmsfile
TMP_FOLDER=args.folderToWrite
MUTATION_FILE=args.mutationfile

#logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

import os
if not os.path.exists(TMP_FOLDER):
    os.makedirs(TMP_FOLDER)

# Function to get ensemble for a given mutation
def getEnsembleForMutation(mutation):
    
    mutID = mutation[0]
    
    # Skip the first line which is for the wildtype
    if mutID=="WT":
        fmd.generateEnsemble(TMP_FOLDER,"WT",WT_FASTA_FILE,SEED,DMS_FILE)
    # Position at last base of exon is -1 and first base of intron is 1
    else:
        # double mutation
        if len(mut)>5:
            data=fmd.obtainSequenceAndDMSdataBasedOnDoubleMutation(mutation,INTRON_EXON_BOUNDARY,WTseq,wt_dms_data)
        # else single mutation
        else:
            data=fmd.obtainSequenceAndDMSdataBasedOnSingleMutation(mutation,INTRON_EXON_BOUNDARY,WTseq,wt_dms_data)
        
        seqToWrite=data[0]
        dmsvalsToWrite=data[1]
        
        MUT_FASTA_FILE=TMP_FOLDER+mutID+"_MutatedSequence.fa"
        
        # Create mutation fasta file
        with open(MUT_FASTA_FILE,"w") as fw:
            fw.write(">MAPT_withMUT"+"\n")
            fw.write(seqToWrite+"\n")
        
        # If no dms file provided, then no mutated dms file given
        if DMS_FILE is not None:
            MUT_DMS_FILE = TMP_FOLDER+mutID+"_MUT_DMSdata.shape"
            # Create new DMS file
            with open(MUT_DMS_FILE,"w") as fw:
                for dmsval in dmsvalsToWrite:
                    fw.write("\t".join(dmsval))
                    fw.write("\n")
        else:
            MUT_DMS_FILE = ""
        
        fmd.generateEnsemble(TMP_FOLDER,mutID,MUT_FASTA_FILE,SEED,MUT_DMS_FILE)

# Open the wildtype fasta file to get WT seq
with open(WT_FASTA_FILE) as f:
    lines=[line.strip() for line in f]
WTseq=lines[1]

# Open the wildtype DMS file to get DMS data if there is a DMS file
if DMS_FILE:
    with open(DMS_FILE) as f:
        wt_dms_data = [line.strip().split("\t") for line in f]
else:
    wt_dms_data=[]

# Open up the mutation file which contains the position of mutations 
with open(MUTATION_FILE) as f:
    mutations = [line.strip().split("\t") for line in f]

# Get number of cores
#num_cores = multiprocessing.cpu_count()

# Run all the mutations in parallel
# Go through every mutation, run parition, followed by stochastic to generate ensemble of 1000, next convert the dot bracket structures to the element representation
#Parallel(n_jobs=num_cores)(delayed(getEnsembleForMutation)(mut) for mut in mutations)
for mut in mutations:
    getEnsembleForMutation(mut)