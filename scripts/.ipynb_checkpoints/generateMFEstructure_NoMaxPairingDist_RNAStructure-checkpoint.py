#import pandas as pd
import argparse
import logging
import Functions_Features.functionsToGetMutationSeqAndDMSfile as fmd
from joblib import Parallel, delayed
import multiprocessing

parser = argparse.ArgumentParser()
parser.add_argument("-f","--wtfastafile", help="Input the wildtype fasta file",type=str)
parser.add_argument("-b","--intronexonboundary", help="Input the location of the last position of the exon in the wildtype fasta file",type=int)
parser.add_argument("-d","--dmsfile",nargs="?",type=str,help="Input a DMS file if you want to use DMS vals to generate MFE structure")
parser.add_argument("-t","--folderToWrite",type=str,help="Input the main level folder path in which you want to write MFE files into, make sure to add / at end")
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")
parser.add_argument("-q","--folderTitle",type=str,help="Input the title of the mutation file, aka is it MUTS or SNPS")
args = parser.parse_args()

WT_FASTA_FILE=args.wtfastafile
INTRON_EXON_BOUNDARY=args.intronexonboundary
DMS_FILE = args.dmsfile
TMP_FOLDER=args.folderToWrite
MUTATION_FILE=args.mutationfile
FOLDERTITLE=args.folderTitle

folderToWriteInto_MFEs = TMP_FOLDER+"MFEstructures_"+FOLDERTITLE+"/"

#logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
import os
if not os.path.exists(folderToWriteInto_MFEs):
    os.makedirs(folderToWriteInto_MFEs)

# Function to get MFE for a given mutation
def getMFEForMutation(mutation):
    
    mutID = mutation[0]
    
    # Skip the first line which is for the wildtype
    if mutID=="WT":
        mfe_val = fmd.generateMFE(folderToWriteInto_MFEs,"WT",WT_FASTA_FILE,DMS_FILE)
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
        
        MUT_FASTA_FILE=folderToWriteInto_MFEs+mutID+"_MutatedSequence.fa"
        
        # Create mutation fasta file
        with open(MUT_FASTA_FILE,"w") as fw:
            fw.write(">MAPT_withMUT"+"\n")
            fw.write(seqToWrite+"\n")
        
        # If no dms file provided, then no mutated dms file given
        if DMS_FILE is not None:
            MUT_DMS_FILE = folderToWriteInto_MFEs+mutID+"_MUT_DMSdata.shape"
            # Create new DMS file
            with open(MUT_DMS_FILE,"w") as fw:
                for dmsval in dmsvalsToWrite:
                    fw.write("\t".join(dmsval))
                    fw.write("\n")
        else:
            MUT_DMS_FILE = ""
        
        mfe_val = fmd.generateMFE(folderToWriteInto_MFEs,mutID,MUT_FASTA_FILE,MUT_DMS_FILE)
    
    return(mfe_val)

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
to_write = []
for mut in mutations:
    mutID = mut[0]
    val = getMFEForMutation(mut)
    to_write.append([mutID,str(val)])

with open(TMP_FOLDER+MUTATION_FILE.split("/")[2].split(".")[0]+"_"+"DeltaGforMFEstructure.tsv","w") as fw:
    for i in to_write:
        fw.write("\t".join(i))
        fw.write("\n")