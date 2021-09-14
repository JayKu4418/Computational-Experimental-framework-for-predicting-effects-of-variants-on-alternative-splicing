# Remember for seqlogo you need to run it using python3
import argparse
import numpy as np
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-r","--folderToRead",type=str,help="Input a folder in which you want to read files from")
parser.add_argument("-w","--folderToWrite",type=str,help="Input a folder in which you want to write files into")
parser.add_argument("-i","--IDfile",type=str,help="Input an ID file to identify the motifs and match them to the RBP")

args = parser.parse_args()

TMPfolder=args.folderToRead
WriteFolder=args.folderToWrite
idFile=args.IDfile

# Format a float number to have 3 decimal places
def format(value):
    return "%.3f" % value


#Get IDs that are relevant for Homo sapiens
with open(idFile) as f:
    idlines=[line.strip().split("\t") for line in f][1:]

relevantIDs = [i[0] for i in idlines if i[3]=="Homo_sapiens"]

# Go through every PFM file in the folder 
for filename in os.listdir(TMPfolder):
    if (".txt" in filename) and (filename.split("_")[0] in relevantIDs):
        data=pd.read_csv(TMPfolder+filename,sep="\t",header=0,index_col=0)
        data_matrix = data.values
    
        PWM_cluster = np.log2(data_matrix/0.25)
        writefilename="_".join(filename.split(".txt")[0].split("_")[0:2])+"_pwm.txt"
        #print(PWM_cluster)
        with open(WriteFolder+"PWMs_RBPs_Compendium/"+writefilename,"w") as fw:
            for line in PWM_cluster:
                for i in line:
                    fw.write(str(format(i)))
                    fw.write("\t")
                fw.write("\n")