import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-f","--folderToWrite",type=str,help="Input a folder in which you want to write files into and access structural ensemble data, make sure to include / at the end")
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")
parser.add_argument("-t","--folderTitle",type=str,help="Input the title of the mutation file, aka is it MUTS or SNPS")
parser.add_argument("-p","--dmsornot",type=str,help="Indicate if Rsample or Partition")
parser.add_argument("-s","--seed", nargs='?',default=50,type=int,help="Input the seed for generating sample structures using stochastic, this is optional, default is 50")
parser.add_argument("-a","--statistic",type=str,help="Input statistic, Mean or Median to calculate")

args = parser.parse_args()

TMPfolder=args.folderToWrite
MUTATION_FILE=args.mutationfile
FOLDERTITLE=args.folderTitle
DMSORNOT=args.dmsornot
SEED=str(args.seed)
STATISTIC=args.statistic

# Open up the mutation file which contains the position of mutations 
with open(MUTATION_FILE) as f:
    mutations = [line.strip().split("\t") for line in f]

to_write=[]
# Go through every mutation, run parition, followed by stochastic to generate ensemble of 1000, next convert the dot bracket structures to the element representation
for mut in mutations:
    mutID = mut[0]
    print(mutID)
    
    efn2filename = TMPfolder+"StructuralEnsemble_"+FOLDERTITLE+"/"+mutID+"_"+DMSORNOT+"Output_NoMaxDist_Seed"+str(SEED)+"_efn2.txt"
    
    with open(efn2filename) as f:
        enfn_vals = [float(line.strip().split(' ')[6]) for line in f]
    
    if STATISTIC=="Mean":
        val = np.mean(enfn_vals)
    else:
        val = np.median(enfn_vals)
        
    to_write.append([mutID,str(val)])

with open(TMPfolder+MUTATION_FILE.split("/")[2].split(".")[0]+"_"+STATISTIC+"DeltaGforEnsemble.tsv","w") as fw:
    for i in to_write:
        fw.write("\t".join(i))
        fw.write("\n")  