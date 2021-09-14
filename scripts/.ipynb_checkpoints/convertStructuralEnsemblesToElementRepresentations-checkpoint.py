import argparse
import Functions_Features.functionsToGetMutationSeqAndDMSfile as fmd

parser = argparse.ArgumentParser()
parser.add_argument("-s","--seed", nargs='?',default=50,type=int,help="Input the seed for generating sample structures using stochastic")
parser.add_argument("-t","--folderToWrite",type=str,help="Input a folder in which you want to write structural ensemble files into")
parser.add_argument("-d","--dmsornot",type=str,help="Indicate if Rsample or Partition output")
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")
args = parser.parse_args()

SEED=args.seed
TMP_FOLDER=args.folderToWrite
DMSORNOT=args.dmsornot
MUTATION_FILE=args.mutationfile

# Open up the mutation file which contains the mutation IDs
with open(MUTATION_FILE) as f:
    mutations = [line.strip().split("\t") for line in f]

for mut in mutations:
    mutID=mut[0]
    fmd.convertEnsembleToElementRepresentation(TMP_FOLDER,mutID,str(SEED),DMSORNOT)