import argparse
import Functions_Features.functionsToDetermineMotifStrength as fdm
import Functions_Features.functionsToGetMutationSeqAndDMSfile as fmd

parser = argparse.ArgumentParser()
parser.add_argument("-t","--folderToWrite",type=str,help="Input a folder in which you want to write files into")
parser.add_argument("-w","--wtfastafile", help="Input the wildtype fasta file",type=str)
parser.add_argument("-b","--intronexonboundary", help="Input the intron exon boundary from which mutation position is calculated",type=int)
parser.add_argument("-m","--mutationfile",type=str,help="Input a mutation file")

args = parser.parse_args()

TMPfolder=args.folderToWrite
WTFASTAFile=args.wtfastafile
INTRON_EXON_BOUNDARY=args.intronexonboundary
MUTATION_FILE=args.mutationfile

with open(WTFASTAFile) as f:
    WTseq=[line.strip() for line in f]
WTseq=WTseq[1]

with open(MUTATION_FILE) as f:
#with open("../data/MAPT_MUTs_ToTest.tsv") as f:
    mutations=[line.strip().split("\t") for line in f]

for mut in mutations:
    mutID = mut[0]
    # Skip the first line which is for the wildtype
    if mutID=="WT":
        seqToWrite=WTseq
    # Position at last base of exon is -1 and first base of intron is 1
    else:
        if len(mut)>5:
            data=fmd.obtainSequenceAndDMSdataBasedOnDoubleMutation(mut,INTRON_EXON_BOUNDARY,WTseq)
        else:
            data=fmd.obtainSequenceAndDMSdataBasedOnSingleMutation(mut,INTRON_EXON_BOUNDARY,WTseq)
    
        seqToWrite=data[0]
    
    newboundary=INTRON_EXON_BOUNDARY
    # Only applicable to insertions and deletions in exon
    # Insertion because WTbase is NA 
    # Need to increase start site and end site by number of bases inserted 
    if mutID!="WT" and int(mut[1]) < 0 and mut[2]=="NA":
        newboundary=INTRON_EXON_BOUNDARY+len(mut[3])
    # Deletion because Mutbase is NA
    # Need to decrease the start site and end site by number of bases deleted
    if mutID!="WT" and int(mut[1]) < 0 and mut[3]=="NA":
        newboundary=INTRON_EXON_BOUNDARY-len(mut[2])
    
    exonseqToAnalyze=seqToWrite[:newboundary]
    intronseqToAnalyze=seqToWrite[newboundary:]
    
   
    fdm.writeRBPMotifScoresForSequence(TMPfolder,mutID,exonseqToAnalyze,"Exon")
    fdm.writeRBPMotifScoresForSequence(TMPfolder,mutID,intronseqToAnalyze,"Intron")