import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-r","--variantBedfile",type=str,help="Input a variant file to read in bed format that needs to be converted")
parser.add_argument("-w","--variantWritefile",type=str,help="Input a variant write file that will be written into")
parser.add_argument("-c","--boundarycoordinate",type=int,help="Input chromosome coordinate that indicates the start of the intron")
parser.add_argument("-s","--varIDStart",nargs='?',default=1,type=int,help="Input at what ID number you want to start at, the default will be 1")
parser.add_argument("-p","--prefix",type=str,help="Indicate what prefix you want for your ID for variant")

args = parser.parse_args()

VARIANT_READ_FILE=args.variantBedfile
VARIANT_WRITE_FILE=args.variantWritefile
EXON_INTRON_BOUNDARY=args.boundarycoordinate
STARTID=args.varIDStart
PREFIX=args.prefix

with open(VARIANT_READ_FILE) as f:
    variants_to_read = [line.strip().split("\t") for line in f]
    
to_write=[]

var_count=STARTID
for var in variants_to_read:
    coord=int(var[2])
    if coord < EXON_INTRON_BOUNDARY:
        posToWrite=coord-EXON_INTRON_BOUNDARY
    else:
        posToWrite=coord-EXON_INTRON_BOUNDARY+1
    to_write.append([PREFIX+str(var_count),str(posToWrite),var[3].split(">")[0],var[3].split(">")[1]])
    var_count += 1
    
with open(VARIANT_WRITE_FILE,"w") as fw:
    for i in to_write:
        fw.write("\t".join(i))
        fw.write("\n")
    
