import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-m","--mapfile",type=str,help="Input the shape file for rewriting")
parser.add_argument("-s","--numNucleotidesAddToStart",type=int,help="Indicate number of nucleotides to add to start, can be zero")
parser.add_argument("-e","--numNucleotidesAddToEnd",type=int,help="Indicate number of nucleotides to add to end, can be zero")
parser.add_argument("-a","--nametoappend",type=str,help="Name to append at end of file")

args = parser.parse_args()

mapFile=args.mapfile
coordStart=args.numNucleotidesAddToStart
coordEnd=args.numNucleotidesAddToEnd
nameToAppend=args.nametoappend

with open(mapFile) as f:
    data=[line.strip().split("\t") for line in f]
    
# New data 
newdata_toWrite=[]

# For additional new coordinates, shape val will be -999.0
for i in range(coordStart):
    newdata_toWrite.append([str(i+1),"-999.0","0.0","N"])

# For existing data, the coordinates need to be modified to reflect new coordinates
for j in data:
    newcoord = int(j[0])+coordStart
    newdata_toWrite.append([str(newcoord),j[1],j[2],j[3]])

# Last coordinate thus far
lastcoord=int(newdata_toWrite[len(newdata_toWrite)-1][0])
print(lastcoord)

for i in range(lastcoord,lastcoord+coordEnd):
    newdata_toWrite.append([str(i+1),"-999.0","0.0","N"])

writefile=mapFile.split(".map")[0] + "_"+nameToAppend+".map"

with open(writefile,"w") as fw:
    for i in newdata_toWrite:
        fw.write("\t".join(i))
        fw.write("\n")
        
writefile_shape=mapFile.split(".map")[0] + "_"+nameToAppend+".shape"

with open(writefile_shape,"w") as fw:
    for i in newdata_toWrite:
        fw.write("\t".join(i[0:2]))
        fw.write("\n")