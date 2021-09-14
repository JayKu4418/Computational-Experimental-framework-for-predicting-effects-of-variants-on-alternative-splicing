#!/bin/python

import sys

# Get the map file
inputfile=str(sys.argv[1])

with open(inputfile) as f:
    data=[line.strip().split("\t") for line in f]
    
# New data 
newdata_toWrite=[]

for i in data:
    if i[1]=="NA" or i[1]=="-999":
        newdata_toWrite.append([i[0],"-999.0"])
    elif i[3]=="A" or i[3]=="C":
        newdata_toWrite.append([i[0],i[1]])
    else:
        newdata_toWrite.append([i[0],"-999.0"])

writefile=inputfile.split(".map")[0] + "_GsAndUsSetToLargeNegativeValue.shape"

with open(writefile,"w") as fw:
    for i in newdata_toWrite:
        fw.write("\t".join(i))
        fw.write("\n")

# New data to write map file
newdata_toWriteMapFile=[]
        
for i in data:
    if i[1]=="NA" or i[1]=="-999":
        newdata_toWriteMapFile.append([i[0],"-999.0","0.0",i[3]])
    elif i[3]=="A" or i[3]=="C":
        newdata_toWriteMapFile.append([i[0],i[1],i[2],i[3]])
    else:
        newdata_toWriteMapFile.append([i[0],"-999.0","0.0",i[3]])
        

writefile_mapFile=inputfile.split(".map")[0] + "_GsAndUsSetToLargeNegativeValue.map"

with open(writefile_mapFile,"w") as fw:
    for i in newdata_toWriteMapFile:
        fw.write("\t".join(i))
        fw.write("\n")