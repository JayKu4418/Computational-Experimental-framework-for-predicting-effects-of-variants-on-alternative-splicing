import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s","--shapefile",type=str,help="Input the shape file for rewriting")
parser.add_argument("-c","--newcoordinatestart",type=int,help="Indicate the coordinate at which existing data should start from")

args = parser.parse_args()

shapeFile=args.shapefile
coord=args.newcoordinatestart

with open(shapeFile) as f:
    data=[line.strip().split("\t") for line in f]
    
# New data 
newdata_toWrite=[]

# For additional new coordinates, shape val will be -999.0
for i in range(coord):
    newdata_toWrite.append([str(i+1),"-999.0"])

# For existing data, the coordinates need to be modified to reflect new coordinates
for j in data:
    newcoord = int(j[0])+coord
    newdata_toWrite.append([str(newcoord),j[1]])

writefile=shapeFile.split(".shape")[0] + "_RewrittenCoordinates.shape"

with open(writefile,"w") as fw:
    for i in newdata_toWrite:
        fw.write("\t".join(i))
        fw.write("\n")