import sys
import subprocess
import numpy as np

def convertCTFile(ctfile,coordsOfInterest,genename,ctfile_unfolded):
    
    # Open CT File and read lines
    with open(ctfile) as f:
        lines = [line.strip() for line in f]
    
    #First line
    first_line=lines[0].split(" ")
    first_line_values = [i for i in first_line if i!=""]
    length_seq=int(first_line_values[0])
    
    # These are the coordinates of interest in string form
    coords = [str(i) for i in range(coordsOfInterest[0],coordsOfInterest[1]+1)]
    
    # To write list
    to_write = []
    
    #Keep count of lines changing 
    lines_changed_count = 0
    
    for line in lines:
        # Let's split line by " " and then remove all the " "
        line_split = line.split(" ")
        line_values = [k for k in line_split if k!=""]
        
        # lines which dont't contain the header line AND contain the coordinates of interest
        if len([i for i in line_values if genename in i])<1 and (line_values[0] in coords or line_values[len(line_values)-2] in coords):
        # then replace the third from the last with zero string
            line_values[len(line_values)-2] = '0'
            lines_changed_count += 1
        
        # Get rid of all the extra spaces and then join remaining values by a space
        to_write.append(" ".join(line_values))
    
    print(lines_changed_count)
    
    # Write converted file
    with open(ctfile_unfolded,"w") as fw:
        for line in to_write:
            fw.write(line)
            fw.write("\n")

def getEnergyOfUnfoldedRegion(efn2file,efn2file_unfolded):
    
    # Calculate differences between original and unfolded
    with open(efn2file) as f:
        prior = [float(line.strip().split(' ')[6]) for line in f]
    with open(efn2file_unfolded) as f:
        post = [float(line.strip().split(' ')[6]) for line in f]
    
    energyToUnfold = list(np.array(post)-np.array(prior))
    
    return(energyToUnfold)