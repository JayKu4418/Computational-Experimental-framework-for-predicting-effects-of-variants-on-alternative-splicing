import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-m","--mapfile",type=str,help="Input the normalized map file to get Gs and Us scaled")
parser.add_argument("-w","--writefilename",type=str,help="Input the write file name to write the new file")

args = parser.parse_args()

mapFile=args.mapfile
namefileToWrite=args.writefilename

mapdata = pd.read_csv(mapFile,sep="\t",header=None)
mapdata.columns = ['Position', 'DMS', 'StdErr', 'Nucleotide']

#print(mapdata.head())

mapdata.loc[(mapdata["DMS"] < -500), "StdErr"] = np.NaN
mapdata.loc[(mapdata["DMS"] < -500), "DMS"] = np.NaN

#print(mapdata.head())

justG = mapdata[mapdata["Nucleotide"]=="G"]["DMS"]
justG_scaled = justG*(0.1/np.nanmax(justG))

justU = mapdata[mapdata["Nucleotide"]=="U"]["DMS"]
justU_scaled = justU*(0.1/np.nanmax(justU))


mapdata.loc[(mapdata["Nucleotide"]=="G"), "DMS"] = justG_scaled
mapdata.loc[(mapdata["Nucleotide"]=="U"), "DMS"] = justU_scaled

#print(mapdata.head())

mapdata.to_csv(namefileToWrite + ".map",sep="\t",header=False,index=False,na_rep='NA')

