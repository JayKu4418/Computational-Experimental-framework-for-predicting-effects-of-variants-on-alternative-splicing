import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p","--profilefile",type=str,help="Input the profile file for renormalization")
parser.add_argument("-w","--writefilename",type=str,help="Input the write file name without any extension to write renormalized data for map and shape files")

args = parser.parse_args()

profileFile=args.profilefile
namefileToWrite=args.writefilename

# Function to normalize reactivities
# Nans have been removed from reactivities
def return_normalizingReactivitiesFactor(reactivities):
    IQR_lim=scipy.stats.iqr(reactivities)
    perc_lim=np.percentile(reactivities, 90)
    threshold_set = max(IQR_lim,perc_lim)
    reactivities_afterThreshold = reactivities[reactivities<threshold_set]
    topTen_reactivities_afterThreshold = reactivities_afterThreshold[reactivities_afterThreshold>np.percentile(reactivities_afterThreshold,90)]
    return(np.mean(topTen_reactivities_afterThreshold))


# Read in profile file
profile_data = pd.read_csv(profileFile,sep="\t",header=0)

# Assign blank columns for new normalized reactivities and std errors
profile_data = profile_data.assign(myNewNormalizedReactivityBy=profile_data["HQ_profile"])
profile_data = profile_data.assign(myNewNormalizedStdErrBy=profile_data["HQ_stderr"])


# The same factor found for reactivities is used to normalize the standard errors. DO NOT NEED TO calculate a separate normalizing factor for std errors
for nucleotide in ["A","C"]:
    
    print(nucleotide)
    
    reactivities_nucleotide = profile_data[profile_data["Sequence"]==nucleotide]["HQ_profile"].dropna().values
    factor = return_normalizingReactivitiesFactor(reactivities_nucleotide)
    #print(factor)
    
    normalized_reactivities = profile_data[profile_data["Sequence"]==nucleotide]["HQ_profile"]/factor
    profile_data.loc[profile_data["Sequence"]==nucleotide,"myNewNormalizedReactivityBy"] = normalized_reactivities
    
    normalized_stderrs = profile_data[profile_data["Sequence"]==nucleotide]["HQ_stderr"]/factor
    profile_data.loc[profile_data["Sequence"]==nucleotide,"myNewNormalizedStdErrBy"] = normalized_stderrs
    
# For nucleotides G and U, we choose to ignore the untreated rate and just normalize on the modified rate, since modification rates are quite low compared to As and Cs
for nucleotide in ["G","U"]:
    
    print(nucleotide)
    
    reactivities_nucleotide_profile = profile_data[(profile_data["Sequence"]==nucleotide)]
    reactivities_nucleotide_profile.loc[reactivities_nucleotide_profile["HQ_profile"].isna(),"Modified_rate"] = np.NaN
    reactivities_nucleotide = reactivities_nucleotide_profile["Modified_rate"].dropna().values
    factor = return_normalizingReactivitiesFactor(reactivities_nucleotide)
    #print(factor)
    
    normalized_reactivities = reactivities_nucleotide_profile["Modified_rate"]/factor
    profile_data.loc[profile_data["Sequence"]==nucleotide,"myNewNormalizedReactivityBy"] = normalized_reactivities
    
    stderrs_nucleotide_profile = profile_data[(profile_data["Sequence"]==nucleotide)]
    stderrs_nucleotide_profile.loc[stderrs_nucleotide_profile["HQ_profile"].isna(),"Modified_rate"] = np.NaN
    stderrs_nucleotide_profile.loc[stderrs_nucleotide_profile["HQ_profile"].isna(),"Modified_effective_depth"] = np.NaN
    stderrs_nucleotide = np.sqrt(stderrs_nucleotide_profile["Modified_rate"])/np.sqrt(stderrs_nucleotide_profile["Modified_effective_depth"])
    normalized_stderrs = stderrs_nucleotide/factor
    profile_data.loc[profile_data["Sequence"]==nucleotide,"myNewNormalizedStdErrBy"] = normalized_stderrs
    
# Write out the reactivities into a map file after setting NA values to -999 and corresponding std errors to 0
mapfiledata_toWrite = profile_data.loc[:,["Nucleotide","myNewNormalizedReactivityBy","myNewNormalizedStdErrBy","Sequence"]]
mapfiledata_toWrite.loc[mapfiledata_toWrite["myNewNormalizedReactivityBy"].isna(),"myNewNormalizedReactivityBy"] = -999.0
mapfiledata_toWrite.loc[mapfiledata_toWrite["myNewNormalizedStdErrBy"].isna(),"myNewNormalizedStdErrBy"] = 0

mapfiledata_toWrite.to_csv(namefileToWrite + ".map",sep="\t",header=False,index=False)
mapfiledata_toWrite.loc[:,["Nucleotide","myNewNormalizedReactivityBy"]].to_csv(namefileToWrite + ".shape",sep="\t",header=False,index=False)