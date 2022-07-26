"""
Script to find the data for missing atlas exposures.
Searches through /atlas/data on the db4 server to retrieve the .ddt or .ddc files.
Run this on the db4 server.
"""

from astropy.time import Time
import pandas as pd
import numpy as np
import sys
import os
import subprocess

def get_data_file(filename):
    """Function to read the data file"""
    with open(filename, "r") as f:
        data = f.readlines()
    data = [d for d in data if d.startswith("#") and "=" in d]
    return data

def data_to_dict(data):
    """Function to turn the data array into a dictionary"""
    dat_dict={}
    for d in data:
        _d = d.split()
        key = _d[1].replace("=","")
        val = _d[-1]
        dat_dict[key] = val
    return dat_dict

# read in list of missing exposures
# df_rA = pd.read_csv("df_rockAtlas_missing_exposures.csv",index_col=0)
df_rA = pd.read_csv("df_rockAtlas_missing_exposures2.csv",index_col=0)

cols = ["OBS","MJD","MAG5SIG"] # columns to get from data files
extension = [".ddt",".ddc"] # possible file extensions

# set up initial file
t = Time.now()
out_file = 'rockAtlas_missing_exposures_{:.3f}.csv'.format(t.mjd)
f = open(out_file, 'w', buffering=1)
f.write(",".join(cols)+"\n")

for i in range(len(df_rA)):

    expname = df_rA.iloc[i]["expname"]
    scope = expname[:3]
    night = expname[3:8]
    file = "/atlas/diff/{}/{}/{}".format(scope,night,expname)
    # file = "data_files/{}".format(expname)

    # try each extension to get data file
    data = []
    for e in extension:
        _file = file+e

        # check file exists
        if os.path.isfile(_file):
            print(i,_file)
            data = get_data_file(_file)

            # make sure there is data in the file, otherwise read the next file
            if len(data)>0:
       	        dat_dict = data_to_dict(data)
                dat = [dat_dict[x] for x in cols]
                f.write(",".join(dat)+"\n")
                break
            else:
                continue
        else:
            print(i,"{} not found".format(_file))

    # if i>10:
    #     break

f.close()
