from calculate_phase import sbpy_phase_fit
import pandas as pd
import mysql.connector
import numpy as np
import time
from datetime import datetime
import multiprocessing
from multiprocessing import Pool
import subprocess
import os
import cProfile
import sys
from contextlib import redirect_stdout
from calculate_phase import atlas_database_connection
from calculate_phase.atlas_SQL_query_df import update_atlas_objects
from astropy.time import Time
import platform

# Set the required flags for the phase fit function
function_flags = {"push_fit_flag":True,"hide_warning_flag":False,"mag_diff_flag":True}
# function_flags = {"push_fit_flag":False,"hide_warning_flag":False,"mag_diff_flag":True}
print(function_flags)
# print(*function_flags)

# set the date beyond which no obs are counted
# end_date=Time(Time.now(),scale='utc',format="iso")
end_date=Time("2021-02-14",scale='utc',format="iso")
end_date_mjd=round(end_date.mjd)
print(end_date,end_date_mjd)

# define the fitting functions
def phase_fit_func_mpc(mpc_number,end_date):
    fit = sbpy_phase_fit.phase_fit(mpc_number=mpc_number,end_date=end_date,**function_flags)
    check=fit.calculate()
    return check

def phase_fit_func_name(name,end_date):
    fit = sbpy_phase_fit.phase_fit(name=name,end_date=end_date,**function_flags)
    check=fit.calculate()
    return check

# wipe the log file
log_file="sbpy_phase_fit.log"
with open(log_file,"w") as f:
    f.close()

# connect to the database
cnx=atlas_database_connection.database_connection().connect()

# perform the query and store results as a dataframe
qry="select mpc_number,name from atlas_objects;"
# qry="select mpc_number,name from atlas_objects limit 10;"
df=pd.read_sql_query(qry,cnx)
print(df)
cnx.close()

# make the mpc number list
df_mpc=df[~np.isnan(df["mpc_number"])]
mpc_number_list=df_mpc['mpc_number'].astype(int)
mpc_number_list=list(mpc_number_list)
print(mpc_number_list[:10])
print(len(mpc_number_list))

# make the name list
df_name=df[np.isnan(df["mpc_number"])]
name_list=list(df_name["name"])
print(name_list[:10])
print(len(name_list))

# do all the mpc_number fits, then do all the name fits
object_lists=[mpc_number_list,name_list]
object_lists=[mpc_number_list[:4],name_list[:4]]
phase_fit_functions=[phase_fit_func_mpc,phase_fit_func_name]

for mpc_number in mpc_number_list:
    print(mpc_number)
    phase_fit_func_mpc(mpc_number,end_date)
    break

for name in name_list:
    print(name)
    phase_fit_func_name(name,end_date)
    break
