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
    with open("tmp", "w") as f:
        sys.stdout = f # redirect output to file (or just suppress? https://stackoverflow.com/questions/8447185/to-prevent-a-function-from-printing-in-the-batch-console-in-python)
        fit = sbpy_phase_fit.phase_fit(mpc_number=mpc_number,end_date=end_date,**function_flags)
        check=fit.calculate()
        sys.stdout = sys.__stdout__ # need to reset the redirect!
    return check

def phase_fit_func_name(name,end_date):
    with open("tmp", "w") as f:
        sys.stdout = f
        fit = sbpy_phase_fit.phase_fit(name=name,end_date=end_date,**function_flags)
        check=fit.calculate()
        sys.stdout = sys.__stdout__
    return check

# # We could just use a single function if we are smart passing name and mpc_number...
# def phase_fit_func(mpc_number=False,name=False,end_date):
#     with open("tmp", "w") as f:
#         sys.stdout = f # redirect output to file (or just suppress? https://stackoverflow.com/questions/8447185/to-prevent-a-function-from-printing-in-the-batch-console-in-python)
#         fit = sbpy_phase_fit.phase_fit(mpc_number=mpc_number,name=name,end_date=end_date,**function_flags)
#         check=fit.calculate()
#         sys.stdout = sys.__stdout__ # need to reset the redirect!
#     return check

# wipe the log file
log_file="sbpy_phase_fit.log"
with open(log_file,"w") as f:
    f.close()

# connect to the database
cnx=atlas_database_connection.database_connection().connect()

# ADD A QUERY TO CALL update_atlas_objects
# qry="CALL update_atlas_objects"
# self.cnx.cursor()
# self.cursor.execute(qry)
# self.cnx.commit()
# with open(log_file,"a") as f:
#     f.write("{} start update_atlas_objects\n".format(Time(Time.now(),scale='utc',format="iso")))
# # update_atlas_objects(cnx)
# with open(log_file,"a") as f:
#     f.write("{} start update_atlas_objects\n".format(Time(Time.now(),scale='utc',format="iso")))
# exit()

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
# object_lists=[mpc_number_list[:4],name_list[:4]]
phase_fit_functions=[phase_fit_func_mpc,phase_fit_func_name]

# create a file that records the date, the start and end mpc numbers in a batch, the length of time to complete a batch, the number of jobs in a batch and the number of objects done per sec
today=datetime.now()
# today=datetime.fromisoformat(end_date.value)
fname="fit_phase_curves_bulk_record_{}-{}-{}.csv".format(today.day,today.month,today.year)
print(fname)
of=open(fname,"w")
of.write("date,mpc1,mpc2,t(s),N_jobs,rate(1/s)\n")

# LOG ALL SETTINGS HERE
fit = sbpy_phase_fit.phase_fit(mpc_number=1,name="Ceres",end_date=end_date_mjd,**function_flags)
# git_status = subprocess.check_output("git show --oneline -s",shell=True).decode()
git_status = subprocess.check_output("git log -1",shell=True).decode()
print(git_status)
conda_env = subprocess.check_output("conda env export",shell=True).decode()
print(conda_env)
with open(log_file,"a") as f:
    f.write("date run:\n{}\n".format(today)) # date being run
    f.write("\nlen(mpc_number_list):{}\nlen(name_list){}\n".format(len(mpc_number_list),len(name_list))) # number of objects to be processed
    f.write("\nscript, function, flags:\n{}\n".format(os.path.basename(__file__))) # name of script run
    f.write("{}\n{}\n".format(fit,fit.__dict__)) # flags for the functions
    f.write("\nplatform:\n{}\n".format(platform.uname())) # OS version
    f.write("\ngit log -1:\n{}\n".format(git_status)) # current git branch
    f.write("\nconda env export:\n{}\n".format(conda_env)) # conda environment and packages
    f.write("\n")

for obj_list,phase_fit_func in zip(object_lists,phase_fit_functions):

    print(phase_fit_func.__name__)

    # set up pool and fit phase
    threads=multiprocessing.cpu_count()

    jobs_done=0
    while len(obj_list)>0:

        t_start=time.time()
        sub_list=obj_list[:threads]
        obj_list=obj_list[threads:]

        print("run parallel, {} threads".format(threads))
        print("fit objects {} - {}".format(sub_list[0],sub_list[-1]))
        pool = Pool(threads)
        multiple_results = [pool.apply_async(phase_fit_func, args=(n,end_date_mjd,)) for n in sub_list]
        pool.close()
        pool.join()
        t_end=time.time()
        print()
        jobs_done+=len(sub_list)
        t_elapsed=t_end-t_start

        print("N jobs done = {}".format(jobs_done))
        print("time elapsed = {}".format(t_elapsed))
        print("{} objects/second".format(float(len(sub_list))/t_elapsed))

        print()

        of.write("{},{},{},{},{},{}\n".format(datetime.now(),sub_list[0],sub_list[-1],t_elapsed,len(sub_list),float(len(sub_list))/t_elapsed))
        of.flush()

of.close()

# add clean up tasks?

cmd="mv sbpy_phase_fit.log sbpy_phase_fit_{}-{}-{}.log".format(today.day,today.month,today.year)
subprocess.Popen(cmd,shell=True)
print(cmd)

# download a csv of all the fits
