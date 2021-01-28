from calculate_phase import sbpy_phase_fit
import pandas as pd
import mysql.connector
import numpy as np
import time
import datetime
import multiprocessing
from multiprocessing import Pool
import subprocess
import os
import cProfile
import sys
from contextlib import redirect_stdout
from calculate_phase import atlas_database_connection

# define the fitting functions
def phase_fit_func_mpc(mpc_number):
    with open("tmp", "w") as f:
        sys.stdout = f # redirect output to file (or just suppress? https://stackoverflow.com/questions/8447185/to-prevent-a-function-from-printing-in-the-batch-console-in-python)
        fit = sbpy_phase_fit.phase_fit(mpc_number=mpc_number,push_fit_flag=True,hide_warning_flag=True)
        check=fit.calculate()
        sys.stdout = sys.__stdout__ # need to reset the redirect!
    # fit = sbpy_phase_fit.phase_fit(mpc_number=mpc_number,push_fit_flag=True,hide_warning_flag=True)
    # check=fit.calculate()

    return check

def phase_fit_func_name(name):
    with open("tmp", "w") as f:
        sys.stdout = f
        fit = sbpy_phase_fit.phase_fit(name=name,push_fit_flag=True,hide_warning_flag=True)
        check=fit.calculate()
        sys.stdout = sys.__stdout__
    return check

# connect to the database
cnx=atlas_database_connection.database_connection().connect()
# perform the query and store results as a dataframe
qry="select mpc_number,name from atlas_objects;"
df=pd.read_sql_query(qry,cnx)
# df=df.dropna()
print(df)
cnx.close()

# make the mpc number list
#mpc_number_list=np.random.choice(df['mpc_number'].astype(int),10,replace=False)
df_mpc=df[~np.isnan(df["mpc_number"])]
mpc_number_list=df_mpc['mpc_number'].astype(int)
# print(mpc_number_list)
# print(len(mpc_number_list))

# mpc_number_list=mpc_number_list[:100]

mpc_number_list=list(mpc_number_list)
print(mpc_number_list[:10])
print(len(mpc_number_list))

# make the name list
df_name=df[np.isnan(df["mpc_number"])]
name_list=list(df_name["name"])
print(name_list[:10])
print(len(name_list))

# do all the mpc_number fits, then do all the name fits
# object_lists=[mpc_number_list,name_list]
object_lists=[mpc_number_list[:4],name_list[:4]]
phase_fit_functions=[phase_fit_func_mpc,phase_fit_func_name]

# create a file that records the date, the start and end mpc numbers in a batch, the length of time to complete a batch, the number of jobs in a batch and the number of objects done per sec
today=datetime.datetime.now()
fname="fit_phase_curves_bulk_record_{}-{}-{}.csv".format(today.day,today.month,today.year)
print(fname)
of=open(fname,"w")
of.write("date,mpc1,mpc2,t(s),N_jobs,rate(1/s)\n")

for obj_list,phase_fit_func in zip(object_lists,phase_fit_functions):
    # print(obj_list[:10])
    print(phase_fit_func.__name__)

    # quick test
    # for j in range(len(obj_list)):
        # print("start {}".format(obj_list[j]))
        # check=phase_fit_func(obj_list[j])

    # print("run in serial")
    # t_start=time.time()
    # count=0
    # for n in mpc_number_list:
    #     phase_fit_func(n)
    #     # if count>3:
    #     #     break
    #     # else:
    #     #     count+=1
    #     count+=1
    # t_end=time.time()
    # print ("N_jobs = {}".format(count))
    # print ("t_pool = {}".format(t_end-t_start))
    # print("{} objects/second".format(float(count)/(t_end-t_start)))
    # exit()

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
        multiple_results = [pool.apply_async(phase_fit_func, args=(n,)) for n in sub_list]
        pool.close()
        pool.join()
        t_end=time.time()
        print()
        jobs_done+=len(sub_list)
        t_elapsed=t_end-t_start

        # # try check if the fit worked
        # for i,r in enumerate(multiple_results):
        #     try:
        #         print(r.get())
        #     except:
        #         print("error {}".format(sub_list[i]))

        print("N jobs done = {}".format(jobs_done))
        print("time elapsed = {}".format(t_elapsed))
        print("{} objects/second".format(float(len(sub_list))/t_elapsed))

        print()

        of.write("{},{},{},{},{},{}\n".format(datetime.datetime.now(),sub_list[0],sub_list[-1],t_elapsed,len(sub_list),float(len(sub_list))/t_elapsed))
        of.flush()

of.close()

# add clean up tasks?

cmd="mv sbpy_phase_fit.log sbpy_phase_fit_{}-{}-{}.log".format(today.day,today.month,today.year)
subprocess.Popen(cmd,shell=True)
print(cmd)

# download a csv of all the fits
