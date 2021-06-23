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

def phase_fit_func(mpc_number):
    # f = io.StringIO()
    # with redirect_stdout(f):
    #     fit = sbpy_phase_fit.phase_fit(mpc_number,push_fit_flag=True,hide_warning_flag=True)
    #     fit.calculate()
    # # del f

    # sys.stdout = open(str(os.getpid()) + ".out", "w")
    with open("tmp", "w") as f:
        sys.stdout = f
        fit = sbpy_phase_fit.phase_fit(mpc_number,push_fit_flag=True,hide_warning_flag=True)
        check=fit.calculate()

    return check

# connect to the database
from calculate_phase import atlas_database_connection

cnx=atlas_database_connection.database_connection().connect()
# perform the query and store results as a dataframe
qry="select mpc_number from atlas_objects where mpc_number is not null;"
df=pd.read_sql_query(qry,cnx)
# df=df.dropna()
print(df)
cnx.close()

#mpc_number_list=np.random.choice(df['mpc_number'].astype(int),10,replace=False)
mpc_number_list=df['mpc_number'].astype(int)
print(mpc_number_list)
print(len(mpc_number_list))

#mpc_number_list=mpc_number_list[:100]

mpc_number_list=list(mpc_number_list)
print(mpc_number_list)
print(len(mpc_number_list))

#exit()

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

# create a file that records the date, the start and end mpc numbers in a batch, the length of time to complete a batch, the number of jobs in a batch and the number of objects done per sec
today=datetime.datetime.now()
fname="fit_phase_curves_bulk_record_{}-{}-{}.csv".format(today.day,today.month,today.year)
print(fname)
of=open(fname,"w")
of.write("date,mpc1,mpc2,t(s),N_jobs,rate(1/s)\n")

jobs_done=0
while len(mpc_number_list)>0:

    t_start=time.time()
    sub_list=mpc_number_list[:threads]
    mpc_number_list=mpc_number_list[threads:]

    print("run parallel, {} threads".format(threads))
    print("fit mpc numbers {} - {}".format(sub_list[0],sub_list[-1]))
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
