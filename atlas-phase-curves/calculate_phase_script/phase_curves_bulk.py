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

def phase_fit_func(mpc_number):
    # !!! use try/except to append to list of errors?
    #py_cmd="python sbpy_phase_fit_data_clip_push.py -n {} > tmp".format(mpc_number)
    # py_cmd="python sbpy_phase_fit.py -n {} -w 1 > tmp".format(mpc_number)
    # py_cmd="python sbpy_phase_fit.py -n {}".format(mpc_number)
    py_cmd="python -m cProfile -s 'tottime' sbpy_phase_fit.py -n {}".format(mpc_number)
    print(py_cmd)
    subprocess.Popen(py_cmd,shell=True).wait()
    return

# # """
# # connect to the database
# config = {
#   'user': 'af',
#   'password': 'afPass',
#   'host': 'localhost',
#   'port': '3308',
#   'database': 'atlas_moving_objects',
#   'raise_on_warnings': True
# }
# cnx = mysql.connector.connect(**config)
#
# # perform the query and store results as a dataframe
# # qry="select mpc_number from atlas_objects;"
# qry="select mpc_number from atlas_objects where mpc_number is not null;"
#
# # of=open("../create_table/atlas_objects_fields.txt","r")
# # fields=[f.rstrip() for f in of.readlines()]
# # of.close()
# # # fields = [f for f in fields if f.startswith("phase_curve_G_B89") or f.startswith("phase_curve_H_B89")]
# # # fields = [f for f in fields if f.startswith("phase_curve_G_B89") or f.startswith("phase_curve_H_B89")]# and "err" not in f]
# # fields = [f for f in fields if f.startswith("phase_curve_G") or f.startswith("phase_curve_H")]# and "err" not in f]
# # # fields_null = " is null or ".join(fields)
# # fields_null = " is null and ".join(fields)
# # fields_check = ",".join(fields)
# # print(fields_null)
# # qry="""select mpc_number,{} from atlas_phase_fits where {} is null;""".format(fields_check,fields_null)
# # print(qry)
#
# df=pd.read_sql_query(qry,cnx)
# # df=df.dropna()
# print(df)
#
# cnx.close()
#
# #mpc_number_list=np.random.choice(df['mpc_number'].astype(int),10,replace=False)
# mpc_number_list=df['mpc_number'].astype(int)
# print(mpc_number_list)
# print(len(mpc_number_list))
# # exit()
#
# #mpc_number_list=mpc_number_list[mpc_number_list>11350]
# #mpc_number_list=mpc_number_list[mpc_number_list>95226]
# #mpc_number_list=mpc_number_list[mpc_number_list>197120]
# # mpc_number_list=mpc_number_list[mpc_number_list>118331]
#
# #exit()
# # """

# mpc_number_list=[339933, 370192, 416678, 238537, 181325, 358275,  68172, 383854,  41777, 215944]
mpc_number_list=[339933]#, 370192, 416678, 238537]
# mpc_number_list=[95227,95228,95229,95230,95231,95232,95233,95234,95235,95236,95237,95238
# ,95239,95240,95241,95242,95243,95244,95245,95246,95247,95248,95249,95250
# ,95251,95252,95253,95254,95255,95256,95257,95258,95259,95260,95261,95262
# ,95263,95264,95265,95266,95267,95268,95269,95270,95271,95272,95273,95274
# ,95275,95276]
#mpc_number_list=mpc_number_list[:1000]
# mpc_number_list=mpc_number_list[:40]

mpc_number_list=list(mpc_number_list)
print(mpc_number_list)
print(len(mpc_number_list))

# exit()

# add mpc_number or date requirements to run?
# run in reverse date order so that objects that have not been calculated are done first?

print("run in serial")
t_start=time.time()
count=0
for n in mpc_number_list:
    phase_fit_func(n)
    # if count>3:
    #     break
    # else:
    #     count+=1
    count+=1
t_end=time.time()
print ("N_jobs = {}".format(count))
print ("t_pool = {}".format(t_end-t_start))
print("{} objects/second".format(float(count)/(t_end-t_start)))
exit()

# set up pool and fit phase
#threads=40
threads=multiprocessing.cpu_count()

# create a file that records the date, the start and end mpc numbers in a batch, the length of time to complete a batch, the number of jobs in a batch and the number of objects done per sec
today=datetime.datetime.now()
fname="phase_curve_bulk_record_{}-{}-{}.csv".format(today.day,today.month,today.year)
print(fname)
of=open(fname,"w")
of.write("date,mpc1,mpc2,t(s),N_jobs,rate(1/s)\n")
# exit()

jobs_done=0
# t_start=time.time()
while len(mpc_number_list)>0:

    t_start=time.time()
    sub_list=mpc_number_list[:threads]
    mpc_number_list=mpc_number_list[threads:]

    print("run parallel, {} threads".format(threads))
    # t_start=time.time()
    pool = Pool(threads)
    multiple_results = [pool.apply_async(phase_fit_func, args=(n,)) for n in sub_list]
    pool.close()
    pool.join()
    t_end=time.time()
    print()
    # print("N_jobs = {}".format(len(sub_list)))
    # print("t_pool = {}".format(t_end-t_start))
    jobs_done+=len(sub_list)
    t_elapsed=t_end-t_start
    print("N jobs done = {}".format(jobs_done))
    print("time elapsed = {}".format(t_elapsed))
    print("{} objects/second".format(float(jobs_done)/t_elapsed))

    print()

    of.write("{},{},{},{},{},{}\n".format(datetime.datetime.now(),sub_list[0],sub_list[-1],t_elapsed,len(sub_list),float(len(sub_list))/t_elapsed))
    of.flush()

of.close()
    # exit()
