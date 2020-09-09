import pandas as pd
import mysql.connector
import numpy as np
import time
import multiprocessing
from multiprocessing import Pool
import subprocess

def phase_fit_func(mpc_number):
    py_cmd="python sbpy_phase_fit_data_clip_push.py -n {}".format(mpc_number)
    print(py_cmd)
    subprocess.Popen(py_cmd,shell=True).wait()
    return

"""
# connect to the database
config = {
  'user': 'af',
  'password': 'afPass',
  'host': 'localhost',
  'port': '3308',
  'database': 'atlas_moving_objects',
  'raise_on_warnings': True
}
cnx = mysql.connector.connect(**config)

# perform the query and store results as a dataframe
qry="select mpc_number from atlas_objects;"
df=pd.read_sql_query(qry,cnx)
df=df.dropna()
print(df)

cnx.close()

mpc_number_list=np.random.choice(df['mpc_number'].astype(int),10,replace=False)
print(mpc_number_list)
"""

mpc_number_list=[339933, 370192, 416678, 238537, 181325, 358275,  68172, 383854,  41777, 215944]
# mpc_number_list=[339933, 370192, 416678, 238537]
print(mpc_number_list)

# add mpc_number or date requirements to run?
# run in reverse date order so that objects that have not been calculated are done first?

# print("run in serial")
# t_start=time.time()
# for n in mpc_number_list:
#     phase_fit_func(n)
# t_end=time.time()
# print ("N_jobs = {}".format(len(mpc_number_list)))
# print ("t_pool = {}".format(t_end-t_start))

# set up pool and fit phase
# threads=4
threads=multiprocessing.cpu_count()

jobs_done=0
t_start=time.time()
while len(mpc_number_list)>0:
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
