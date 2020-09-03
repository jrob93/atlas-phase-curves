import pandas as pd
import mysql.connector
import subprocess
import numpy as np

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

# mpc_number_list=[339933, 370192, 416678, 238537, 181325, 358275,  68172, 383854,  41777, 215944]
mpc_number_list=[339933, 370192, 416678, 238537]#, 181325, 358275,  68172, 383854,  41777, 215944]

# cmds=[]
for N_mpc in mpc_number_list:
    cmd="python sbpy_phase_fit_data_clip_push.py -n {}".format(N_mpc)
    print(cmd)
    # cmds.append(cmd)
    subprocess.Popen(cmd,shell=True).wait()
    # exit()

#     # exit()
# # exit()
# print(cmds)
#
# import multiprocessing
# import subprocess
#
# def work(cmd):
#     return subprocess.Popen(cmd, shell=True)
#
# if __name__ == '__main__':
#     count = multiprocessing.cpu_count()
#     pool = multiprocessing.Pool(processes=count)
#     print(count, pool)
#     pool.map(work, cmds)

# # set up pool and find orbits
# # works best if the command is functionalised!
# t_start=time.time()
# pool = Pool(threads)
# multiple_results = [pool.apply_async(bin_orb_func, args=("/".join(d.split("/")[:-1]),save_path,d.split("/")[-1])) for d in all_run_dirs]
# pool.close()
# pool.join()
# t_end=time.time()
# print "N_jobs = {}".format(len(all_run_dirs))
# print "t_pool = {}".format(t_end-t_start)
