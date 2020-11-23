import pandas as pd
import mysql.connector
import numpy as np
import time
import multiprocessing
from multiprocessing import Pool
import subprocess

def phase_fit_func(mpc_number):
    # !!! use try/except to append to list of errors?
    py_cmd="python sbpy_phase_fit_data_clip_push.py -n {} > tmp".format(mpc_number)
    print(py_cmd)
    subprocess.Popen(py_cmd,shell=True).wait()
    return

# """
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

#mpc_number_list=np.random.choice(df['mpc_number'].astype(int),10,replace=False)
mpc_number_list=df['mpc_number'].astype(int)
#mpc_number_list=mpc_number_list[mpc_number_list>11350]
mpc_number_list=mpc_number_list[mpc_number_list>95226]
# print(np.array(mpc_number_list[:50]))
# exit()
# """

# mpc_number_list=[339933, 370192, 416678, 238537, 181325, 358275,  68172, 383854,  41777, 215944]
# mpc_number_list=[95227,95228,95229,95230,95231,95232,95233,95234,95235,95236,95237,95238
# ,95239,95240,95241,95242,95243,95244,95245,95246,95247,95248,95249,95250
# ,95251,95252,95253,95254,95255,95256,95257,95258,95259,95260,95261,95262
# ,95263,95264,95265,95266,95267,95268,95269,95270,95271,95272,95273,95274
# ,95275,95276]
#print(mpc_number_list)

mpc_number_list=mpc_number_list[:120]
N_batch=60

# from multiprocessing import Lock, Process, Queue, current_process
# import time
# import queue # imported for using queue.Empty exception
#
#
# def do_job(tasks_to_accomplish, tasks_that_are_done):
#     while True:
#         try:
#             '''
#                 try to get task from the queue. get_nowait() function will
#                 raise queue.Empty exception if the queue is empty.
#                 queue(False) function would do the same task also.
#             '''
#             task = tasks_to_accomplish.get_nowait()
#             # print(task)
#             phase_fit_func(task)
#
#         except queue.Empty:
#
#             break
#         else:
#             '''
#                 if no exception has been raised, add the task completion
#                 message to task_that_are_done queue
#             '''
#             # print(task)
#             # tasks_that_are_done.put(task + ' is done by ' + current_process().name)
#             tasks_that_are_done.put("{} is done by {}".format(task,current_process().name))
#             # time.sleep(.5)
#     return True
#
#
# def main():
#
#     t_start=time.time()
#
#     number_of_task = 10
#     number_of_processes = 4
#     tasks_to_accomplish = Queue()
#     tasks_that_are_done = Queue()
#     processes = []
#
#     # for i in range(number_of_task):
#     #     tasks_to_accomplish.put("Task no " + str(i))
#     for i in mpc_number_list:
#         tasks_to_accomplish.put(i)
#
#     # creating processes
#     for w in range(number_of_processes):
#         p = Process(target=do_job, args=(tasks_to_accomplish, tasks_that_are_done))
#         processes.append(p)
#         p.start()
#
#     # completing process
#     for p in processes:
#         p.join()
#
#     # print the output
#     jobs_done=0
#     while not tasks_that_are_done.empty():
#         print(tasks_that_are_done.get())
#         jobs_done+=1
#
#     t_end=time.time()
#     t_elapsed=t_end-t_start
#
#     print("N jobs done = {}".format(jobs_done))
#     print("time elapsed = {}".format(t_elapsed))
#     print("{} objects/second".format(float(jobs_done)/t_elapsed))
#
#     return True
#
#
# if __name__ == '__main__':
#     main()

""" based on https://www.journaldev.com/15631/python-multiprocessing-example"""

from multiprocessing import Lock, Process, Queue, current_process
import time
import queue # imported for using queue.Empty exception


def do_job(tasks_to_accomplish):#, tasks_that_are_done):
    while True:
        try:
            '''
                try to get task from the queue. get_nowait() function will
                raise queue.Empty exception if the queue is empty.
                queue(False) function would do the same task also.
            '''
            task = tasks_to_accomplish.get_nowait()
            phase_fit_func(task)

        except queue.Empty:
            break
    return True

jobs_done=0
t_start=time.time()
while len(mpc_number_list)>0:
    sub_list=mpc_number_list[:N_batch]
    mpc_number_list=mpc_number_list[N_batch:]

    number_of_processes = multiprocessing.cpu_count()
    tasks_to_accomplish = Queue()
    processes = []

    # for i in mpc_number_list:
    for i in sub_list:
        tasks_to_accomplish.put(i)

    # creating processes
    for w in range(number_of_processes):
        p = Process(target=do_job, args=(tasks_to_accomplish,))#, tasks_that_are_done))
        processes.append(p)
        p.start()

    # completing process
    for p in processes:
        p.join()

    t_end=time.time()
    t_elapsed=t_end-t_start
    # jobs_done=len(mpc_number_list)
    jobs_done+=len(sub_list)

    print("N jobs done = {}".format(jobs_done))
    print("time elapsed = {}".format(t_elapsed))
    print("{} objects/second".format(float(jobs_done)/t_elapsed))
