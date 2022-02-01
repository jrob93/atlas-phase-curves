""" Calls the atlas_SQL_query function to get an object's observations from ATLAS """

import pandas as pd
from calculate_phase import atlas_SQL_query_df
from calculate_phase import atlas_database_connection
from optparse import OptionParser
import time
import numpy as np

parser = OptionParser()
parser.add_option( "-n", "--mpc-number", dest="mpc_number", help="mpc_number", metavar="MPC_NUMBER" ) # mpc number of object to fit
parser.add_option( "-N", "--name", dest="name", help="name", metavar="NAME" ) # NAME of object to fit
(options,args)=parser.parse_args()

if options.mpc_number:
    mpc_number=int(options.mpc_number)
else:
    # mpc_number="4986"
    mpc_number=False
if options.name:
    name=options.name
else:
    name=False

# connect to database
cnx=atlas_database_connection.database_connection().connect()

# primid,orbid=atlas_SQL_query_df.get_unique_ids(cnx,mpc_number,name)
# print(primid,orbid)
ids=atlas_SQL_query_df.get_unique_ids(cnx,mpc_number,name)
print(ids)
print(ids["primaryId"])
print(ids["orbital_elements_id"])
# exit()

# # run the SQL query using mpc number
# start = time.process_time()
# df_data=atlas_SQL_query_df.atlas_SQL_query(mpc_number=mpc_number,cnx=cnx)
# time1=time.process_time() - start

# # time single function for orbid and data
# start = time.process_time()
# df_data=atlas_SQL_query_df.atlas_SQL_query_test(mpc_number=mpc_number,cnx=cnx,name=name)
# time2=time.process_time() - start

# time separate funcs
start = time.process_time()
orbid=atlas_SQL_query_df.get_orb_elements_id(cnx,mpc_number,name)
# df_data=atlas_SQL_query_df.atlas_SQL_query_orbid(cnx,orbid)
df_data=atlas_SQL_query_df.atlas_SQL_query_orbid_expname(cnx,orbid)
time3=time.process_time() - start

print(df_data[np.absolute(df_data["galactic_latitude"])<10][["expname",'ra_deg','dec_deg']])
# save dataframe
print(df_data)
print(list(df_data))
print(df_data[["expname",'ra_deg','dec_deg']])

print(df_data[["m",'merr','galactic_latitude']])

# df_data.to_csv("results_analysis/obs/df_data_{}.csv".format(mpc_number))

# print(time1/time2,time2/time2, time3/time2)

# err_cut=0.01
err_cut=0.005
N_merr = len(df_data[df_data["merr"]<err_cut])
print("N merr<{} : {}".format(err_cut,N_merr))


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
fig = plt.figure()
gs = gridspec.GridSpec(1, 1)
ax1 = plt.subplot(gs[0,0])

# ax1.hist(df_data["merr"],bins="auto")
# ax1.axvline(err_cut,c="r")
# ax1.set_xlabel("N")
# ax1.set_ylabel("merr")

# ax1.scatter(df_data["m"],df_data["merr"],c=df_data["mjd"])
# ax1.axhline(0.01,c="r",label="merr=0.01")
# ax1.axvline(np.median(df_data["m"]),label="median m")
# ax1.legend()
# ax1.set_xlabel("m")
# ax1.set_ylabel("merr")

ax1.errorbar(df_data['phase_angle'],df_data['reduced_mag'],df_data['merr'], fmt='ko',label="data",zorder=0,markersize="2")
ax1.legend()
ax1.set_xlabel("phase_angle")
ax1.set_ylabel("reduced_mag")

plt.show()


print()
