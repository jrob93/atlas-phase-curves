
import numpy as np
from sbpy.photometry import HG, HG1G2, HG12, HG12_Pen16
import pandas as pd
import os
import mysql.connector
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime
import os.path

# table="test_table2"
# fname="{}_12_10_20.csv".format(table)
#
# # connect to the database to write to
# config = {
#   'user': 'af',
#   'password': 'afPass',
#   'host': 'localhost',
#   'port': '3308',
#   'database': 'jamie_test_db',
#   'raise_on_warnings': True
# }

# table="atlas_phase_fits"
table="atlas_phase_fits_orbs"
# fname="{}_16_2_2021.csv".format(table)
utc_date_now=datetime.datetime.utcnow()#.strftime("%Y-%m-%d %H:%M:%S")
fname="{}_{}_{}_{}.csv".format(table,utc_date_now.day,utc_date_now.month,utc_date_now.year)

# check to prevent overwriting of file
if os.path.isfile(fname):
    print("{} already exists, stop".format(fname))

    # df=pd.read_csv(fname)

    print(df)
    df=pd.read_csv(fname,index_col=0)
    print(list(df))
    print(df[["mpc_number","mpc_number.1","name","name.1","primaryId","primaryId.1"]])

    exit()

# CHANGE FNAME TO BE DATE THAT FIT WAS RUN, NOT DATE IT WAS DOWNLOADED
# get date from fit_phase_curves_bulk_record file
print(fname)
# exit()

# connect to the database to write to
config = {
  'user': 'af',
  'password': 'afPass',
  'host': 'localhost',
  'port': '3308',
  'database': 'atlas_moving_objects',
  'raise_on_warnings': True
}

cnx = mysql.connector.connect(**config)
cursor =cnx.cursor()

# qry="select * from {};".format(table)
# qry="""select a.*,
#  o.mpc_number,o.name,o.a_semimajor_axis,o.e_eccentricity,o.i_inclination_deg
#  from atlas_phase_fits a, orbital_elements o where a.orbital_elements_id=o.primaryId;""" # limit 10;"""
qry="""select a.*,
 o.a_semimajor_axis,o.e_eccentricity,o.i_inclination_deg
 from atlas_phase_fits a, orbital_elements o where a.orbital_elements_id=o.primaryId;""" # limit 10;"""

df=pd.read_sql_query(qry,cnx)
print(len(df))
# df=df.dropna() # If any NA values are present, drop that row or column.
print(len(df))
df.to_csv(fname)
exit()


# # drop high G
# df = df[df['phase_curve_G_B89_o']<2.0]

print(df.to_string())
# print(len(np.isnan(df['detection_count'])))
# print(df[np.isnan(df['detection_count'])][['detection_count','detection_count_o','detection_count_c','phase_curve_G_B89_o','phase_curve_G_B89_c']])

# # Plot the distribution
# fig = plt.figure()
# gs = gridspec.GridSpec(2,2)
# ax1 = plt.subplot(gs[0,0])
# ax2 = plt.subplot(gs[1,0])
# ax3 = plt.subplot(gs[0,1])
# ax4 = plt.subplot(gs[1,1])
#
# ax1.hist(df['phase_curve_H_B89_o'],bins="auto")
# ax2.hist(df['phase_curve_G_B89_o'],bins="auto")
# ax3.hist(df['phase_curve_H_err_B89_o'],bins="auto")
# ax4.hist(df['phase_curve_G_err_B89_o'],bins="auto")
#
# fig.suptitle('HG B89, 2-sigma clip, N_tot={}'.format(len(df)))
#
# ax1.set_xlabel('phase_curve_H_B89_o')
# ax2.set_xlabel('phase_curve_G_B89_o')
# ax3.set_xlabel('phase_curve_H_err_B89_o')
# ax4.set_xlabel('phase_curve_G_err_B89_o')
#
# plt.show()
