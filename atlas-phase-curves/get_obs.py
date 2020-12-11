""" Calls the atlas_SQL_query function to get an object's observations from ATLAS """

import pandas as pd
from calculate_phase import atlas_SQL_query_df
from calculate_phase import atlas_database_connection
from optparse import OptionParser

parser = OptionParser()
parser.add_option( "-n", "--mpc-number", dest="mpc_number", help="mpc_number", metavar="MPC_NUMBER" ) # mpc number of object to fit
(options,args)=parser.parse_args()

if options.mpc_number:
    mpc_number=int(options.mpc_number)
else:
    mpc_number="4986"

# connect to database
cnx=atlas_database_connection.database_connection().connect()

# run the SQL query
df_data=atlas_SQL_query_df.atlas_SQL_query(mpc_number=mpc_number,cnx=cnx)

# save dataframe
print(df_data)
print(list(df_data))
df_data.to_csv("results_analysis/obs/df_data_{}.csv".format(mpc_number))
# print(pd.read_csv("df_data_{}.csv".format(mpc),index_col=0))
