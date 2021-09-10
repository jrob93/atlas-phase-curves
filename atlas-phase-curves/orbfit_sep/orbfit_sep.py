"""
Investigate the properties of the orbfit separation, large separations tend to occur for poor photometry
"""

import pandas as pd
import numpy as np
import os
import sys
sys.path.append("/Users/jrobinson/atlas-phase-curves/atlas-phase-curves")
from calculate_phase import atlas_SQL_query_df
from calculate_phase import atlas_database_connection

orbfit_file = "df_orbfit_sep.csv"

# connect to database
cnx=atlas_database_connection.database_connection().connect()

# select all objects
name_file = "df_names.csv"
if os.path.isfile(name_file):
    df_names = pd.read_csv(name_file,index_col=0)
else:
    qry="select name from atlas_objects;"
    df_names=pd.read_sql_query(qry,cnx)
    cnx.close()
    df_names.to_csv(name_file)

print(df_names)

# med_sep_list = []
# name_list = []

f = open(orbfit_file,"w")
f.write(",name,med_orbfit_sep(arcsec)\n")

for i in range(len(df_names)):

    name = df_names.iloc[i]["name"]
    mpc_number=False

    print(i,name)

    # get unique object id
    orbid=atlas_SQL_query_df.get_orb_elements_id(cnx,mpc_number,name)

    # query to get just the orbfit separation
    qry = """SELECT d.orbfit_separation_arcsec
            FROM
                dophot_photometry d,
                orbfit_positions o,
                atlas_exposures a
            WHERE
                d.orbfit_postions_id = o.primaryId
                AND a.expname = d.expname
                AND o.orbital_elements_id = {};
    """.format(orbid)
    # print(qry)

    df=pd.read_sql_query(qry,cnx)
    # print(df)

    med_sep = np.median(df["orbfit_separation_arcsec"])
    # print(med_sep)

    # med_sep_list.append(med_sep)
    # name_list.append(name)

    f.write("{},{},{}\n".format(i,name,med_sep))
    f.flush()

    # if i>5:
    #     break

f.close()

# df_seps = pd.DataFrame({'name': np.array(name_list), 'orbfit_separation_arcsec': np.array(med_sep_list)})
# df_seps.to_csv("df_orbfit_sep.csv")
df_seps = pd.read_csv(orbfit_file,index_col=0)
print(df_seps)
