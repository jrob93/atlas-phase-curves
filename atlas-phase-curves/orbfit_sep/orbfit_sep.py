"""
Investigate the properties of the orbfit separation, large separations tend to occur for poor photometry
"""

import pandas as pd
from calculate_phase import atlas_SQL_query_df
from calculate_phase import atlas_database_connection
import numpy as np

# select all objects
qry="select name from atlas_objects;"
df_names=pd.read_sql_query(qry,cnx)
print(df_names)
cnx.close()


for i in range(len(df_names)):

    name = df_names.iloc[i]["name"]
    mpc_number=False

    # connect to database
    cnx=atlas_database_connection.database_connection().connect()

    # get unique object id
    orbid=atlas_SQL_query_df.get_orb_elements_id(cnx,mpc_number,name)

    # query to get just the orbfit separation
    qry = """SELECT d.orbfit_separation_arcsec
            FROM
                dophot_photometry d,
                orbfit_positions o,
                atlas_exposures a
            WHERE
                AND d.orbfit_postions_id = o.primaryId
                AND a.expname = d.expname
                AND o.orbital_elements_id = {};
    """.format(orbid)
