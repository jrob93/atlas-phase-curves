"""
Use this script to query the atlas db, assuming you have set up a connection to dormammu via port 3308
(follow Dave's connection instructions)
The query is parsed by pandas to create a handy dataframe

Install mysql package with:
pip install mysql-connector-python

Add option for object number or name

This function version returns the dataframe
"""

# compare Alan and Dave's queries, e.g. Dave uses orbital_elements_id to select data, is this faster?
# write function to accept name or mpc number

from pandas import read_sql_query

def atlas_SQL_query(cnx,mpc_number=4986,filter="all"):

    # print(cnx)
    # orbital_elements_id=mpc_number

    # write out the sql query

    # REVIEW THIS QUERY, E.G. WHY d.expname?
    if filter=="all":
        sqlQuery_alan = u"""
        SELECT d.mjd, d.m,  d.dfitmag as merr, a.filter, o.observer_distance, o.heliocentric_distance, o.phase_angle, d.m - 5*log10(o.heliocentric_distance*o.observer_distance) as reduced_mag, o.apparent_mag, o.galactic_latitude
        FROM
            dophot_photometry d,
            orbfit_positions o,
            atlas_exposures a
        WHERE
            d.orbfit_postions_id = o.primaryId
                AND a.expname = d.expname
                AND o.mpc_number = %(mpc_number)s;
        """ % locals()
    else:
        sqlQuery_alan = u"""
        SELECT d.mjd, d.m,  d.dfitmag as merr, a.filter, o.observer_distance, o.heliocentric_distance, o.phase_angle, d.m - 5*log10(o.heliocentric_distance*o.observer_distance) as reduced_mag, o.apparent_mag, o.galactic_latitude
        FROM
            dophot_photometry d,
            orbfit_positions o,
            atlas_exposures a
        WHERE
            filter = '%(filter)s'
                AND d.orbfit_postions_id = o.primaryId
                AND a.expname = d.expname
                AND o.mpc_number = %(mpc_number)s;
        """ % locals()

    # perform the query and store results as a dataframe
    sqlQuery=sqlQuery_alan
    # df=pd.read_sql_query(sqlQuery,cnx)
    df=read_sql_query(sqlQuery,cnx)

    # cnx.close()

    return df

# if __name__ == "__main__":
#
#     import pandas as pd
#     from atlas_database_connection import database_connection
#     cnx=database_connection().connect()
#     print(cnx)
#     # mpc=4985
#     # mpc=121514
#     mpc=136199
#     df_data=atlas_SQL_query(mpc_number=mpc,cnx=cnx)
#     print(df_data)
#     print(list(df_data))
#     df_data.to_csv("df_data_{}.csv".format(mpc))
#     print(pd.read_csv("df_data_{}.csv".format(mpc),index_col=0))
#     # print("atlas_SQL_query_df!!!!")
#     exit()
