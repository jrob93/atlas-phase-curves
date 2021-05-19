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

def update_atlas_objects(cnx):

    """ Call the update_atlas_objects procedure written by Dave. Does what it says on the tin - updates fields such as detection_count
    Takes a long time to run!
    """

    qry="CALL update_atlas_objects"
    # qry="SELECT * from atlas_objects limit 100"
    print(qry)
    cursor=cnx.cursor()
    cursor.execute(qry)
    cnx.commit()
    # print(cursor)
    # for c in cursor:
    #     print(c)

    return

def get_orb_elements_id(cnx,mpc_number=False,name=False):
    """ Function to retrieve the orbital_elements_id from either mpc_number or name.
     name case should not matter for SQL """

    print("mpc_number={}".format(mpc_number))
    print("name={}".format(name))

    # use mpc_number or name to get orbital_elements_id
    if name:
        qry="SELECT orbital_elements_id FROM atlas_objects WHERE name=\"{}\";".format(name)
    else:
        qry="SELECT orbital_elements_id FROM atlas_objects WHERE mpc_number=\"{}\";".format(mpc_number)

    print(qry)
    cursor=cnx.cursor()
    cursor.execute(qry)
    row = cursor.fetchone()
    orbital_elements_id=int(row[0])

    print("orbital_elements_id={}".format(orbital_elements_id))

    return orbital_elements_id

def get_unique_ids(cnx,mpc_number=False,name=False):
    """ Function to retrieve the orbital_elements_id from either mpc_number or name.
     name case should not matter for SQL """

    print("mpc_number={}".format(mpc_number))
    print("name={}".format(name))

    # use mpc_number or name to get orbital_elements_id
    if name:
        qry="SELECT primaryId,orbital_elements_id FROM atlas_objects WHERE name=\"{}\";".format(name)
    else:
        qry="SELECT primaryId,orbital_elements_id FROM atlas_objects WHERE mpc_number=\"{}\";".format(mpc_number)

    print(qry)
    cursor=cnx.cursor()
    cursor.execute(qry)
    row = cursor.fetchone()
    print(row)
    primaryId=int(row[0])
    orbital_elements_id=int(row[1])

    print("primaryId={}".format(primaryId))
    print("orbital_elements_id={}".format(orbital_elements_id))

    return {"primaryId":primaryId,"orbital_elements_id":orbital_elements_id}

def atlas_SQL_query_orbid(cnx,orbital_elements_id,filter="all"):
    """ Query that accept orb orbital_elements_id """

    # # use mpc_number or name to get orbital_elements_id
    # if name:
    #     qry="SELECT orbital_elements_id FROM atlas_objects WHERE name=\"{}\";".format(name)
    # else:
    #     qry="SELECT orbital_elements_id FROM atlas_objects WHERE mpc_number=\"{}\";".format(mpc_number)
    #
    # cursor=cnx.cursor()
    # cursor.execute(qry)
    # row = cursor.fetchone()
    # orbital_elements_id=row[0]
    #
    # print(cnx)
    # orbital_elements_id=mpc_number

    # write out the sql query

    # REVIEW THIS QUERY, E.G. WHY d.expname?
    # remove unneccesary data, e.g. apparent_mag

    if filter=="all":

        # sqlQuery_alan = u"""
        # SELECT d.mjd, d.m,  d.dfitmag as merr, a.filter, o.observer_distance, o.heliocentric_distance, o.phase_angle, d.m - 5*log10(o.heliocentric_distance*o.observer_distance) as reduced_mag, o.apparent_mag, o.galactic_latitude
        # FROM
        #     dophot_photometry d,
        #     orbfit_positions o,
        #     atlas_exposures a
        # WHERE
        #     d.orbfit_postions_id = o.primaryId
        #         AND a.expname = d.expname
        #         AND o.mpc_number = %(mpc_number)s;
        # """ % locals()

        # based on dave's query that uses orbital_elements_id
        sqlQuery_dave = u"""
        SELECT d.mjd, d.m,  d.dfitmag as merr, a.filter, o.observer_distance, o.heliocentric_distance, o.phase_angle, d.m - 5*log10(o.heliocentric_distance*o.observer_distance) as reduced_mag, o.apparent_mag, o.galactic_latitude
        FROM
            dophot_photometry d,
            orbfit_positions o,
            atlas_exposures a
        WHERE
            d.orbfit_postions_id = o.primaryId
                AND a.expname = d.expname
                AND o.orbital_elements_id = %(orbital_elements_id)s;
        """ % locals()

    else:
        # sqlQuery_alan = u"""
        # SELECT d.mjd, d.m,  d.dfitmag as merr, a.filter, o.observer_distance, o.heliocentric_distance, o.phase_angle, d.m - 5*log10(o.heliocentric_distance*o.observer_distance) as reduced_mag, o.apparent_mag, o.galactic_latitude
        # FROM
        #     dophot_photometry d,
        #     orbfit_positions o,
        #     atlas_exposures a
        # WHERE
        #     filter = '%(filter)s'
        #         AND d.orbfit_postions_id = o.primaryId
        #         AND a.expname = d.expname
        #         AND o.mpc_number = %(mpc_number)s;
        # """ % locals()

        sqlQuery_dave = u"""
        SELECT d.mjd, d.m,  d.dfitmag as merr, a.filter, o.observer_distance, o.heliocentric_distance, o.phase_angle, d.m - 5*log10(o.heliocentric_distance*o.observer_distance) as reduced_mag, o.apparent_mag, o.galactic_latitude
        FROM
            dophot_photometry d,
            orbfit_positions o,
            atlas_exposures a
        WHERE
            filter = '%(filter)s'
                AND d.orbfit_postions_id = o.primaryId
                AND a.expname = d.expname
                AND o.orbital_elements_id = %(orbital_elements_id)s;
        """ % locals()

    # perform the query and store results as a dataframe
    # sqlQuery=sqlQuery_alan
    sqlQuery=sqlQuery_dave
    print(sqlQuery)
    # df=pd.read_sql_query(sqlQuery,cnx)
    df=read_sql_query(sqlQuery,cnx)

    # cnx.close()

    return df # should also return orbital_elements_id!!!

def atlas_SQL_query_orbid_expname(cnx,orbital_elements_id,filter="all"):
    """ Query that accept orb orbital_elements_id
    Grabs the exposure name of each detection, useful for creating postage stamps"""

    # based on dave's query that uses orbital_elements_id
    sqlQuery_dave = u"""
    SELECT d.expname,o.dec_deg,o.ra_deg,d.mjd, d.m,  d.dfitmag as merr, a.filter, o.observer_distance, o.heliocentric_distance, o.phase_angle, d.m - 5*log10(o.heliocentric_distance*o.observer_distance) as reduced_mag, o.apparent_mag, o.galactic_latitude
    FROM
        dophot_photometry d,
        orbfit_positions o,
        atlas_exposures a
    WHERE
        d.orbfit_postions_id = o.primaryId
            AND a.expname = d.expname
            AND o.orbital_elements_id = %(orbital_elements_id)s;
    """ % locals()

    # perform the query and store results as a dataframe
    sqlQuery=sqlQuery_dave
    print(sqlQuery)
    df=read_sql_query(sqlQuery,cnx)

    return df # should also return orbital_elements_id!!!

def atlas_SQL_query_test(cnx,mpc_number=False,filter="all",name=False):
    """ Test grabbing orb id and the data in the same function """

    # use mpc_number or name to get orbital_elements_id
    if name:
        qry="SELECT orbital_elements_id FROM atlas_objects WHERE name=\"{}\";".format(name)
    else:
        qry="SELECT orbital_elements_id FROM atlas_objects WHERE mpc_number=\"{}\";".format(mpc_number)

    cursor=cnx.cursor()
    cursor.execute(qry)
    row = cursor.fetchone()
    orbital_elements_id=row[0]

    print(cnx)

    if filter=="all":

        # based on dave's query that uses orbital_elements_id
        sqlQuery_dave = u"""
        SELECT d.mjd, d.m,  d.dfitmag as merr, a.filter, o.observer_distance, o.heliocentric_distance, o.phase_angle, d.m - 5*log10(o.heliocentric_distance*o.observer_distance) as reduced_mag, o.apparent_mag, o.galactic_latitude
        FROM
            dophot_photometry d,
            orbfit_positions o,
            atlas_exposures a
        WHERE
            d.orbfit_postions_id = o.primaryId
                AND a.expname = d.expname
                AND o.orbital_elements_id = %(orbital_elements_id)s;
        """ % locals()

    else:

        sqlQuery_dave = u"""
        SELECT d.mjd, d.m,  d.dfitmag as merr, a.filter, o.observer_distance, o.heliocentric_distance, o.phase_angle, d.m - 5*log10(o.heliocentric_distance*o.observer_distance) as reduced_mag, o.apparent_mag, o.galactic_latitude
        FROM
            dophot_photometry d,
            orbfit_positions o,
            atlas_exposures a
        WHERE
            filter = '%(filter)s'
                AND d.orbfit_postions_id = o.primaryId
                AND a.expname = d.expname
                AND o.orbital_elements_id = %(orbital_elements_id)s;
        """ % locals()

    # perform the query and store results as a dataframe
    # sqlQuery=sqlQuery_alan
    sqlQuery=sqlQuery_dave
    print(sqlQuery)
    # df=pd.read_sql_query(sqlQuery,cnx)
    df=read_sql_query(sqlQuery,cnx)

    # cnx.close()

    return df # should also return orbital_elements_id!!!

def atlas_SQL_query(cnx,mpc_number=4986,filter="all"):
    """ Original query that searches by mpc number """

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
    print(sqlQuery)
    # df=pd.read_sql_query(sqlQuery,cnx)
    df=read_sql_query(sqlQuery,cnx)

    # cnx.close()

    return df # should also return orbital_elements_id!!!

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
