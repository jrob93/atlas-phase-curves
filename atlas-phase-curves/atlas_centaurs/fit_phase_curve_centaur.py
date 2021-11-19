import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import importlib
import os
from astroquery.jplhorizons import Horizons
from sbpy.data import Names
from astropy.time import Time

sys.path.append("/Users/jrobinson/atlas-phase-curves/atlas-phase-curves")
from calculate_phase import sbpy_phase_fit
importlib.reload(sbpy_phase_fit)

name_list = ["Hidalgo",
"2016 ND21",
"Chiron",
"Chariklo",
"Narcissus",
"Bienor",
"Echeclus",
"2012 TW236",
"2009 YF7",
"2013 XZ8",
"2014 QA43"]

data_path = "data"
fname_save = "df_results.csv"

df_results = pd.DataFrame()

for name in name_list:

    # if name!="Bienor":
    #     continue

    # get the mpc_number and orbital elements (could also get galactic latitudes)
    obj = Horizons(id=name)
    el = obj.elements()
    # print(el)
    # print(list(el))
    des = Names.parse_asteroid(str(el['targetname'][0]))
    try:
        mpc_number = des["number"]
    except:
        mpc_number=False
    # print(name,mpc_number)

    # create the object dataframe
    df_obj = pd.DataFrame()
    df_obj["name"] = [name]
    df_obj["mpc_number"] = [mpc_number]
    df_obj["H_abs_mag"] = float(el["H"])
    df_obj["G_slope"] = float(el["G"])
    df_obj["a_semimajor_axis"] = float(el["a"])
    df_obj["e_eccentricity"] = float(el["e"])
    df_obj["i_inclination_deg"] = float(el["incl"])
    print(df_obj)

    # load and set up the photometry dataframe
    if mpc_number:
        load_name = "{} {}".format(mpc_number,name)
    else:
        load_name = name
    df_o = pd.read_csv("{}/{}_ATLAS_Forced_Photometry_o_filter.csv".format(data_path,load_name),index_col=None)
    df_c = pd.read_csv("{}/{}_ATLAS_Forced_Photometry_c_filter.csv".format(data_path,load_name),index_col=None)
    df_o["filter"] = ["o"]*len(df_o)
    df_c["filter"] = ["c"]*len(df_c)
    data_all_filt = df_o.append(df_c)
    data_all_filt = data_all_filt.rename(columns={"MJD":"mjd","alpha":"phase_angle","reduced_dmag":"merr"})
    print(data_all_filt)

    JPL_file = "{}/{}_JPL_Horizons.csv".format(data_path,load_name)
    if os.path.isfile(JPL_file):
        print("load {}".format(JPL_file))
        df_eph = pd.read_csv(JPL_file,index_col=0)
    else:
        t1 = Time(int(np.amin(data_all_filt["mjd"])), format='mjd').iso
        t2 = Time(int(np.amax(data_all_filt["mjd"])), format='mjd').iso
        obj = Horizons(id=name,epochs={"start":t1,"stop":t2,"step":"1d"})
        eph = obj.ephemerides()
        df_eph = eph.to_pandas()
        df_eph = df_eph.rename(columns={"GlxLat":"galactic_latitude","elong":"sun_obs_target_angle"})
        df_eph["night(mjd)"] = df_eph["datetime_jd"] - 2400000.5
        df_eph = df_eph[["night(mjd)","sun_obs_target_angle","galactic_latitude"]]
        df_eph.to_csv(JPL_file)
    print(df_eph)

    # combine the JPL horizons nightly elongation and galactic latitude with forced photometry
    data_all_filt["night(mjd)"] = np.floor(data_all_filt["mjd"])
    data_all_filt = data_all_filt.merge(df_eph,on="night(mjd)")
    # option to drop/pass galactic latitude if you want to drop points near plane
    data_all_filt = data_all_filt.drop("galactic_latitude",1)
    print(data_all_filt)
    # exit()

    fit = sbpy_phase_fit.phase_fit(mpc_number,name,
    plot_fig_flag=True,save_fig_flag=True,
    # show_fig_flag=True,
    save_path="figs",
    # start_date=start_date, end_date=end_date,
    # model_list=["HG"],
    # filter_list=["c"],
    mag_diff_flag=True,save_file_suffix="_mag_diff",
    connection=False
    )
    df=fit.calculate_forced_phot(data_all_filt,df_obj)
    print(df.iloc[0].to_string())
    if name:
        objid="_".join(name.split(" "))
    else:
        objid = int(mpc_number)
    print("open figs/sbpy_phase_fit_{}_HG*".format(objid))

    df_results = df_results.append(df)

    # break

df_results.to_csv(fname_save)
