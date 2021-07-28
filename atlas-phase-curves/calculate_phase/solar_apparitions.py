"""
Calls the atlas_SQL_query function to get an object's observations from ATLAS
Uses the solar elongation to analyse the apparitions of an object.
Looks for drops in solar elongation angle to define epochs.

N.B. will save queried data and try reload, these files may need updated over time!
Delete the loaded files to get fresh ones
"""

import pandas as pd
from calculate_phase import atlas_SQL_query_df
from calculate_phase import atlas_database_connection
from optparse import OptionParser
import time
import numpy as np
from astropy.time import Time
from astroquery.jplhorizons import Horizons
import os
from scipy.signal import argrelextrema

class solar_apparitions():
    # fixed parameters
    loc = "T05" # which telescope to use for Horizons query, T05 or T08?
    mjd_jd = 2400000.5 # conversion from MJD to JD

    def __init__(self,
        mpc_number=False,
        name=False,
        reload=False,
        save_path=False,
        elong_mod = 360,
        data_load_path="/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/results_analysis/obs",
        eph_load_path="solar_apparitions_files"):

        # set up data path names
        if name:
            self.obj_id = name
            self.obj_id_save = "_".join(name.split())
        else:
            self.obj_id = mpc_number
            self.obj_id_save = mpc_number
        self.data_fname = "{}/df_data_{}.csv".format(data_load_path,self.obj_id_save)
        self.eph_fname = "{}/df_eph_{}.csv".format(eph_load_path,self.obj_id_save)
        print(self.data_fname)
        print(self.eph_fname)

        # path to save a figure
        self.save_path = save_path

        # delete existing data to force refresh
        if reload:
            for f in [self.data_fname,self.eph_fname]:
                try:
                    os.remove(f)
                    print("delete {}".format(f))
                except:
                    continue

        # Try load an existing data file
        if os.path.isfile(self.data_fname):
            print("load data")
            df_data = pd.read_csv(self.data_fname,index_col=0)
        else:
            print("query data")
            # connect to database
            cnx=atlas_database_connection.database_connection().connect()

            # get orbit_id
            orbid=atlas_SQL_query_df.get_orb_elements_id(cnx,mpc_number,name)

            # do the query
            df_data=atlas_SQL_query_df.atlas_SQL_query_orbid_expname(cnx,orbid)

            # save df
            df_data.to_csv(data_fname)

        # ADD OPTION TO PASS DATAFRAME INSTEAD OF LOADING

        # drop any nans and sort by date
        df_data = df_data.dropna()
        df_data = df_data.sort_values(['mjd'])
        self.df_data = df_data

        # Solar elongation retrieved from rockAtlas
        # N.B. rockAtlas solar elongation angle cycles from -180 to 180 deg, use mod to get 0 to 360
        # or abs to force 0 -> 180
        self.xdata=np.array(df_data["mjd"])
        if elong_mod:
            self.ydata=np.array(df_data["sun_obs_target_angle"])%elong_mod
            self.ydata_lab="solar elongation % {}".format(elong_mod)
        else:
            self.ydata=np.array(df_data["sun_obs_target_angle"])
            self.ydata_lab="solar elongation"

    def solar_elongation_JPL(self):

        # Get the JPL Horizons ephemerides to search for the turning points of solar elongation
        # This is a relatively slow check, only do this for NEAs
        # TRY LOAD AN EXISTING FILE
        if os.path.isfile(self.eph_fname):
            print("load Horizons")
            df_eph = pd.read_csv(self.eph_fname,index_col=0)
        else:
            print("query Horizons")
            t1=Time(int(np.amin(df_data["mjd"])),format="mjd")
            t2=Time(int(np.amax(df_data["mjd"])),format="mjd")
            epoch_list = {'start':t1.iso, 'stop':t2.iso, 'step':'1d'} # a range of epochs in Horizons format is FAST!
            obj = Horizons(id=self.obj_id, location=self.loc, epochs=self.epoch_list)
            eph = obj.ephemerides()
            df_eph = eph.to_pandas()
            df_eph["mjd"] = df_eph["datetime_jd"] - self.mjd_jd
            # save df
            df_eph.to_csv(eph_fname)

        df_eph = df_eph.sort_values("mjd")

        # find the minima from JPL data
        elon_min_mjd = np.array(df_eph.iloc[argrelextrema(np.array(df_eph["elong"]), np.less)[0]]["mjd"])

        # select by the exact JPL ephemerides
        turning_points = elon_min_mjd

        # add the first and last mjds in the dataset to fully bound all apparitions
        turning_points = np.insert(turning_points,0,self.xdata[0])
        turning_points = np.insert(turning_points,len(turning_points),self.xdata[-1])
        # print(turning_points)
        # N_app = len(turning_points)-1
        return np.array(turning_points)

    def solar_elongation(self,sol_elong_diff):
        # find the turning points of the solar elongation data
        # by searching for drops in the solar_elongation angle (which is generally monotonic when %360)
        ydiff = np.diff(self.ydata)

        # select only drops of at least a certain size
        turning_points = self.xdata[1:][ydiff<sol_elong_diff]

        # add the first and last mjds in the dataset to fully bound all apparitions
        turning_points = np.insert(turning_points,0,self.xdata[0])
        turning_points = np.insert(turning_points,len(turning_points),self.xdata[-1])
        # print(turning_points)
        # N_app = len(turning_points)-1
        return np.array(turning_points)

    def plot_solar_elongation(self,turning_points=[]):

        if self.save_path:
            # non-interactive
            import matplotlib
            matplotlib.use('Agg')

        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        fig = plt.figure()
        gs = gridspec.GridSpec(1, 1)
        ax1 = plt.subplot(gs[0,0])

        ax2 = ax1.twinx()
        ax2.scatter(self.df_data["mjd"],self.df_data["reduced_mag"])
        ax2.set_ylabel("reduced_mag")

        ax1.scatter(self.xdata,self.ydata,c="k",s=5,label=self.ydata_lab)
        ax1.set_xlabel("mjd")
        ax1.set_ylabel("angle(degrees)")

        ax1.legend()

        for x in turning_points:
            ax1.axvline(x,c="r")

        title = "{}_{}".format("solar_apparitions",self.obj_id_save)
        fig.suptitle(title)
        plt.tight_layout()

        if self.save_path:
            fname="{}/{}.png".format(self.save_path,title)
            print(fname)
            plt.savefig(fname, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

        return

if __name__ == "__main__":

    sol = solar_apparitions(mpc_number=4986)
    tp1 = sol.solar_elongation_JPL()
    tp2 = sol.solar_elongation(-1.0)
    # sol.plot_solar_elongation(tp1)
    sol.plot_solar_elongation(tp2)
