"""
Calls the atlas_SQL_query function to get an object's observations from ATLAS
Uses the solar elongation to analyse the apparitions of an object.
Looks for drops in solar elongation angle to define epochs.

N.B. will save queried data and try reload, these files may need updated over time!
Delete the loaded files to get fresh ones
"""

import pandas as pd
from optparse import OptionParser
import time
import numpy as np
from astropy.time import Time
from astroquery.jplhorizons import Horizons
import os
from scipy.signal import argrelextrema

# importing our classes depends on whether this file as a script or module
if __name__ == "__main__": # import other modules as a script
    import atlas_SQL_query_df
    import atlas_database_connection
else: # import as a package
    import calculate_phase.atlas_SQL_query_df as atlas_SQL_query_df
    import calculate_phase.atlas_database_connection as atlas_database_connection

class solar_apparitions():
    # fixed parameters
    loc = "T05" # which telescope to use for Horizons query, T05 or T08?
    mjd_jd = 2400000.5 # conversion from MJD to JD
    # G = 6.693e-11 # Newton's Gravitational Constant, m3⋅kg−1⋅s−2
    # M_sun = # Mass of the sun
    year_days = 365.25

    def __init__(self,
        mpc_number=False,
        name=False,
        df_data=None,
        reload=False,
        save_path=False,
        elong_mod = 360,
        data_load_path=None,
        eph_load_path=None):

        # path to save a figure
        self.save_path = save_path

        # object identifiers. Use mpc_number if available to avoid name confusion (e.g. asteroid vs comet Fitzsimmons)
        if mpc_number:
            self.obj_id = mpc_number
            self.obj_id_save = mpc_number
        else:
            self.obj_id = name
            self.obj_id_save = "_".join(name.split())

        # set up data path names
        self.data_load_path = data_load_path
        self.eph_load_path = eph_load_path
        self.data_fname = "{}/df_data_{}.csv".format(data_load_path,self.obj_id_save)
        self.eph_fname = "{}/df_eph_{}.csv".format(eph_load_path,self.obj_id_save)
        # print(self.data_fname)
        # print(self.eph_fname)

        # Option to delete existing data files
        if reload:
            for f in [self.data_fname,self.eph_fname]:
                try:
                    os.remove(f)
                    print("delete {}".format(f))
                except:
                    continue

        # Check if dataframe has been passed, otherwise load
        if df_data is None:

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
                if data_load_path is not None:
                    print("save {}".format(self.data_fname))
                    df_data.to_csv(self.data_fname)
        else:
            print("use passed df_data")

        # print(df_data)
        print(np.median(df_data["reduced_mag"]))
        print(np.median(df_data["merr"]))

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

        # initialise ephemerides dataframe
        self.df_eph = None

    def solar_elongation_JPL(self,JPL_step='1d'):

        # Get the JPL Horizons ephemerides to search for the turning points of solar elongation
        # This is a relatively slow check, only do this for NEAs
        # TRY LOAD AN EXISTING FILE
        if os.path.isfile(self.eph_fname):
            print("load Horizons")
            df_eph = pd.read_csv(self.eph_fname,index_col=0)
        else:
            print("query Horizons")
            t1=Time(int(np.amin(self.df_data["mjd"])),format="mjd")
            t2=Time(int(np.amax(self.df_data["mjd"])),format="mjd")

            if t1==t2:
                print("dates equal, just use upper and lower mjds")
                # add the first and last mjds in the dataset to fully bound all apparitions
                turning_points = np.array([])
                turning_points = np.insert(turning_points,0,self.xdata[0])
                turning_points = np.insert(turning_points,len(turning_points),self.xdata[-1])
                return np.array(turning_points)

            epoch_list = {'start':t1.iso, 'stop':t2.iso, 'step':JPL_step} # a range of epochs in Horizons format is FAST!
            # print(epoch_list)
            obj = Horizons(id=self.obj_id, location=self.loc, epochs=epoch_list)
            eph = obj.ephemerides()
            df_eph = eph.to_pandas()
            df_eph["mjd"] = df_eph["datetime_jd"] - self.mjd_jd
            # save df
            if self.eph_load_path is not None:
                print("save {}".format(self.eph_fname))
                df_eph.to_csv(self.eph_fname)

        df_eph = df_eph.sort_values("mjd")
        self.df_eph = df_eph

        # find the minima from JPL data
        elon_min_mjd = np.array(df_eph.iloc[argrelextrema(np.array(df_eph["elong"]), np.less)[0]]["mjd"])

        # NOTE that searching for minima will not work when less than one apparition has been observed
        # if time_span<synodic_period/2.0:
            # elon_min_mjd = np.array([])

        # select by the exact JPL ephemerides
        turning_points = elon_min_mjd

        # add the first and last mjds in the dataset to fully bound all apparitions
        turning_points = np.insert(turning_points,0,self.xdata[0])
        turning_points = np.insert(turning_points,len(turning_points),self.xdata[-1])
        # print(turning_points)
        # N_app = len(turning_points)-1
        turning_points = np.array(turning_points)

        # ONLY COUNT APPARITIONS WITH DATA! see e.g. Kridsadaporn
        # MAKE THIS FASTER, LIST comprehension?
        n_data_app = []
        for i in range(1,len(turning_points)):
            # print(i-1,i)
            n_data_app.append(len(self.xdata[(self.xdata>=turning_points[i-1]) &
            (self.xdata<=turning_points[i])])) # equal on upper and lower bounds to catch cases with one detection e.g. Erinleeryan
        n_data_app = np.array(n_data_app)
        # print(n_data_app)
        # print(turning_points[1:][n_data_app>0])
        turning_points = np.insert(turning_points[1:][n_data_app>0],0,turning_points[0])
        return turning_points

    def solar_elongation(self,sol_elong_diff,period=False):
        # find the turning points of the solar elongation data
        # by searching for drops in the solar_elongation angle (which is generally monotonic when %360)
        ydiff = np.diff(self.ydata)

        # we can check based on gaps in time as well as drops in elongation
        if period:
            xdiff = np.diff(self.xdata)
            synodic_period = (period / ((period / 1.0) - 1.0))*self.year_days # CHECK THIS https://en.wikipedia.org/wiki/Elongation_(astronomy)
            print("synodic period = {}".format(synodic_period))
            # mask = (ydiff<sol_elong_diff) | (xdiff>synodic_period)
            # mask = (ydiff<sol_elong_diff) & (xdiff>synodic_period/2.0)

            # mask = (ydiff<sol_elong_diff) | (xdiff>synodic_period/2.0)
            mask = ((ydiff<sol_elong_diff) & (xdiff>14.0)) | (xdiff>synodic_period/2.0)

        else:
            mask = (ydiff<sol_elong_diff)

        # select only drops of at least a certain size
        turning_points = self.xdata[1:][mask]

        # # IF PERIOD
        # # reject any turning points that have a time gap shorter than period
        # if period and (len(turning_points)>0):
        #     print(turning_points)
        #     print(np.diff(turning_points))
        #     mask_xdiff = list(np.diff(turning_points)<synodic_period)
        #     # print(mask_xdiff)
        #     indices = np.array(range(len(turning_points)))
        #     # print(indices)
        #     # print(indices[1:])
        #     drop = indices[1:][mask_xdiff]
        #     print(drop)
        #     turning_points = np.delete(turning_points,drop)

        # add the first and last mjds in the dataset to fully bound all apparitions
        turning_points = np.insert(turning_points,0,self.xdata[0])
        turning_points = np.insert(turning_points,len(turning_points),self.xdata[-1])
        print(turning_points)
        # N_app = len(turning_points)-1


        return np.array(turning_points)

    def plot_solar_elongation(self,turning_points=[],label=None):

        if self.save_path:
            # non-interactive
            import matplotlib
            matplotlib.use('Agg')

        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        fig = plt.figure()
        gs = gridspec.GridSpec(1, 1)
        ax1 = plt.subplot(gs[0,0])

        # plot the solar elongation from rockAtlas
        ax1.scatter(self.df_data["mjd"],np.absolute(self.df_data["sun_obs_target_angle"]),
        c="k",s=5,
        label="solar elongation")

        # plot the modded elongation angle
        ax1.scatter(self.xdata,self.ydata,c="b",s=5,label=self.ydata_lab)
        ax1.set_xlabel("mjd")
        ax1.set_ylabel("angle(degrees)")

        # plot the actual JPL ephemerides solar elongation if available
        if self.df_eph is not None:
            ax1.scatter(self.df_eph["mjd"],self.df_eph["elong"],
            c="r",s=5,alpha=0.5,zorder=0,
            label="solar elongation JPL")

        # mark the edges of the epochs with vertical lines
        for x in turning_points:
            ax1.axvline(x,c="r")

        # add a blank line for the label
        if label is not None:
            ax1.axvline(np.nan,c="r",label = label)

        # ax1.legend()

        title = "{}_{}".format("solar_apparitions",self.obj_id_save)
        fig.suptitle(title)
        plt.tight_layout()

        # ax2 = ax1.twinx()
        # ax2.scatter(self.df_data["mjd"],self.df_data["reduced_mag"])
        # ax2.set_ylabel("reduced_mag")
        #
        # # Try reorder axes
        # ax1.set_zorder(ax2.get_zorder()+1)
        # ax1.patch.set_visible(False)

        if self.save_path:
            fname="{}/{}.png".format(self.save_path,title)
            print(fname)
            plt.savefig(fname, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

        return

if __name__ == "__main__":

    # mpc_number=4986
    # mpc_number=3103
    # name=False
    mpc_number = False

    # objects where ydiff does not work
    name = "Asher"
    period = 2.58 # yrs
    # name = "2002 CW205"
    # period = 3.59
    # name = "2001 UG"
    # period = 3.26
    # name = "2008 RQ98"
    # period = 4.19
    # name = "2002 XQ7"
    # period = 4.27
    # name = "2003 OH12"
    # period = 5.56
    # name = "2001 TG99"
    # period = 7.94

    # Load/query data
    # sol = solar_apparitions(mpc_number=mpc_number,name=name)
    data_load_path="/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/results_analysis/obs"
    eph_load_path="/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/solar_apparitions_files"
    sol = solar_apparitions(mpc_number=mpc_number,name=name,data_load_path=data_load_path,eph_load_path=eph_load_path)
    tp1 = sol.solar_elongation_JPL(JPL_step="7d")
    print(tp1)
    tp2 = sol.solar_elongation(-1.0)
    # tp2 = sol.solar_elongation(-1.0,period=period)
    print(tp2)
    # sol.plot_solar_elongation(tp1)
    sol.plot_solar_elongation(tp2)
    print(len(tp1)-1,len(tp2)-1)

    # # Pass data directly, otherwise the code will load/query
    # data_fname = "/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/results_analysis/obs/df_data_4986.csv"
    # df_data = pd.read_csv(data_fname,index_col=0)
    # sol = solar_apparitions(mpc_number=mpc_number,df_data=df_data)
