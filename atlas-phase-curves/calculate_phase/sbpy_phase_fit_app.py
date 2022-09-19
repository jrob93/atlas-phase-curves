"""
Module to perform ATLAS asteroid phase fits with various models and parameters.
Fit each apparition separately
"""

import numpy as np
import astropy.units as u
from astropy.modeling.fitting import LevMarLSQFitter
from sbpy.photometry import HG, HG1G2, HG12, HG12_Pen16
from sbpy.photometry import LinearPhaseFunc
import pandas as pd
import os
import mysql.connector
import datetime
from pathlib import Path
import logging
from astropy.utils.exceptions import AstropyUserWarning
from astropy.utils.exceptions import AstropyWarning
import warnings
from sklearn.preprocessing import StandardScaler

# # importing our classes depends on whether this file as a script or module
# if __name__ == "__main__": # import other modules as a script
#     from atlas_SQL_query_df import atlas_SQL_query
#     from atlas_database_connection import database_connection
# else: # import as a package
#     from calculate_phase.atlas_SQL_query_df import atlas_SQL_query
#     from calculate_phase.atlas_database_connection import database_connection

from calculate_phase.atlas_SQL_query_df import atlas_SQL_query,get_orb_elements_id,get_unique_ids,atlas_SQL_query_orbid_expname
from calculate_phase.atlas_database_connection import database_connection
import calculate_phase.solar_apparitions as sa

class phase_fit():

    # default fixed values
    low_alpha_cut=5.0*u.deg # we want to quantify how many data points are fit at low phase angles, alpha < low_alpha_cut
    param_converge_check=0.01 # the model is fit until the change in parameters (e.g. H and G) is less than param_converge_check (or max iterations is reached)
    max_iters=30 # maximum number of attempts at fitting and cutting
    std=2 # standard deviation of the sigma data clip
    mag_err_threshold = 0.1 # limit for the error of "good" data, we record N_mag_err number of data points with error < mag_err_threshold
    mag_err_small = 0.005 # we discount observations with error less than this
    gal_lat_cut=10 # galatic latitude cut in degrees
    # mag_med_cut=2 # initial magnitude difference cut on initial HG model
    mag_med_cut=5 # initial magnitude difference cut on initial HG model - most lcdb lightcurves have (peak to peak) amplitude < 2.5
    phase_lin_min=5 # minimum phase angle for linear fit - MAKE OPTIONAL?
    phase_lin_max=25 # maximum phase angle for linear fit
    orbfit_sep_cut=1.0 # maximum allowed orbfit separation for matching dophot_photometry (arcsec)
    min_app_data=10 # minimum number of data points to consider when fitting individual apparitions
    lim_mag = 17.5 # minimum limiting_magnitude of an atlas exposure to be used

    # set the clipping method for the phase_fit class
    def data_clip_sigma(self,data,data_predict,low=std,high=std):
        """ function to cut outliers by standard deviation, i.e. sigma clipping
        returns the mask of data points to cut"""
        std=np.std(data) # !!!SHOULD THIS BE THE STD OF DATA-DATA PREDICT???
        print("std = {}".format(std))
        print("std = {}".format(np.std(data-data_predict)))
        print("low = {}, high = {}".format(low,high))
        clip_mask=((data < (data_predict-(std*low))) | (data > (data_predict+(std*high))))
        return clip_mask

    def data_clip_diff(self,data,data_predict,diff=1):
        """cut outliers by diff (this function doesn't like astropy units, use np arrays)
        returns the mask of data points to cut"""
        x=np.array(np.absolute(data_predict-data))
        clip_mask=(x>diff)
        return clip_mask

    # define which clipping function to use
    # SET THIS IN __init__?
    # data_clip_func=data_clip_sigma
    # clip_label="{}-sigma_clip".format(std)

    # set up the astropy/scipy fitter and sbpy models
    # https://sbpy.readthedocs.io/en/latest/sbpy/photometry.html#disk-integrated-phase-function-models
    fitter = LevMarLSQFitter()

    # These are all the sbpy models that can be selected (use model list in __init__)
    all_models={
    "HG":{"model_name_short":"_B89","model_parameters":["H","G"],"model_function":HG()},
    "HG1G2":{"model_name_short":"_3M10","model_parameters":["H","G1","G2"],"model_function":HG1G2()},
    "HG12":{"model_name_short":"_2M10","model_parameters":["H","G12"],"model_function":HG12()},
    "HG12_Pen16":{"model_name_short":"_P16","model_parameters":["H","G12"],"model_function":HG12_Pen16()},
    "LinearPhaseFunc":{"model_name_short":"_Lin","model_parameters":["H","S"],"model_function":LinearPhaseFunc(H=15,S=0.04)}
    }

    # IS THERE A FASTER WAY, E.G. ONLY LOAD THIS ONCE AND IF PUSH FIT IS TRUE?
    # load the columns used to make the db table (might need to update the path)
    # fname_path = Path(__file__).parent / "../create_table/atlas_objects_fields.txt"
    fname_path = "{}/{}".format(Path(__file__).parent,"../create_table/atlas_objects_app_fields.txt")
    with open(fname_path,"r") as f:
        db_columns=f.readlines()
    db_columns=[d.rstrip() for d in db_columns]

    def __init__(self,
        mpc_number=False, # mpc number of object to be fit (if available)
        name=False, # name/designation of object
        save_path=".", # location to save figures
        save_file_suffix="", # add a suffix to saved files to customise their filenames if required
        save_file_type="png", # file type for saving figures
        push_fit_flag=False,plot_fig_flag=False,show_fig_flag=False,save_fig_flag=False,hide_warning_flag=False, # flags controlling plotting/saving of figures
        plot_elong_fig=False, # flag to make solar elongation plot
        start_date=False,end_date=False, # set a start date and end date to control selection of observations
        mag_diff_flag=True, # Flag to perform an initial cut of observations based on magnitude difference from expected values. DEFAULT THIS TO BE TRUE?
        model_list=["HG", "HG1G2", "HG12", "HG12_Pen16"], # Which models to fit. ADD LinearPhaseFunc here as default?
        filter_list=["o","c"], # which filters to fit
        tab_name="atlas_phase_fits_app", # name of the sql table to save results
        save_data_flag=False, # save the data that was used in the fit?
        connection=True # flag to control if we establish a connection to the sql database or not
        ):

        # set up the class
        # define object and some flags

        # the object will have EITHER both mpc number and name OR just a name
        self.mpc_number=mpc_number
        self.name=name

        # set the either the name or mpc_number for naming output files - N.B. what if name is passed but object does have mpc_number?
        if name:
            self.file_identifier="_".join((self.name).split())
        else:
            self.file_identifier=self.mpc_number

        self.push_fit=push_fit_flag # flag to push fit to database
        self.plot_fig=plot_fig_flag # flag to generate plot for each object
        self.plot_elong_fig=plot_elong_fig
        self.show_fig=show_fig_flag # flag to display interactive plot
        self.save_fig=save_fig_flag # flag to save the figure
        self.utc_date_now=datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S") # time at which the script is run (UTC)
        self.save_path=save_path # where to save figures
        self.save_file_suffix=save_file_suffix # option to add a suffix to the png of plot saved
        self.save_file_type=save_file_type # option to choose type of file saved, e.g. png or pdf
        self.mag_diff_flag=mag_diff_flag # flag to perform an additional mag diff cut on the initial HG model
        self.save_data_flag=save_data_flag

        self.selected_models= {key: value for key, value in self.all_models.items() if key in model_list} # create the dictionary of models to use

        self.start_date=start_date # only use observations after this date (mjd)
        self.end_date=end_date # only use observations before this date (mjd)
        # self.date_hard_cap=? # set a fixed date beyond which we do not consider data - SET BY end_date
        self.filters=filter_list # list of filters to use: ["o"], ["c"], ["o","c"]

        # we do not want to set up a conenction if we are supplying the observations directly, i.e. using calculate_forced_phot
        if connection:
            # cnx1 is the connection where we retrieve the data from
            self.cnx1=database_connection().connect()
            self.cursor1 = self.cnx1.cursor()

            # cnx2 is the connection where we push the fit data to
            self.cnx2 = self.cnx1
            self.cursor2 = self.cursor1

        self.tab_name=tab_name # table name within the database

        # create an empty pandas series to hold all object info
        self.df_obj_datafit=pd.DataFrame(data=[],columns=self.db_columns) # create an empty series to store all values for the object: metadata and fit data
        print(self.df_obj_datafit)

        # all objects in the db have a unique orbital_elements_id - retrieve this from db using the mpc_number or name
        try:
            self.orbital_elements_id=get_orb_elements_id(self.cnx1,self.mpc_number,self.name)
        except:
            # print("Cannot find object {}, {} in database".format(self.mpc_number,self.name))
            self.orbital_elements_id = None
            # return

        # option to suppress warnings being shown to the screen
        if hide_warning_flag==1:
            print("hide warnings")
            warnings.filterwarnings('ignore')

        # set up a log file - This will append to an existing log file of the same name, delete the file if you need a new log
        logging.basicConfig(filename='{}.log'.format(os.path.basename(__file__).split('.')[0]), level=logging.INFO)

    def get_obj_data(self,cnx,orbid,t_start=False,t_end=False):
        """load data to be fitted, loads both filters (o & c)"""

        # query the database for photometry
        data_all_filt=atlas_SQL_query_orbid_expname(cnx,orbid)

        N_data1 = len(data_all_filt)
        print("data before date cut = {}".format(N_data1))
        # cut by date range. Note that different filters may have different dates/epochs!
        if t_start:
            data_all_filt=data_all_filt[data_all_filt['mjd']>float(t_start)]
        if t_end:
            data_all_filt=data_all_filt[data_all_filt['mjd']<float(t_end)]
        # else: # if t_end is not defined default to the hard cap for date
        #     data_all_filt=data_all_filt[data_all_filt['mjd']<self.date_hard_cap]

        # Number of data points outside of date range
        N_data2 = len(data_all_filt)
        print("data after date cut = {}".format(N_data2))
        N_data_date = N_data1-N_data2
        print("CUT N_data_date = {}".format(N_data_date))

        # DROP ALL ROWS WITH A NAN for m,merr,mjd,phase_angle
        N_data3 = len(data_all_filt)
        print("data before nan cut = {}".format(N_data3))
        # data_all_filt=data_all_filt.dropna()
        data_all_filt=data_all_filt.dropna(subset = ["m", "merr", "mjd","phase_angle"])
        N_data4 = len(data_all_filt)
        print("data after nan cut = {}".format(N_data4))
        N_data_nan = N_data3 - N_data4
        print("CUT N_data_nan = {}".format(N_data_nan))

        return data_all_filt

    # def get_obj_HG(self,cnx,orbid):
    #     # get the astorb H and G (B89)
    #     # qry_HG="select G_slope, H_abs_mag from orbital_elements where primaryId='{}';".format(self.mpc_number)
    #     # qry_HG="SELECT G_slope, H_abs_mag FROM orbital_elements WHERE mpc_number='{}';".format(mpc_num)
    #     qry_HG="SELECT G_slope, H_abs_mag FROM orbital_elements WHERE orbital_elements_id='{}';".format(orbid)
    #     # print(qry_HG)
    #     df_HG=pd.read_sql_query(qry_HG,cnx)
    #     return df_HG

    def get_obj_meta(self,cnx,orbid,mpc_num,name):
        """ Get the object metadata.
        Might be out of date if rockatlas isn't running or hasn't has CALL update_atlas_objects in a while...?"""

        # use either the mpc_number or object name to query atlas_objects AND orbital_elements tables
        if mpc_num:
            qry_obj1=u"""SELECT
            a.dateLastModified,
            a.detection_count,
            a.detection_count_c,
            a.detection_count_o,
            a.last_detection_mjd,
            a.last_photometry_update_date_c,
            a.last_photometry_update_date_o,
            a.mpc_number,
            a.name,
            a.orbital_elements_id,
            a.primaryId,
            a.updated,
            o.G_slope,
            o.H_abs_mag,
            o.a_semimajor_axis,
            o.e_eccentricity,
            o.i_inclination_deg
            FROM atlas_objects a, orbital_elements o WHERE a.orbital_elements_id=%(orbid)s AND o.mpc_number=%(mpc_num)s;
            """ % locals()
        else:
            qry_obj1=u"""SELECT
            a.dateLastModified,
            a.detection_count,
            a.detection_count_c,
            a.detection_count_o,
            a.last_detection_mjd,
            a.last_photometry_update_date_c,
            a.last_photometry_update_date_o,
            a.mpc_number,
            a.name,
            a.orbital_elements_id,
            a.primaryId,
            a.updated,
            o.G_slope,
            o.H_abs_mag,
            o.a_semimajor_axis,
            o.e_eccentricity,
            o.i_inclination_deg
            FROM atlas_objects a, orbital_elements o WHERE a.orbital_elements_id=%(orbid)s AND o.name="%(name)s";
            """ % locals()
        # print(qry_obj1)

        # send the query to the rockAtlas db
        df_obj=pd.read_sql_query(qry_obj1,cnx)
        # print(df_obj[["primaryId","orbital_elements_id","mpc_number","name"]])

        # fix the date strings returned by mysql
        dates=["dateLastModified","last_photometry_update_date_c","last_photometry_update_date_o"]
        for d in dates:
            if df_obj[d].iloc[0] is not None:
                df_obj.loc[0,d]="'{}'".format(str(df_obj[d].iloc[0]))

        # for some reason None needs changed to Null?
        df_obj=df_obj.fillna(value="NULL")

        return df_obj

    def push_obj_db(self,df_obj):
        """ This function will take the dataframe from the fit and only push the columns in self.db_columns
        These are the columns in the db, other columns will break the sql query.
        G_slope, H_abs_mag and the linear phase fits will not be pushed
        """

        print("push fit to db")

        tab_name=self.tab_name
        utc_date_now=self.utc_date_now
        mpc_number=self.mpc_number

        cols=self.db_columns
        vals=df_obj[cols].iloc[0]
        vals=np.array(vals).astype(str)

        # for c,v in zip(cols,vals):
        #     print(c,v)

        N_cols=len(cols)
        N_vals=len(vals)

        # create a qry that sql can handle!
        col_vals_update=""
        for i in range(N_cols):
            # print("update {} {}".format(cols[i],vals[i]))
            # get rid of if statement if possible? better coding? list comprehension?

            # catch any values that are nan
            if vals[i].lower()=="nan": # mysql doesn't understand nan, needs NULL
                vals[i]="NULL"

            # catch any names/dates that will need fixed
            if cols[i] in ["name","phase_curve_refresh_date_o","phase_curve_refresh_date_c","dateLastModified","last_photometry_update_date_c","last_photometry_update_date_o"]: # these fields are strings and need quotation marks
                if vals[i]!="NULL":
                    col_vals_update+="{}=\"{}\",".format(cols[i],vals[i])
                    vals[i]="\"{}\"".format(vals[i])
                else:
                    col_vals_update+="{}={},".format(cols[i],vals[i])
            else:
                col_vals_update+="{}={},".format(cols[i],vals[i])

        col_vals_update=col_vals_update[:-1] # drop the last comma

        # sanity check to make sure no values are dropped
        N_col_vals_update=len(col_vals_update.split(","))
        if N_col_vals_update!=N_vals or N_col_vals_update!=N_cols:
            print("check lengths: {} {} {}".format(len(col_vals_update.split(",")),len(vals),len(cols)))
            warning_message="{} - {} - SQL query lengths do not add up".format(self.mpc_number,self.name)
            print(warning_message)
            logging.warning(warning_message)

        qry=u"""INSERT INTO {} ({}) values ({}) ON DUPLICATE KEY UPDATE {};""".format(self.tab_name,",".join(cols), ",".join(vals), col_vals_update)
        # print(qry)

        self.cursor2.execute(qry)
        self.cnx2.commit()

        return

# migrate all these long plotting functions above to a different file?

    def calculate(self):
        """calculate the phase curves on the phase_fit object"""

        if self.orbital_elements_id is None:
            print("object not in db, nothing to fit")
            return

        # get the object observation data, cutting on date range and dropping rows with nan values for certain columns
        dataAllFilt=self.get_obj_data(self.cnx1,self.orbital_elements_id,self.start_date,self.end_date)

        # get the object metadata and combine with the phase fit dataframe structure
        df_obj_main=self.get_obj_meta(self.cnx1,self.orbital_elements_id,self.mpc_number,self.name)
        print(df_obj_main)
        d1=df_obj_main
        d2=self.df_obj_datafit
        df_obj_main=d2.append(d1) # add the values to the df

        # retrieve astorb H and G values for the predicted fit
        G_slope=float(df_obj_main.iloc[0]['G_slope'])
        H_abs_mag=float(df_obj_main.iloc[0]['H_abs_mag']) * u.mag
        print("G_slope = {}\nH_abs_mag = {}".format(G_slope,H_abs_mag))

        # set the mpc_number and name from df_obj
        df_obj_mpc_num = df_obj_main.iloc[0]['mpc_number']
        if df_obj_mpc_num!="NULL": # sql might lead to a NULL value overwriting the mpc_number which will interfere with cases of "if self.mpc_number"
            self.mpc_number=df_obj_main.iloc[0]['mpc_number']
        self.name=df_obj_main.iloc[0]['name']
        print(self.mpc_number,self.name)

        # get primaryId from atlas_objects
        primaryId = df_obj_main.iloc[0]["primaryId"]

        #-----
        # Find the solar apparitions from elongation
        # do this before cuts on observations so that apparitions in the galactic plane are recorded

        print(df_obj_main[["a_semimajor_axis","e_eccentricity","i_inclination_deg"]])
        q_perihelion = df_obj_main.iloc[0]["a_semimajor_axis"] * (1.0 - df_obj_main.iloc[0]["e_eccentricity"])
        epochs = []
        # if an object is an NEO, accurate JPL ephem query is required
        if q_perihelion<=1.3:
            print("NEO, JPL apparitions required")
            # Check if the JPL query is out of date. Note, if end_date is not set it is False and old queries will be wiped
            if np.amax(dataAllFilt["mjd"]) > self.end_date:
                # wipe the old JPL query if it exists, query will run again
                sol = sa.solar_apparitions(mpc_number = self.mpc_number, name = self.name, df_data = dataAllFilt,
                                            reload = True, eph_load_path = "NEO_JPL_eph")
            else:
                # use JPL apparition method, will load an existing eph file if available
                sol = sa.solar_apparitions(mpc_number = self.mpc_number, name = self.name, df_data = dataAllFilt,
                                            eph_load_path = "NEO_JPL_eph")

            epochs = sol.solar_elongation_JPL(JPL_step="7d")

        # Simple and fast apparition finder works for non-NEO objects
        else:
            # just use the simple apparition finder
            orbital_period_yrs = df_obj_main.iloc[0]["a_semimajor_axis"]**1.5
            sol = sa.solar_apparitions(mpc_number = self.mpc_number, name = self.name, df_data = dataAllFilt)
            epochs = sol.solar_elongation(-1.0,period = orbital_period_yrs)

        if len(epochs)==0:
            print("Error determining apparitions")
            logging.warning("{} - {} Error determining apparitions".format(self.mpc_number,self.name))

        # make the epoch plot?
        if self.plot_elong_fig:
            sol.save_path = self.save_path
            print(self.save_path)
            sol.plot_solar_elongation(epochs)

        print(epochs)
        N_app = len(epochs)-1 # number of apparitions detected in both filters
        df_obj_main["N_apparitions"]=N_app

        #-----
        # cuts on observations

        data_all_filt = dataAllFilt.copy()

        # Record number of detections in the filter AFTER nan and date cuts but BEFORE other cuts are made
        # This updates the detection_count value from rockatlas
        N_data_start = len(data_all_filt)
        df_obj_main["detection_count"]=N_data_start # update the number of detections in the df

        # cut starting data for this filter
        print("{} starting data".format(len(data_all_filt)))

        # drop measurements with large orbfit separation
        mask_orbfit = np.absolute(data_all_filt["orbfit_separation_arcsec"])>self.orbfit_sep_cut
        data_orbfit=data_all_filt[mask_orbfit]
        data_all_filt=data_all_filt[~mask_orbfit]
        print("{} after orbfit sep cut".format(len(data_all_filt)))
        # drop measurements with bright limiting_magnitude
        mask_lim_mag = np.absolute(data_all_filt["limiting_magnitude"])<self.lim_mag
        data_lim_mag=data_all_filt[mask_lim_mag]
        data_all_filt=data_all_filt[~mask_lim_mag]
        print("{} after lim mag cut".format(len(data_all_filt)))
        # drop any measurements with zero uncertainty
        mask_zero = data_all_filt['merr']==0
        data_zero_err=data_all_filt[mask_zero]
        data_all_filt=data_all_filt[~mask_zero]
        print("{} after zero error cut".format(len(data_all_filt)))
        # drop measurements with very small uncertainty
        mask_err = data_all_filt['merr']<self.mag_err_small
        data_small_err=data_all_filt[mask_err]
        data_all_filt=data_all_filt[~mask_err]
        print("{} after small error cut".format(len(data_all_filt)))
        # drop measurements near galactic plane
        mask_gal = np.absolute(data_all_filt["galactic_latitude"])<self.gal_lat_cut
        data_gal=data_all_filt[mask_gal]
        data_all_filt=data_all_filt[~mask_gal]
        print("{} after galactic plane cut".format(len(data_all_filt)))

        print("{} data after cuts".format(len(data_all_filt)))

        # RECORD THE NUMBER OF DATA POINTS THAT HAVE BEEN CUT
        N_data_orbfit = len(data_orbfit)
        N_data_lim_mag = len(data_lim_mag)
        N_data_zero_err = len(data_zero_err)
        N_data_small_err = len(data_small_err)
        N_data_gal = len(data_gal)
        data_all_cut = pd.concat([data_orbfit,data_lim_mag,data_zero_err,data_small_err,data_gal]) # store all data that is cut
        N_data_cut = N_data_orbfit + N_data_lim_mag + N_data_zero_err + N_data_small_err + N_data_gal
        print("CUT data_orbfit = {}\nCUT data_lim_mag = {}\nCUT data_zero_err = {}\nCUT data_small_err = {}\nCUT data_gal = {}".format(
        N_data_orbfit,N_data_lim_mag,N_data_zero_err,N_data_small_err,N_data_gal))
        print("TOTAL CUT N_data_cut = {}".format(N_data_cut))

        # if no data remains after loading, e.g. all data outside of date range, then nothing can be fit
        if len(data_all_filt)==0:
            print("no data to start with, cannot proceed")
            # push df_obj to the database
            if self.push_fit==True:
                self.push_obj_db(df_obj_main)
            return df_obj_main # return the fit df

        #------
        # select the primary epoch of the o filter data - best number of data points, phase angle range and number of obs at low phase
        # for each epoch count the number of data points, phase angle range and number of obs at <5 degrees phase
        N_data_epoch = []
        phase_range = []
        phase_low_alpha_cut = []

        for i in range(len(epochs)-1):
            mask = ((data_all_filt["mjd"]>=epochs[i]) & (data_all_filt["mjd"]<epochs[i+1])) & (data_all_filt["filter"]=="o")
            data = data_all_filt[mask]
            N_data_epoch.append(len(data))
            if len(data)==0:
                phase_range.append(np.nan)
                phase_low_alpha_cut.append(np.nan)
            else:
                phase_range.append(np.ptp(data["phase_angle"]))
                phase_low_alpha_cut.append(len(data[data["phase_angle"]<self.low_alpha_cut]))
            # print("N_data = {}".format(N_data_epoch[-1]))
            # print("phase range = {} degrees".format(phase_range[-1]))
            # print("N_data<5degrees = {}".format(phase_lt_5[-1]))

        # create dataframe to hold epoch info
        df_epoch = pd.DataFrame()
        df_epoch["epoch"] = np.arange(len(epochs[:-1]))
        df_epoch["mjd_start"] = epochs[:-1]
        df_epoch["N_data"] = N_data_epoch
        df_epoch["phase_angle_range"] = phase_range
        df_epoch["phase_angle_<low_alpha_cut"] = phase_low_alpha_cut

        # use a standard scaler to allow us to choose the best combination of number of data points, phase angle range and number of data points <5 deg
        X = df_epoch[["N_data","phase_angle_range","phase_angle_<low_alpha_cut"]]
        X = StandardScaler().fit_transform(X)
        X = np.sum(X, axis = 1)
        df_epoch["scaled"] = X

        # chose the largest not nan scaled parameter
        df_epoch.sort_values("scaled", inplace = True)
        epoch_prime_ind = int(df_epoch[~np.isnan(df_epoch["scaled"])].iloc[-1]["epoch"])
        print(df_epoch)
        print("best epoch: ",epoch_prime_ind,epochs[epoch_prime_ind], epochs[epoch_prime_ind+1])

        #-----

        df_obj_all = pd.DataFrame()

        # fit each epoch - in the "best" order as sorted above
        for i in range(len(df_epoch))[::-1]:

            # make a copy of the object dataframe frame containing all common metadata
            # we need a separate df for each apparition
            df_obj = df_obj_main.copy()

            print("\n")
            print(df_epoch.iloc[i])
            # continue
            epoch_ind = int(df_epoch.iloc[i]["epoch"])
            df_obj["app_start_mjd"]=epochs[epoch_ind]
            df_obj["app_ind"]=epoch_ind

            # set a unique primaryId for each row: original primaryId + apparition index
            df_obj["primaryId"] = str(primaryId) + str(epoch_ind)

            # define the data for this apparition
            mask_date = (((data_all_filt["mjd"]>=epochs[epoch_ind]) & (data_all_filt["mjd"]<epochs[epoch_ind+1])))

            # do a separate fit for data in each filter
            for filt in self.filters:

                df_obj["phase_curve_refresh_date_{}".format(filt)]=self.utc_date_now

                # do filter correction from V band (Heinze et al. 2020) - see also Erasmus et al 2020 for the c-o colours of S and C types (0.388 and 0.249 respectively)
                if filt=="o":
                    H_abs_mag+=(-0.332*u.mag)
                if filt=="c":
                    H_abs_mag+=(0.054*u.mag)

                # define the data for this apparition and filter
                mask_filt =  (data_all_filt["filter"]==filt)
                data = data_all_filt[mask_date & mask_filt]
                print(len(data))

                # record original number of datapoints in apparition - from original dataAllFilt dataframe
                mask = (((dataAllFilt["mjd"]>=epochs[epoch_ind]) & (dataAllFilt["mjd"]<epochs[epoch_ind+1])) & (dataAllFilt["filter"]==filt))
                N_data_app = len(dataAllFilt[mask])

                # Skip any epochs were all data points have been dropped - sometimes the first epoch in the dataframe might be empty, skip it
                if float(df_epoch.iloc[i]['N_data'])==0:
                    # the apparition should be pushed to the database, but with nans for missing columns
                    continue

                # define arrays for fitting
                alpha = np.array(data["phase_angle"]) * u.deg
                mag = np.array(data["reduced_mag"]) * u.mag
                mag_err = np.array(data["merr"]) * u.mag

                if self.mag_diff_flag:

                    # check there is sufficient data points to try fit
                    if len(data)<=2:
                        print("less data to fit than parameters")
                        continue

                    # For each epoch do an initial H fit and mag residual cut to remove remaining extreme outliers
                    model_HG = HG(H = H_abs_mag, G = G_slope)
                    model_HG.G.fixed = True # keep G fixed
                    print(model_HG)
                    # model_HG = self.fitter(model_HG, alpha, mag, weights=1.0/np.array(mag_err))
                    model_HG = self.fitter(model_HG, alpha, mag) # do not use weights to avoid fit being pulled by extremely low error outliers
                    print(model_HG)

                    # !!! Make this section neater, how to do two cuts in a row and store clipped data in one thing?

                    # do a first mag diff clip to get rid of the most ludicrous outliers
                    # e.g. the sigma clip for 1999 XN102 fails because of a really bad outliers
                    reduced_mag = model_HG(alpha)
                    mask_residual = self.data_clip_diff(np.array(mag),np.array(reduced_mag),diff=self.mag_med_cut)
                    data_residual = data[mask_residual] # drop these obs from data_all_filt as well
                    data = data[~mask_residual]
                    N_data_residual = len(data_residual)
                    print("mag diff residual cut")
                    print(len(data))
                    print(N_data_residual)
                    N_data_cut += N_data_residual

                    # define arrays for fitting (again)
                    alpha = np.array(data["phase_angle"]) * u.deg
                    mag = np.array(data["reduced_mag"]) * u.mag
                    mag_err = np.array(data["merr"]) * u.mag

                    # now we can do a sigma clip
                    reduced_mag = model_HG(alpha)
                    mask_residual = self.data_clip_sigma(np.array(mag),np.array(reduced_mag),low=self.std,high=self.std)
                    data_residual = pd.concat([data_residual,data[mask_residual]]) # drop these obs from data_all_filt as well
                    data_all_filt = data_all_filt.drop(data_residual.index)
                    data_all_cut = pd.concat([data_all_cut,data_residual]) # store the cut datapoints

                    data = data[~mask_residual]
                    N_data_residual = len(data_residual)

                    print("mag sig residual cut")
                    print(len(data))
                    print(N_data_residual)
                    N_data_cut += N_data_residual

                    # define arrays for fitting (again)
                    alpha = np.array(data["phase_angle"]) * u.deg
                    mag = np.array(data["reduced_mag"]) * u.mag
                    mag_err = np.array(data["merr"]) * u.mag

                    N_data_fit=len(data) # number of data points fit after clipping

                # retrieve metrics for the data
                N_data_fit=len(data) # number of data points fit after clipping
                alpha_min=np.amin(data["phase_angle"]) # minimum phase angle in data that was fit
                alpha_range = np.amax(data["phase_angle"]) - np.amin(data["phase_angle"])
                N_alpha_low=len(data["phase_angle"][data["phase_angle"]<self.low_alpha_cut]) # Number of data points at low phase angle
                N_mag_err=len(data["merr"][data["merr"]<self.mag_err_threshold]) # Use leq? Number of data points with error below some threshold

                # apparition data properties
                df_obj["phase_curve_N_fit_{}".format(filt)]=N_data_fit
                df_obj["phase_curve_alpha_min_{}".format(filt)]=alpha_min
                df_obj["phase_angle_range_{}".format(filt)] = alpha_range
                df_obj["phase_curve_N_alpha_low_{}".format(filt)]=N_alpha_low
                df_obj["phase_curve_N_mag_err_{}".format(filt)]=N_mag_err
                df_obj["phase_curve_N_data_app_{}".format(filt)]=N_data_app

                # iterate over all models for this apparition
                for model_name,model_values in self.selected_models.items():

                    # cut on phase angle range only if using the Linear phase function
                    if model_name=="LinearPhaseFunc":
                        mask_alpha = ((data_filt["phase_angle"]>=self.phase_lin_min) & (data_filt["phase_angle"]<=self.phase_lin_max))
                        N_alpha_lin=len(data_filt[~mask_alpha])
                        data_filt=data_filt[mask_alpha]
                        print("CUT N_alpha_lin = {}".format(N_alpha_lin))
                    else:
                        N_alpha_lin=0

                    print(filt,model_name,model_values)

                    # retrieve model names etc
                    ms=model_values["model_name_short"]
                    pc=model_values["model_parameters"]
                    model=model_values["model_function"]

                    # check there is sufficient data points to try fit
                    if N_data_fit<=len(pc):
                        print("less data to fit than parameters")
                        continue

                    if (epoch_ind==epoch_prime_ind) and filt=="o":
                        # set up model for fitting the primary o filter apparition
                        fit_slope = True
                        model.H = H_abs_mag # use the initial guess for H
                        # initial guess for slope parameters
                        for x in pc[1:]:
                            # choose G1 and G2 to be similar to usual default of G = 0.15
                            if x=="G1":
                                _G_slope = 0.5
                            elif x=="G2":
                                _G_slope = 0.2
                            else:
                                _G_slope = G_slope
                            setattr(model, "{}".format(x), _G_slope)

                    else:
                        # set up the model with fixed slope for every other apparition
                        fit_slope = False

                    if fit_slope: # store the first HG fit used for cut to plotting later if needed
                        _model_HG = model_HG
                        _data_residual = data_residual

                    # record whether slope is fitted. skip if filter == c as o is done first
                    if filt=="o":
                        df_obj["fit_slope"]=fit_slope
                        print("fit_slope={}".format(df_obj["fit_slope"]))

                    # set whether slope params are fixed
                    for x in pc[1:]:
                        getattr(model, '{}'.format(x)).fixed = not fit_slope
                        print("slope fixed = {}".format(getattr(model, '{}'.format(x)).fixed))

                    # fit the model to the data
                    # !!! ADD ANY NEW FITTING CODE HERE

                    model_fit = self.fitter(model, alpha, mag, weights=1.0/np.array(mag_err))
                    params = model_fit.parameters

                    # set the new H_abs_mag and slope parameters for subsequent fits
                    for j,x in enumerate(pc):
                        if x=="H":
                            setattr(model, "{}".format(x), params[j]*u.mag)
                        else:
                            setattr(model, "{}".format(x), params[j])

                    # store the fit and apparition parameters

                    # retrieve the fit metrics
                    x_vals=params
                    param_cov=self.fitter.fit_info['param_cov'] # see notes of https://docs.astropy.org/en/stable/api/astropy.modeling.fitting.LevMarLSQFitter.html for difference between param_cov and cov_x
                    if param_cov is None:
                        print("A value of None indicates a singular matrix, which means the curvature in parameters x is numerically flat")
                        # What should I do here?
                        break

                    # retrieve errors in the parameters: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
                    param_err_x = np.sqrt(np.diag(param_cov))
                    param_shape=param_cov.shape
                    correl_coef=np.zeros(param_shape)
                    for l in range(param_shape[0]):
                        for m in range(param_shape[1]):
                            correl_coef[l,m]=param_cov[l,m]/np.sqrt(param_cov[l,l]*param_cov[m,m]) # Hughes AND Hase eqn. 7.3

                    # calculate the residual properties
                    residuals = mag - model(alpha)
                    OC_mean = np.mean(residuals)
                    OC_std = np.std(residuals)
                    OC_range = np.absolute(np.amax(residuals)-np.amin(residuals))

                    # ALSO check for any nans in metrics? e.g. OC?

                    print(self.fitter.fit_info['message'])

                    print(df_obj["detection_count_{}".format(filt)])

                    # populate the df_obj dataframe: add all fit params/metrics to df_obj

                    # phase curve model fit properties
                    df_obj["phase_curve_OC_mean{}_{}".format(ms,filt)]=OC_mean
                    df_obj["phase_curve_OC_std{}_{}".format(ms,filt)]=OC_std
                    df_obj["phase_curve_OC_range{}_{}".format(ms,filt)]=OC_range
                    # if slope is fixed there is no error on it, need to add nan to error list
                    if not fit_slope:
                        extra_err = np.array([np.nan]*len(pc[1:]))
                        param_err_x = np.concatenate((param_err_x,extra_err))
                    # print(params)
                    # print(param_err_x)

                    for p in range(len(pc)):
                        df_obj["phase_curve_{}{}_{}".format(pc[p],ms,filt)]=params[p]
                        df_obj["phase_curve_{}_err{}_{}".format(pc[p],ms,filt)]=param_err_x[p]

                    # # make a plot?
                    # if self.plot_fig:
                    #
                    #     if not self.show_fig:
                    #         import matplotlib
                    #         print("use agg")
                    #         matplotlib.use('agg') # use agg backend to stop python stealing focus when plotting
                    #
                    #     import matplotlib.pyplot as plt
                    #     import matplotlib.gridspec as gridspec
                    #
                    #     fig = plt.figure()
                    #     gs = gridspec.GridSpec(1,1)
                    #     ax1 = plt.subplot(gs[0,0])
                    #
                    #     # plot all the apparition data
                    #     ax1.errorbar(data['phase_angle'],data['reduced_mag'],data['merr'], fmt='ko',label="fitted data",zorder=0,markersize="2")
                    #
                    #     # highlight rejected data
                    #     alpha_reject = np.concatenate([np.array(data_zero_err['phase_angle']),
                    #     np.array(data_small_err['phase_angle']),np.array(data_gal['phase_angle']),
                    #     np.array(data_orbfit['phase_angle']),np.array(data_lim_mag["phase_angle"])])
                    #     mag_reject = np.concatenate([np.array(data_zero_err['reduced_mag']),
                    #     np.array(data_small_err['reduced_mag']),np.array(data_gal['reduced_mag']),
                    #     np.array(data_orbfit['reduced_mag']),np.array(data_lim_mag["reduced_mag"])])
                    #
                    #     # plot phase curve fit
                    #     alpha_fit=np.linspace(0,np.amax(alpha),250)
                    #     ax1.plot(alpha_fit,model_fit(alpha_fit))
                    #
                    #     # set y lims to better show the phase curve, not rejected data
                    #     yshift = 0.5
                    #     plt.ylim(np.amin(data['reduced_mag'])-yshift, np.amax(data['reduced_mag'])+yshift)
                    #
                    #     ax1.set_xlabel('phase angle (degrees)')
                    #     ax1.set_ylabel('reduced magnitude')
                    #     ax1.invert_yaxis()
                    #     ax1.legend()
                    #
                    #     ax1.set_title("{}_{}_{}_{}_{}_{}".format(os.path.basename(__file__).split('.')[0],self.file_identifier,epoch_ind,model_name,self.clip_label,filt))
                    #     plt.tight_layout()
                    #
                    #     if self.save_fig:
                    #         fname="{}/{}_{}_{}_{}_{}_{}_fancy{}.{}".format(self.save_path,os.path.basename(__file__).split('.')[0],self.file_identifier,epoch_ind,model_name,self.clip_label,filt,self.save_file_suffix,self.save_file_type)
                    #         print(fname)
                    #         plt.savefig(fname, bbox_inches='tight')
                    #
                    #     if self.show_fig:
                    #         plt.show()
                    #     else:
                    #         plt.close()
                    #
                    #     plt.style.use('default')

            # push the parameters for this epoch to database
            if self.push_fit==True:
                # push_obj_db SHOULD only push selected columns in db_columns
                # ADD LinearPhaseFunc fields to db!
                self.push_obj_db(df_obj)

            df_obj_all = pd.concat([df_obj_all,df_obj])
            # print(df_obj_all.iloc[0].to_string())
            # print(len(list(df_obj_all)))

        # close database connections
        self.cnx1.disconnect()
        self.cnx2.disconnect()

        # save the data that was fit?
        if self.save_data_flag:
            save_file = "{}/df_data_fit_{}.csv".format(self.save_path,self.file_identifier)
            print("save data: {}".format(save_file))
            data_all_filt.to_csv(save_file)

        # make a plot?
        if self.plot_fig:

            if not self.show_fig:
                import matplotlib
                print("use agg")
                matplotlib.use('agg') # use agg backend to stop python stealing focus when plotting
            else:
                import matplotlib
                print("use TkAgg")
                matplotlib.use("TkAgg")

            import matplotlib.pyplot as plt
            import matplotlib.gridspec as gridspec

            fig = plt.figure()
            gs = gridspec.GridSpec(1,2)
            ax1 = plt.subplot(gs[0,0])
            ax2 = plt.subplot(gs[0,1], sharey = ax1)

            # define the phase angle range
            alpha_fit=np.linspace(0,np.amax(data_all_filt["phase_angle"]),250) * u.deg

            # # set y lims to better show the phase curve, not rejected data
            # yshift = 0.5
            # ax1.set_ylim(np.amin(data_all_filt['reduced_mag'])-yshift, np.amax(data_all_filt['reduced_mag'])+yshift)
            # ax2.set_ylim(np.amin(data_all_filt['reduced_mag'])-yshift, np.amax(data_all_filt['reduced_mag'])+yshift)

            # iterate over all models for this apparition
            for model_name,model_values in self.selected_models.items():

                for ax,filt in zip([ax1,ax2],["o","c"],):

                    ax.set_title("{}_{}_{}".format(self.file_identifier,model_name,filt))

                    for k in range(len(df_epoch)):

                        epoch_ind = int(df_epoch.iloc[len(df_epoch)-1-k]["epoch"])

                        # retrieve model names etc
                        ms=model_values["model_name_short"]
                        pc=model_values["model_parameters"]
                        model=model_values["model_function"]

                        # get obs data
                        mask_date = (((data_all_filt["mjd"]>=epochs[epoch_ind]) & (data_all_filt["mjd"]<epochs[epoch_ind+1])))
                        mask_filt =  (data_all_filt["filter"]==filt)
                        data = data_all_filt[mask_date & mask_filt]

                        # data_cut = data_all_cut[((data_all_cut["mjd"]>=epochs[epoch_ind]) & (data_all_cut["mjd"]<epochs[epoch_ind+1])) &
                        #                         (data_all_cut["filter"]==filt)]
                        # ax.scatter(data_cut['phase_angle'],data_cut['reduced_mag'],zorder=0,c="r",s=2)

                        # plot all the apparition data
                        ax.errorbar(data['phase_angle'],data['reduced_mag'],data['merr'], fmt='k.',zorder=0,markersize="2",alpha=0.3)
                        ax.scatter(data['phase_angle'],data['reduced_mag'],zorder=1,c="C{}".format(k),s=2)

                        # get the fit data
                        _df_obj = df_obj_all[df_obj_all["app_ind"]==epoch_ind]

                        for p in range(len(pc)):
                            x = _df_obj.iloc[0]["phase_curve_{}{}_{}".format(pc[p],ms,filt)]
                            if pc[p]=="H":
                                x *= u.mag
                            setattr(model, "{}".format(pc[p]), x)

                        # no fit params
                        if np.isnan(model.H):
                            continue

                        # plot phase curve fit
                        fit_label = "{},".format(epoch_ind)
                        for p in range(len(pc)):
                            x=getattr(model, "{}".format(pc[p]))
                            fit_label += " {}={:.2f}".format(pc[p],x.value)

                        ax.plot(alpha_fit,model(alpha_fit), c = "C{}".format(k), label = fit_label)

                    # if filt=="o": # plot the first fit used to clip data
                    #     ax.plot(alpha_fit,_model_HG(alpha_fit), c = "k", label = "initial cut")
                    #     # ax.scatter(_data_residual["phase_angle"],_data_residual["reduced_mag"],edgecolor="r",facecolor="none")

                break

            # # set y lims to better show the phase curve, not rejected data
            # yshift = 0.5
            # plt.ylim(np.median(data['reduced_mag'])-yshift, np.median(data['reduced_mag'])+yshift)

            ax1.set_xlabel('phase angle (degrees)')
            ax2.set_xlabel('phase angle (degrees)')
            ax1.set_ylabel('reduced magnitude')
            ax1.legend()
            ax2.legend()
            ax1.invert_yaxis()

            plt.tight_layout()

            if self.save_fig:
                fname="{}/{}_{}_{}{}.{}".format(self.save_path,os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.save_file_suffix,self.save_file_type)
                print(fname)
                plt.savefig(fname, bbox_inches='tight')

            if self.show_fig:
                plt.show()
            else:
                plt.close()

            plt.style.use('default')

        # !!! set dtypes of dataframe before returning?

        return df_obj_all # return the fit df of all epochs
