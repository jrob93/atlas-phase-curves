""" Module to perform ATLAS asteroid phase fits with various models and parameters """

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

# # importing our classes depends on whether this file as a script or module
# if __name__ == "__main__": # import other modules as a script
#     from atlas_SQL_query_df import atlas_SQL_query
#     from atlas_database_connection import database_connection
# else: # import as a package
#     from calculate_phase.atlas_SQL_query_df import atlas_SQL_query
#     from calculate_phase.atlas_database_connection import database_connection

from calculate_phase.atlas_SQL_query_df import atlas_SQL_query,get_orb_elements_id,get_unique_ids,atlas_SQL_query_orbid
from calculate_phase.atlas_database_connection import database_connection

class phase_fit():

    # default values
    low_alpha_cut=5.0*u.deg # we want to quantify how many data points are fit at low phase angles, alpha < low_alpha_cut
    param_converge_check=0.01 # the model is fit until the change in parameters (e.g. H and G) is less than param_converge_check (or max iterations is reached)
    max_iters=30 # maximum number of attempts at fitting and cutting
    std=2 # standard deviation of the sigma data clip
    mag_err_threshold = 0.1 # limit for the error of "good" data, we record N_mag_err number of data points with error < mag_err_threshold
    mag_err_small = 0.01 # we discount observations with error less than this
    # mag_err_small = 0.005 # we discount observations with error less than this
    gal_lat_cut=10 # galatic latitude cut in degrees
    mag_med_cut=2 # initial magnitude difference cut on initial HG model

    # set the clipping method for the phase_fit class
    def data_clip_sigma(self,data,data_predict,low=std,high=std):
        """ function to cut outliers by standard deviation, i.e. sigma clipping"""
        std=np.std(data)
        clip_mask=((data < (data_predict-(std*low))) | (data > (data_predict+(std*high))))
        return clip_mask

    def data_clip_diff(self,data,data_predict,diff=1):
        # cut outliers by diff (this function doesn't like astropy units, use np arrays)
        x=np.array(np.absolute(data_predict-data))
        clip_mask=(x>diff)
        return clip_mask

    # define which clipping function to use
    data_clip_func=data_clip_sigma
    clip_label="{}-sigma_clip".format(std)

    # set up the sbpy fitter and models
    # https://sbpy.readthedocs.io/en/latest/sbpy/photometry.html#disk-integrated-phase-function-models
    fitter = LevMarLSQFitter()

    # add functionality to choose models!!!
    all_models={
    "HG":{"model_name_short":"_B89","model_parameters":["H","G"],"model_function":HG()},
    "HG1G2":{"model_name_short":"_3M10","model_parameters":["H","G1","G2"],"model_function":HG1G2()},
    "HG12":{"model_name_short":"_2M10","model_parameters":["H","G12"],"model_function":HG12()},
    "HG12_Pen16":{"model_name_short":"_P16","model_parameters":["H","G12"],"model_function":HG12_Pen16()},
    "LinearPhaseFunc":{"model_name_short":"_Lin","model_parameters":["H","S"],"model_function":LinearPhaseFunc(H=15,S=0.04)}
    }

    # # use only the phase functions
    # model_names_str = ["HG", "HG1G2", "HG12", "HG12_Pen16"]
    # model_short = ["_B89","_3M10","_2M10","_P16"] # use shorthand for all models
    # phase_curve = [["H","G"],["H","G1","G2"],["H","G12"],["H","G12"]]
    # model_names = [HG(), HG1G2(), HG12(), HG12_Pen16()]

    # model_names_str = ["HG"]
    # model_short = ["_B89"] # use shorthand for all models
    # phase_curve = [["H","G"]]
    # model_names = [HG()]

    # # only linear fit
    # model_names_str = ["LinearPhaseFunc"]
    # model_short = ["_Lin"] # use shorthand for all models
    # phase_curve = [["H","S"]]
    # model_names = [LinearPhaseFunc(H=15,S=0.04)]

    # # do phase functions and linear fit
    # model_names_str = ["HG", "HG1G2", "HG12", "HG12_Pen16","LinearPhaseFunc"]
    # model_short = ["_B89","_3M10","_2M10","_P16","_Lin"] # use shorthand for all models
    # phase_curve = [["H","G"],["H","G1","G2"],["H","G12"],["H","G12"],["H","S"]]
    # model_names = [HG(), HG1G2(), HG12(), HG12_Pen16(),LinearPhaseFunc(H=15,S=0.04)] # NB that the sbpy LinearPhaseFunc class requires initial values, we use the ~mean H_B89_o and a sensible S, the other classes have built in default values

    phase_lin_min=5 # minimum phase angle for linear fit
    phase_lin_max=25 # maximum phase angle for linear fit

    # load the columns used to make the db table (might need to update the path)
    # fname_path = Path(__file__).parent / "../create_table/atlas_objects_fields.txt"
    fname_path = "{}/{}".format(Path(__file__).parent,"../create_table/atlas_objects_fields.txt")
    # print(fname_path)
    with open(fname_path,"r") as f:
        db_columns=f.readlines()
    db_columns=[d.rstrip() for d in db_columns]
    # print(db_columns)

    def __init__(self,
        mpc_number=False,
        name=False,
        save_path=".",
        save_file_suffix="",
        save_file_type="png",
        push_fit_flag=False,plot_fig_flag=False,show_fig_flag=False,save_fig_flag=False,hide_warning_flag=False,
        start_date=False,end_date=False,
        mag_diff_flag=False,
        H_abs_mag_o=False,H_abs_mag_c=False,
        model_list=["HG", "HG1G2", "HG12", "HG12_Pen16"], # ADD LinearPhaseFunc here as default?
        filter_list=["o","c"],
        tab_name="atlas_phase_fits"):

        # set up the class
        # define object and some flags

        # the object will have EITHER both mpc number and name OR just a name
        self.mpc_number=mpc_number
        self.name=name

        # set the variable for naming output files - N.B. what if name is passed but object does have mpc_number?
        if name:
            self.file_identifier="_".join((self.name).split())
        else:
            self.file_identifier=self.mpc_number

        self.push_fit=push_fit_flag # flag to push fit to database
        self.plot_fig=plot_fig_flag # flag to generate plot for each object
        self.show_fig=show_fig_flag # flag to display interactive plot
        self.save_fig=save_fig_flag # flag to save the figure
        self.utc_date_now=datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S") # time at which the script is run (UTC)
        self.save_path=save_path # where to save figures
        self.save_file_suffix=save_file_suffix # option to add a suffix to the png of plot saved
        self.save_file_type=save_file_type # option to choose type of file saved, e.g. png or pdf
        self.mag_diff_flag=mag_diff_flag # flag to perform an additional mag diff cut on the initial HG model

        self.selected_models= {key: value for key, value in self.all_models.items() if key in model_list} # create the dictionary of models to use

        self.start_date=start_date # only use observations after this date (mjd)
        self.end_date=end_date # only use observations before this date (mjd)
        # self.date_hard_cap=? # set a fixed date beyond which we do not consider data - SET BY end_date
        self.filters=filter_list # list of filters to use: ["o"], ["c"], ["o","c"]

        # pass initial guesses for the H_abs_mag (overrides astorb value)
        self.H_abs_mag_o=H_abs_mag_o
        self.H_abs_mag_c=H_abs_mag_c

        # cnx1 is the connection where we retrieve the data from
        self.cnx1=database_connection().connect()
        self.cursor1 = self.cnx1.cursor()

        # cnx2 is the connection where we push the fit data to
        self.cnx2 = self.cnx1
        self.cursor2 = self.cursor1

        self.tab_name=tab_name # table name within the database

        # create an empty pandas series to hold all object info
        # self.df_obj_datafit=pd.Series(data=np.zeros(len(self.db_columns))+np.nan,index=self.db_columns) # create an empty series to store all values for the object: metadata and fit data
        self.df_obj_datafit=pd.DataFrame(data=[],columns=self.db_columns) # create an empty series to store all values for the object: metadata and fit data
        print(self.df_obj_datafit)

        # all objects in the db have a unique orbital_elements_id - retrieve this from db using the mpc_number or name
        self.orbital_elements_id=get_orb_elements_id(self.cnx1,self.mpc_number,self.name)

        # get both primaryId and orbital_elements_id
        # unique_ids=get_unique_ids(self.cnx1,self.mpc_number,self.name)
        # self.primaryId=unique_ids["primaryId"]
        # self.orbital_elements_id=unique_ids["orbital_elements_id"]

        if hide_warning_flag==1:
            print("hide warnings")
            warnings.filterwarnings('ignore')

        # set up a log file - This will append to an existing log file of the same name, delete the file if you need a new log
        logging.basicConfig(filename='{}.log'.format(os.path.basename(__file__).split('.')[0]), level=logging.INFO)
        # logging.info('start log file')
        # logging.warning('{} warning test!'.format(self.orbital_elements_id))

    # def get_obj_data(self,cnx,mpc_num,t_start=False,t_end=False):
    def get_obj_data(self,cnx,orbid,t_start=False,t_end=False):
        # load data to be fitted, loads both filters (o & c)
        # data_all_filt=atlas_SQL_query(cnx=cnx,mpc_number=mpc_num)
        data_all_filt=atlas_SQL_query_orbid(cnx,orbid)

        print("data before date cut = {}".format(len(data_all_filt)))
        # cut by date range. Note that different filters may have different dates/epochs!
        if t_start:
            data_all_filt=data_all_filt[data_all_filt['mjd']>float(t_start)]
        if t_end:
            print(data_all_filt[~(data_all_filt['mjd']<float(t_end))])
            data_all_filt=data_all_filt[data_all_filt['mjd']<float(t_end)]
            print(data_all_filt)
        # else: # if t_end is not defined default to the hard cap for date
        #     data_all_filt=data_all_filt[data_all_filt['mjd']<self.date_hard_cap]
        print("data after date cut = {}".format(len(data_all_filt)))

        # exit()

        # DROP ALL ROWS WITH A NAN? !!!
        print("data before nan cut = {}".format(len(data_all_filt)))
        data_all_filt=data_all_filt.dropna()
        print("data after nan cut = {}".format(len(data_all_filt)))

        return data_all_filt

    # def get_obj_metadata(self,cnx,mpc_num):
    def get_obj_metadata(self,cnx,orbid):

        # get the last bits of data. might be out of date if rockatlas isn't running...?
        # qry_obj1=u"""SELECT
        # dateLastModified,
        # detection_count,
        # detection_count_c,
        # detection_count_o,
        # last_detection_mjd,
        # last_photometry_update_date_c,
        # last_photometry_update_date_o,
        # mpc_number,
        # name,
        # orbital_elements_id,
        # primaryId,
        # updated
        # FROM atlas_objects WHERE mpc_number=%(mpc_num)s
        # """ % locals()
        qry_obj1=u"""SELECT
        dateLastModified,
        detection_count,
        detection_count_c,
        detection_count_o,
        last_detection_mjd,
        last_photometry_update_date_c,
        last_photometry_update_date_o,
        mpc_number,
        name,
        orbital_elements_id,
        primaryId,
        updated
        FROM atlas_objects WHERE orbital_elements_id=%(orbid)s
        """ % locals()
        # print(qry_obj1)

        df_obj=pd.read_sql_query(qry_obj1,cnx)
        # print(pd.Series(df_obj.iloc[0]))
        # print(df_obj.to_string())
        # print(df_obj[['detection_count',  'detection_count_c',  'detection_count_o']])

        # fix the date strings
        dates=["dateLastModified","last_photometry_update_date_c","last_photometry_update_date_o"]
        for d in dates:
            if df_obj[d].iloc[0] is not None:
                df_obj.loc[0,d]="'{}'".format(str(df_obj[d].iloc[0]))
        # for some reason None needs changed to Null?
        df_obj=df_obj.fillna(value="NULL")
        # print(df_obj.to_string())
        # print(df_obj[['detection_count',  'detection_count_c',  'detection_count_o']])
        return df_obj

    # def get_obj_HG(self,cnx,mpc_num):
    def get_obj_HG(self,cnx,orbid):
        # get the astorb H and G (B89)
        # qry_HG="select G_slope, H_abs_mag from orbital_elements where primaryId='{}';".format(self.mpc_number)
        # qry_HG="SELECT G_slope, H_abs_mag FROM orbital_elements WHERE mpc_number='{}';".format(mpc_num)
        qry_HG="SELECT G_slope, H_abs_mag FROM orbital_elements WHERE orbital_elements_id='{}';".format(orbid)
        # print(qry_HG)
        df_HG=pd.read_sql_query(qry_HG,cnx)
        return df_HG

    # def get_obj_meta(self,cnx,mpc_num):
    def get_obj_meta(self,cnx,orbid,mpc_num,name):

        # get the last bits of data. might be out of date if rockatlas isn't running...?
        # qry_obj1=u"""SELECT
        # a.dateLastModified,
        # a.detection_count,
        # a.detection_count_c,
        # a.detection_count_o,
        # a.last_detection_mjd,
        # a.last_photometry_update_date_c,
        # a.last_photometry_update_date_o,
        # a.mpc_number,
        # a.name,
        # a.orbital_elements_id,
        # a.primaryId,
        # a.updated,
        # o.G_slope,
        # o.H_abs_mag
        # FROM atlas_objects a, orbital_elements o WHERE a.mpc_number=%(mpc_num)s AND o.mpc_number=%(mpc_num)s
        # """ % locals()

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
            o.H_abs_mag
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
            o.H_abs_mag
            FROM atlas_objects a, orbital_elements o WHERE a.orbital_elements_id=%(orbid)s AND o.name="%(name)s";
            """ % locals()

        # qry_obj1=u"""SELECT
        # atlas_objects.dateLastModified,
        # atlas_objects.detection_count,
        # atlas_objects.detection_count_c,
        # atlas_objects.detection_count_o,
        # atlas_objects.last_detection_mjd,
        # atlas_objects.last_photometry_update_date_c,
        # atlas_objects.last_photometry_update_date_o,
        # atlas_objects.mpc_number,
        # atlas_objects.name,
        # atlas_objects.orbital_elements_id,
        # atlas_objects.primaryId,
        # atlas_objects.updated,
        # orbital_elements.G_slope,
        # orbital_elements.H_abs_mag
        # FROM atlas_objects
        # INNER JOIN orbital_elements ON atlas_objects.primaryId=orbital_elements.primaryId
        # WHERE atlas_objects.orbital_elements_id=%(orbid)s;
        # """ % locals()

        print(qry_obj1)

        df_obj=pd.read_sql_query(qry_obj1,cnx)
        print(df_obj[["primaryId","orbital_elements_id","mpc_number","name"]])
        # print(list(df_obj))

        # fix the date strings
        dates=["dateLastModified","last_photometry_update_date_c","last_photometry_update_date_o"]
        for d in dates:
            if df_obj[d].iloc[0] is not None:
                df_obj.loc[0,d]="'{}'".format(str(df_obj[d].iloc[0]))
        # for some reason None needs changed to Null?
        df_obj=df_obj.fillna(value="NULL")
        # print(df_obj.to_string())
        # print(df_obj[['detection_count',  'detection_count_c',  'detection_count_o']])
        return df_obj

    # def database_obj_entry(self,mpc_num,df_obj,mpc_check):
    def database_obj_entry(self,orbid,df_obj,mpc_check):

        if mpc_check==0:

            qry_obj = u"""INSERT INTO {}
            (mpc_number,
            dateLastModified,
            detection_count,
            last_detection_mjd,
            last_photometry_update_date_c,
            last_photometry_update_date_o,
            name,
            orbital_elements_id,
            primaryId)
            VALUES
            (%(mpc_num)s,{},{},{},{},{},"{}",{},{});""".format(self.tab_name,
            str(df_obj['dateLastModified'].iloc[0]),
            int(df_obj['detection_count'].iloc[0]),
            float(df_obj['last_detection_mjd']),
            str(df_obj['last_photometry_update_date_c'].iloc[0]),
            str(df_obj['last_photometry_update_date_o'].iloc[0]),
            str(df_obj['name'].iloc[0]),
            int(df_obj['orbital_elements_id']),
            int(df_obj['primaryId'])) % locals()

        else:
            qry_obj = u"""UPDATE {} SET
            dateLastModified={},
            detection_count={},
            last_detection_mjd={},
            last_photometry_update_date_c={},
            last_photometry_update_date_o={},
            name="{}",
            orbital_elements_id={},
            primaryId={}
            WHERE mpc_number=%(mpc_num)s;""".format(self.tab_name,
            str(df_obj['dateLastModified'].iloc[0]),
            int(df_obj['detection_count'].iloc[0]),
            float(df_obj['last_detection_mjd']),
            str(df_obj['last_photometry_update_date_c'].iloc[0]),
            str(df_obj['last_photometry_update_date_o'].iloc[0]),
            str(df_obj['name'].iloc[0]),
            int(df_obj['orbital_elements_id']),
            int(df_obj['primaryId'])) % locals()

        return qry_obj

    def database_obj_entry2(self,mpc_num,df_obj,mpc_check):

        insert_obj = u"""(mpc_number,
        dateLastModified,
        detection_count,
        last_detection_mjd,
        last_photometry_update_date_c,
        last_photometry_update_date_o,
        name,
        orbital_elements_id,
        primaryId)
        VALUES
        (%(mpc_num)s,{},{},{},{},{},"{}",{},{})\n""".format(str(df_obj['dateLastModified'].iloc[0]),
        int(df_obj['detection_count'].iloc[0]),
        float(df_obj['last_detection_mjd']),
        str(df_obj['last_photometry_update_date_c'].iloc[0]),
        str(df_obj['last_photometry_update_date_o'].iloc[0]),
        str(df_obj['name'].iloc[0]),
        int(df_obj['orbital_elements_id']),
        int(df_obj['primaryId'])) % locals()
        # insert_obj="insert"
        update_obj = u"""dateLastModified={},
        detection_count={},
        last_detection_mjd={},
        last_photometry_update_date_c={},
        last_photometry_update_date_o={},
        name="{}",
        orbital_elements_id={},
        primaryId={}""".format(str(df_obj['dateLastModified'].iloc[0]),
        int(df_obj['detection_count'].iloc[0]),
        float(df_obj['last_detection_mjd']),
        str(df_obj['last_photometry_update_date_c'].iloc[0]),
        str(df_obj['last_photometry_update_date_o'].iloc[0]),
        str(df_obj['name'].iloc[0]),
        int(df_obj['orbital_elements_id']),
        int(df_obj['primaryId'])) % locals()
        # update_obj="update"
        qry_obj=u"""INSERT INTO {} {} ON DUPLICATE KEY UPDATE {};""".format(self.tab_name,insert_obj,update_obj) % locals()

        return qry_obj

    def push_fit_db(self,params,param_err_x,pc,ms,filt,detection_count_filt,N_data_fit,alpha_min,alpha_max,phase_angle_range,N_alpha_low,N_nights,N_iter,nfev,ier,N_mag_err):

        # how to pass large numbers of parameters? Define a dataframe in the class?
        print("push fit to db")

        HG_params_str=""
        for p in range(len(pc)):
            HG_params_str+="phase_curve_{}%(ms)s_%(filt)s={},phase_curve_{}_err%(ms)s_%(filt)s={},".format(
            pc[p],params[p],
            pc[p],param_err_x[p]
            ) % locals()

        tab_name=self.tab_name
        utc_date_now=self.utc_date_now
        mpc_number=self.mpc_number

        qry = u"""UPDATE %(tab_name)s SET
        detection_count_%(filt)s=%(detection_count_filt)s,
        %(HG_params_str)s
        phase_curve_N_fit%(ms)s_%(filt)s=%(N_data_fit)s,
        phase_curve_alpha_min%(ms)s_%(filt)s=%(alpha_min)s,
        phase_curve_alpha_max%(ms)s_%(filt)s=%(alpha_max)s,
        phase_angle_range_%(filt)s=%(phase_angle_range)s,
        phase_curve_N_alpha_low%(ms)s_%(filt)s=%(N_alpha_low)s,
        phase_curve_N_nights%(ms)s_%(filt)s=%(N_nights)s,
        phase_curve_N_iter%(ms)s_%(filt)s=%(N_iter)s,
        phase_curve_refresh_date_%(filt)s='%(utc_date_now)s',
        phase_curve_nfev%(ms)s_%(filt)s=%(nfev)s,
        phase_curve_ier%(ms)s_%(filt)s=%(ier)s,
        phase_curve_N_mag_err%(ms)s_%(filt)s=%(N_mag_err)s
        WHERE mpc_number=%(mpc_number)s;""" % locals()

        qry=qry.replace('=nan', '=NULL') # mysql doesn't understand nan, needs NULL

        # print(qry)
        # exit()

        self.cursor2.execute(qry)
        self.cnx2.commit()

    def push_obj_db(self,df_obj):
        """ This function will take the dataframe from the fit and only push the columns in self.db_columns
        These are the columns in the db, other columns will break the sql query.
        G_slope, H_abs_mag and the linear phase fits will not be pushed
        """

        # how to pass large numbers of parameters? Define a dataframe in the class?
        print("push fit to db")

        tab_name=self.tab_name
        utc_date_now=self.utc_date_now
        mpc_number=self.mpc_number

        cols=self.db_columns
        vals=df_obj[cols].iloc[0]
        vals=np.array(vals).astype(str)

        for c,v in zip(cols,vals):
            print(c,v)

        # print(cols)
        # print(list(df_obj))
        # print(vals)
        # print(len(cols),len(vals))
        # print(len(list(df_obj)))

        N_cols=len(cols)
        N_vals=len(vals)

        col_vals_update=""
        for i in range(N_cols):
            print("update {} {}".format(cols[i],vals[i]))
            # get rid of if statement if possible? better coding? list comprehension?

            # catch any values that are nan
            if vals[i].lower()=="nan": # mysql doesn't understand nan, needs NULL
                vals[i]="NULL"

            # catch any dates that will need fixed
            if cols[i] in ["name","phase_curve_refresh_date_o","phase_curve_refresh_date_c"]: # these fields are strings and need quotation marks
                col_vals_update+="{}=\"{}\",".format(cols[i],vals[i])
                vals[i]="\"{}\"".format(vals[i])
            else:
                col_vals_update+="{}={},".format(cols[i],vals[i])

        col_vals_update=col_vals_update[:-1] # drop the last comma
        # print(col_vals_update)

        # sanity check to make sure no values are dropped
        # print(col_vals_update.split(","))
        N_col_vals_update=len(col_vals_update.split(","))
        if N_col_vals_update!=N_vals or N_col_vals_update!=N_cols:
            print("check lengths: {} {} {}".format(len(col_vals_update.split(",")),len(vals),len(cols)))
            warning_message="{} - {} - SQL query lengths do not add up".format(self.mpc_number,self.name)
            print(warning_message)
            logging.warning(warning_message)

        # qry = u"""UPDATE %(tab_name)s SET {} WHERE mpc_number=%(mpc_number)s;""".format(col_vals_update) % locals()

        # UPDATE TO BE UPDATE/INSERT!
        qry=u"""INSERT INTO {} ({}) values ({}) ON DUPLICATE KEY UPDATE {};""".format(self.tab_name,",".join(cols), ",".join(vals), col_vals_update)

        # qry=qry.replace('=nan', '=NULL') # mysql doesn't understand nan, needs NULL

        print(qry)
        # exit()

        # print("\n!!!!!!!!!\nUNCOMMENT TO PUSH!\n!!!!!!!!!\n")
        self.cursor2.execute(qry)
        self.cnx2.commit()

        return

    def plot_phase_fit(self,model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
    data_filt,data_zero_err,data_small_err,data_gal):
        # plot a figure

        # how to pass large numbers of parameters?

        if not self.show_fig:
            import matplotlib
            print("use agg")
            matplotlib.use('agg') # use agg backend to stop python stealing focus when plotting

        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        # extract asteroid phase data from data, with units
        alpha = np.array(data['phase_angle']) * u.deg
        mag = np.array(data["reduced_mag"]) * u.mag
        mag_err = np.array(data["merr"]) * u.mag

        fig = plt.figure()
        # gs = gridspec.GridSpec(1,1)
        # ax1 = plt.subplot(gs[0,0])
        gs = gridspec.GridSpec(3,1,height_ratios=[1,1,0.1])
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[1,0])
        ax3 = plt.subplot(gs[2,0])

        # plot all the data from the SQL that goes into the fitting process
        ax1.errorbar(data_filt['phase_angle'],data_filt['reduced_mag'],data_filt['merr'], fmt='ko',label="data",zorder=0,markersize="2")

        # highlight any measurements with zero uncertainty
        ax1.scatter(data_zero_err['phase_angle'],data_zero_err['reduced_mag'],c='r',marker="+",s=50)
        ax1.scatter(data_small_err['phase_angle'],data_small_err['reduced_mag'],edgecolor='r',facecolor="none",marker="o",s=50)
        # highlight low galactic latitude
        ax1.scatter(data_gal['phase_angle'],data_gal['reduced_mag'],c='r',marker="x",s=50)

        # plot iterative fits and cuts
        alpha_fit=np.linspace(np.amin(alpha),np.amax(alpha),100)
        print(label_iter_list)
        print(model_iter_list)
        # print(k)
        for j in range(len(model_iter_list)):
            print(j,label_iter_list[j])
            ax1.plot(alpha_fit,model_iter_list[j](alpha_fit),label=label_iter_list[j])
            ax1.scatter(alpha_cut_iter_list[j],mag_cut_iter_list[j],marker="x",zorder=3)

        # plot the phase data after the cuts
        ax2.errorbar(np.array(alpha),np.array(mag),np.array(mag_err),fmt='k.',alpha=0.2,label="data error",zorder=0,markersize="2")
        s2=ax2.scatter(np.array(alpha),np.array(mag),c=np.array(data["mjd"]),label="data time",s=10)
        cbar2=fig.colorbar(s2,ax3,use_gridspec=True, orientation='horizontal')

        # plot the final fit
        ax2.plot(alpha_fit,model(alpha_fit),label=label)

        ax1.set_xlabel('alpha(degrees)')
        ax1.set_ylabel('mag')
        ax1.invert_yaxis()
        ax1.legend(prop={'size': 6})

        ax2.set_xlabel('alpha(degrees)')
        ax2.set_ylabel('mag')
        ax2.invert_yaxis()
        ax2.legend(prop={'size': 6})

        # fig.set_size_inches(4,10)

        ax3.set_xlabel("MJD")

        ax1.set_title("{}_{}_{}_{}_{}".format(os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt))
        plt.tight_layout()

        if self.save_fig:
            fname="{}/{}_{}_{}_{}_{}{}.{}".format(self.save_path,os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt,self.save_file_suffix,self.save_file_type)
            print(fname)
            plt.savefig(fname, bbox_inches='tight')

        if self.show_fig:
            plt.show()
        else:
            plt.close()

        return

    def plot_phase_fit_iteration(self,model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
    data_filt,data_zero_err,data_small_err,data_gal,data_diff):
        # plot a figure

        # how to pass large numbers of parameters?

        if not self.show_fig:
            import matplotlib
            print("use agg")
            matplotlib.use('agg') # use agg backend to stop python stealing focus when plotting

        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        # extract asteroid phase data from data, with units
        alpha = np.array(data['phase_angle']) * u.deg
        mag = np.array(data["reduced_mag"]) * u.mag
        mag_err = np.array(data["merr"]) * u.mag

        fig = plt.figure()
        gs = gridspec.GridSpec(1,1)
        ax1 = plt.subplot(gs[0,0])

        # plot all the data from the SQL that goes into the fitting process
        ax1.errorbar(data_filt['phase_angle'],data_filt['reduced_mag'],data_filt['merr'], fmt='ko',label="data",zorder=0,markersize="2")

        # highlight any measurements with zero uncertainty
        ax1.scatter(data_zero_err['phase_angle'],data_zero_err['reduced_mag'],edgecolor='r',facecolor="none",marker="^",s=50,label="{} mag_err=0".format(len(data_zero_err)))
        ax1.scatter(data_small_err['phase_angle'],data_small_err['reduced_mag'],edgecolor='r',facecolor="none",marker="s",s=50,label="{} mag_err<{}".format(len(data_small_err),self.mag_err_small))
        # highlight low galactic latitude
        ax1.scatter(data_gal['phase_angle'],data_gal['reduced_mag'],edgecolor='r',facecolor="none",marker="o",s=50,label="{} galactic_latitude<{}".format(len(data_gal),self.gal_lat_cut))

        if self.mag_diff_flag:
            #plot objects dropped in initial cut
            ax1.scatter(data_diff['phase_angle'],data_diff['reduced_mag'],edgecolor='r',facecolor="none",marker="p",s=50,label="{} HG model diff>{}".format(len(data_diff),self.mag_med_cut))

        # plot iterative fits and cuts
        alpha_fit=np.linspace(np.amin(alpha),np.amax(alpha),100)
        print(label_iter_list)
        print(model_iter_list)
        # print(k)
        for j in range(len(model_iter_list)):
            print(j,label_iter_list[j])
            ax1.plot(alpha_fit,model_iter_list[j](alpha_fit),label=label_iter_list[j])
            ax1.scatter(alpha_cut_iter_list[j],mag_cut_iter_list[j],marker="x",zorder=3)

        ax1.plot(alpha_fit,model(alpha_fit),label=label)

        ax1.set_xlabel('alpha(degrees)')
        ax1.set_ylabel('reduced mag')
        ax1.invert_yaxis()
        ax1.legend(prop={'size': 6})

        ax1.set_title("{}_{}_{}_{}_{}".format(os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt))
        plt.tight_layout()

        # ax1.set_ylim(14.5,8.3)

        if self.save_fig:
            # fname="{}/{}_{}_{}_{}_{}_iter{}.png".format(self.save_path,os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt,self.save_file_suffix)
            fname="{}/{}_{}_{}_{}_{}_iter{}.{}".format(self.save_path,os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt,self.save_file_suffix,self.save_file_type)
            print(fname)
            plt.savefig(fname, bbox_inches='tight')

        if self.show_fig:
            plt.show()
        else:
            plt.close()

        return

    def calculate(self):
        # calculate the phase curves

        # get the object observation data
        # data_all_filt=self.get_obj_data(self.cnx1,self.mpc_number,self.start_date,self.end_date)
        data_all_filt=self.get_obj_data(self.cnx1,self.orbital_elements_id,self.start_date,self.end_date)
        # print(data_all_filt)
        print(data_all_filt[np.isnan(data_all_filt["mjd"])])

        # get the object metadata and combine with the phase fit dataframe structure
        # df_obj=self.get_obj_metadata(self.cnx1,self.mpc_number)
        # df_obj=self.get_obj_meta(self.cnx1,self.mpc_number)
        df_obj=self.get_obj_meta(self.cnx1,self.orbital_elements_id,self.mpc_number,self.name)
        print(df_obj)
        d1=df_obj
        d2=self.df_obj_datafit
        # print(d1)
        # print(d2)
        df_obj=d2.append(d1) # add the values to the df
        # print(df_obj.iloc[0].to_string)
        # print(df_obj['phase_curve_H_2M10_o'])

        # update the detection_count
        df_obj["detection_count"]=len(data_all_filt) # note that atlas_objects may not have up to date detection count...

        # set the mpc_number and name from df_obj !!!
        self.mpc_number=df_obj.iloc[0]['mpc_number']
        self.name=df_obj.iloc[0]['name']

        # # get the astorb H and G (B89)
        # # qry_HG="select G_slope, H_abs_mag from orbital_elements where primaryId='{}';".format(self.mpc_number)
        # qry_HG="select G_slope, H_abs_mag from orbital_elements where mpc_number='{}';".format(self.mpc_number)
        # print(qry_HG)
        # df_HG=pd.read_sql_query(qry_HG,self.cnx1)
        # df_HG=self.get_obj_HG(self.cnx1,self.mpc_number)
        # print(df_HG)
        # exit()

        # # is the object already in the tab_name at the storage database?
        # self.cursor2.execute("SELECT COUNT(1) FROM {} where mpc_number={}".format(self.tab_name,self.mpc_number))
        # mpc_check=self.cursor2.fetchone()[0]
        # print(mpc_check)

        # DO THIS PUSH AT END
        # if self.push_fit==True: # only update the entry if we are going to push results
        #     # MOVE THIS TO PUSH FIT VALUES
        #     qry_obj = self.database_obj_entry(self.mpc_number,df_obj,mpc_check)
        #     # qry_obj = self.database_obj_entry2(self.mpc_number,df_obj,mpc_check)
        #     print(qry_obj)
        #     # self.cursor2.execute(qry_obj)
        #     # self.cnx2.commit()
        #     # exit()

        for filt in self.filters:

            # H and G values for the predicted fit
            # G_slope=float(df_HG.iloc[0]['G_slope'])
            # H_abs_mag=float(df_HG.iloc[0]['H_abs_mag'])
            G_slope=float(df_obj.iloc[0]['G_slope'])
            H_abs_mag=float(df_obj.iloc[0]['H_abs_mag'])
            print("G_slope = {}\nH_abs_mag = {}".format(G_slope,H_abs_mag))

            # do filter correction from V band (Heinze et al. 2020) - see also Erasmus et al 2020 for the c-o colours of S and C types (0.388 and 0.249 respectively)
            if filt=="o":
                if self.H_abs_mag_o:
                    print("override shifted astorb value of {}".format(H_abs_mag-0.332))
                    H_abs_mag=self.H_abs_mag_o
                    print(H_abs_mag)
                else:
                    H_abs_mag+=-0.332
            if filt=="c":
                if self.H_abs_mag_c:
                    print("override shifted astorb value of {}".format(H_abs_mag+0.054))
                    H_abs_mag=self.H_abs_mag_c
                    print(H_abs_mag)
                else:
                    H_abs_mag+=0.054

            # use our own guess for H, e.g. a previous fit?
            # if H_abs_mag_o is passed use this: # median of phase_curve_H_2M10_o,phase_curve_H_3M10_o,phase_curve_H_B89_o,phase_curve_H_P16_o
            # if H_abs_mag_c is passed use this: # median of phase_curve_H_2M10_c,phase_curve_H_3M10_c,phase_curve_H_B89_c,phase_curve_H_P16_c
            # else use the astorb guess:

            # use the updated field to track how many times this iteration has been done?

            # H_abs_mag=14.6833
            # H_vals=np.array([14.5801,np.nan,14.6833,14.4651])
            # H_abs_mag=21.9497
            # H_vals=np.array([np.nan,14.4716,21.9497,16.6153])
            # H_abs_mag=np.median(H_vals[~np.isnan(H_vals)])

            # select all data from a certain filter
            data_filt=data_all_filt[data_all_filt['filter']==filt]
            # print(data_filt)
            detection_count_filt=len(data_filt) # Record number of detections BEFORE cuts are made
            df_obj["detection_count_{}".format(filt)]=detection_count_filt # update the number of detections in the df

            # print(detection_count_filt)
            # continue

            # iterate over all models
            # for i,model_name in enumerate(self.model_names):
            for model_name,model_values in self.selected_models.items():

                print(model_name,model_values)

                # ms=self.model_short[i]
                # pc=self.phase_curve[i]
                ms=model_values["model_name_short"]
                pc=model_values["model_parameters"]
                model_func=model_values["model_function"]

                print(ms,pc,model_func)

                old_params=[999]*len(model_func.parameters)

                # store models/labels/cut data for plotting
                label_iter_list=[]
                model_iter_list=[]
                mag_cut_iter_list=[]
                alpha_cut_iter_list=[]

                # print("{}, {}: fit {}, filter {}".format(self.mpc_number,self.mpc_number,self.model_names_str[i],filt))
                print("{}, {}: fit {}, filter {}".format(self.name,self.mpc_number,model_name,filt))

                # initialise the data that we will iteratively fit and cut
                data=data_filt
                data=data.sort_values("phase_angle") # ensure that the dataframe is in order for plotting

                # print("DATA WITH NAN")
                # print(data[data.isna().any(axis=1)])

                print("{} starting data".format(len(data)))

                # # do an initial mag diff cut to remove extreme outliers
                # mag_med_cut=2
                # data=data[np.absolute(np.array(data["reduced_mag"]-np.nanmedian(data["reduced_mag"])))<mag_med_cut]
                # print("{} {} mag diff".format(len(data),mag_med_cut))

                # drop any measurements with zero uncertainty
                data_zero_err=data[data['merr']==0]
                data=data[data['merr']!=0]
                print("{} zero error".format(len(data)))
                # drop measurements with small (or zero) uncertainty
                data_small_err=data[data['merr']<self.mag_err_small]
                data=data[~(data['merr']<self.mag_err_small)]
                print("{} small error".format(len(data)))

                # drop measurements near galactic plane
                data_gal=data[np.absolute(data["galactic_latitude"])<self.gal_lat_cut]
                data=data[~(np.absolute(data["galactic_latitude"])<self.gal_lat_cut)]
                print("{} galactic plane".format(len(data)))

                # phase angle cut for linear fit
                if model_name=="LinearPhaseFunc":
                    data=data[(data["phase_angle"]>self.phase_lin_min) & (data["phase_angle"]<self.phase_lin_max)]

                print("{} data after cuts".format(len(data)))

                if len(data)==0:
                    print("no data, cannot fit")
                    break

                # iteratively fit and cut data
                k=0
                while k<self.max_iters:

                    # print("iteration: {}".format(k))
                    print("iteration: {}, N_data={}".format(k,len(data)))

                    # print("DATA WITH NAN")
                    # print(data[data.isna().any(axis=1)])

                    # extract asteroid phase data from data, with units
                    alpha = np.array(data['phase_angle']) * u.deg
                    mag = np.array(data["reduced_mag"]) * u.mag
                    mag_err = np.array(data["merr"]) * u.mag

                    if k==0:
                        # for first iteration start with the predicted HG mag
                        model=HG(H = H_abs_mag * u.mag, G = G_slope)
                        model_str="predicted mag HG"

                        param_names=model.param_names
                        params=model.parameters

                        # label the first fit
                        labels=[model_str]
                        for l in range(len(model.parameters)):
                            labels.append("{}={:.2f} ".format(param_names[l],params[l]))
                        label=", ".join(labels)

                        if self.mag_diff_flag:
                            # do a mag diff cut based on the initial assumed HG
                            mask=self.data_clip_diff(mag, model(alpha),self.mag_med_cut)
                            data_diff=data[mask]
                            data=data[~mask]
                            alpha = np.array(data['phase_angle']) * u.deg
                            mag = np.array(data["reduced_mag"]) * u.mag
                            mag_err = np.array(data["merr"]) * u.mag
                        else:
                            data_diff=[]


                    else:
                        # apply the fitter to the data
                        if len(alpha)<=len(pc): # check that there is enough data to fit
                            print("less data to fit than parameters")
                            break

                        # !!! RECORD ANY WARNINGS ETC? see fitter.fit_info
                        # logging will record these warnings

                        # model = self.fitter(model_name, alpha, mag, weights=1.0/np.array(mag_err)) # fit using weights by uncertainty
                        # if k==1:
                        #     model = fitter(model_name, alpha, mag) # drop the weights for the first fit
                        # else:
                        #     model = fitter(model_name, alpha, mag, weights=1.0/np.array(mag_err))

                        # warnings.filterwarnings("error")
                        # print("start fit")
                        # try:
                        #     model = self.fitter(model_name, alpha, mag, weights=1.0/np.array(mag_err)) # fit using weights by uncertainty
                        # except Warning:
                        #     warning_list=warnings.catch_warnings(*, record=False, module=None)
                        #     print ('Warning was raised as an exception!')
                        # except RuntimeWarning:
                        #     print ('RuntimeWarning was raised as an exception!')
                        # except AstropyWarning:
                        #     print ('AstropyWarning was raised as an exception!')
                        # except AstropyUserWarning:
                        #     print ('AstropyUserWarning was raised as an exception!')
                        # print("end fit")

                        with warnings.catch_warnings(record=True) as w:
                            model = self.fitter(model_func, alpha, mag, weights=1.0/np.array(mag_err)) # fit using weights by uncertainty
                            if len(w)>0:
                                warning_message="{} - {} - {} - {} - {}".format(self.mpc_number,self.name,model_name,filt,w[-1].message)
                                print(warning_message)
                                logging.warning(warning_message)

                        # model_name=self.model_names_str[i]

                        param_names=model.param_names
                        params=model.parameters
                        # print(model.__dict__)
                        # print()
                        # print(fitter.fit_info['cov_x']) # The scipy.optimize.leastsq result for the most recent fit https://docs.astropy.org/en/stable/api/astropy.modeling.fitting.LevMarLSQFitter.html
                        # print(err_x)

                        # label each fit
                        labels=["{}. {}".format(k,model_name)]
                        for l in range(len(model.parameters)):
                            labels.append("{}={:.2f} ".format(param_names[l],params[l]))
                        label=", ".join(labels)

                        # test for convergence
                        # find difference in params
                        delta_params = np.absolute(old_params-params)
                        # print(old_params)
                        # print(fitter.fit_info)
                        # print("nfev:{}".format(fitter.fit_info['nfev']))
                        # print("ierr:{}".format(fitter.fit_info['ierr']))
                        # print("message:{}".format(fitter.fit_info['message']))
                        # print("N_data:{}".format(len(data)))
                        # exit()
                        # print(params)
                        # print("delta params = {}".format(delta_params))

                        if np.sum(delta_params<self.param_converge_check)==len(delta_params):
                            # print("converged")
                            print(params)

                            # retrieve the fit metrics
                            x_vals=params
                            param_cov=self.fitter.fit_info['param_cov'] # see notes of https://docs.astropy.org/en/stable/api/astropy.modeling.fitting.LevMarLSQFitter.html for difference between param_cov and cov_x
                            if param_cov is None:
                                print("A value of None indicates a singular matrix, which means the curvature in parameters x is numerically flat")
                                # What should I do here?
                                break

                            # cov_x=self.fitter.fit_info['cov_x']
                            # print("cov_x:\n{}".format(cov_x))
                            # err_x = np.sqrt(np.diag(cov_x))
                            # print("err_x: {}".format(err_x))
                            # print(param_cov/cov_x)

                            param_err_x = np.sqrt(np.diag(param_cov)) # retrieve errors in the parameters: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
                            param_shape=param_cov.shape
                            correl_coef=np.zeros(param_shape)
                            for l in range(param_shape[0]):
                                for m in range(param_shape[1]):
                                    correl_coef[l,m]=param_cov[l,m]/np.sqrt(param_cov[l,l]*param_cov[m,m]) # Hughes AND Hase eqn. 7.3
                            N_data_fit=len(data)
                            # nfev=self.fitter.fit_info['nfev']
                            # print(np.array(data["mjd"]).astype(int))
                            N_nights=len(np.unique(np.array(data["mjd"]).astype(int)))
                            alpha_min=np.amin(alpha).value
                            alpha_max=np.amax(alpha).value
                            phase_angle_range=alpha_max-alpha_min
                            # N_alpha_low=sum(alpha<low_alpha_cut)
                            N_alpha_low=len(alpha[alpha<self.low_alpha_cut])
                            N_iter=k
                            nfev=self.fitter.fit_info['nfev']
                            ier=self.fitter.fit_info['ierr']
                            N_mag_err=len(mag_err[np.array(mag_err)<self.mag_err_threshold]) # leq?

                            # calculate the residual properties
                            residuals = mag - model(alpha)
                            # print(mag)
                            # print(residuals)
                            # if len(residuals)!=N_data_fit:
                            #     print("ERROR")
                            #     exit()
                            OC_mean = np.mean(residuals)
                            OC_std = np.std(residuals)
                            OC_range = np.absolute(np.amax(residuals)-np.amin(residuals))

                            print("std of residuals = {}".format(OC_std))
                            print("RMS of residuals = {}".format(np.sqrt(np.mean(residuals**2))))
                            # exit()

                            print("N_mag_err={}".format(N_mag_err))
                            if N_mag_err>N_data_fit:
                                # ERROR
                                logging.warning("{} - {} - {} - {} - N_mag_err>N_data_fit".format(self.mpc_number,self.name,model_name,filt))

                            # check for any nans?

                            # # time len() vs sum()
                            # import time
                            # start = time.process_time()
                            # sum(alpha<low_alpha_cut)
                            # sum(np.array(mag_err)<mag_err_threshold)
                            # print(time.process_time() - start)
                            #
                            # start = time.process_time()
                            # len(alpha[alpha<low_alpha_cut])
                            # len(mag_err[np.array(mag_err)<mag_err_threshold])
                            # print(time.process_time() - start)

                            # exit()

                            # print(self.fitter.fit_info)
                            print(self.fitter.fit_info['message'])

                            # if self.push_fit==True:
                            # STORE FIT DAT IN DF

                            # DO SOMETHING HERE TO EXTRACT THE LINEAR PHASE FIT?
                            if model_name=="LinearPhaseFunc":
                                print("record linear fit data")

                            print(df_obj["detection_count_{}".format(filt)])

                            # populate the df_obj dataframe
                            # add fit to df_obj
                            # %(HG_params_str)s
                            df_obj["phase_curve_N_fit{}_{}".format(ms,filt)]=N_data_fit
                            df_obj["phase_curve_alpha_min{}_{}".format(ms,filt)]=alpha_min
                            df_obj["phase_curve_alpha_max{}_{}".format(ms,filt)]=alpha_max
                            df_obj["phase_angle_range_{}".format(filt)]=phase_angle_range
                            df_obj["phase_curve_N_alpha_low{}_{}".format(ms,filt)]=N_alpha_low
                            df_obj["phase_curve_N_nights{}_{}".format(ms,filt)]=N_nights
                            df_obj["phase_curve_N_iter{}_{}".format(ms,filt)]=N_iter
                            df_obj["phase_curve_refresh_date_{}".format(filt)]=self.utc_date_now
                            df_obj["phase_curve_nfev{}_{}".format(ms,filt)]=nfev
                            df_obj["phase_curve_ier{}_{}".format(ms,filt)]=ier
                            df_obj["phase_curve_N_mag_err{}_{}".format(ms,filt)]=N_mag_err
                            df_obj["phase_curve_OC_mean{}_{}".format(ms,filt)]=OC_mean
                            df_obj["phase_curve_OC_std{}_{}".format(ms,filt)]=OC_std
                            df_obj["phase_curve_OC_range{}_{}".format(ms,filt)]=OC_range

                            for p in range(len(pc)):
                                df_obj["phase_curve_{}{}_{}".format(pc[p],ms,filt)]=params[p]
                                df_obj["phase_curve_{}_err{}_{}".format(pc[p],ms,filt)]=param_err_x[p]

                            # print(df_obj["detection_count_{}".format(filt)])
                            # print(df_obj.iloc[0].to_string())
                            # exit()

                            # self.push_fit_db(params,param_err_x,pc,ms,filt,detection_count_filt,N_data_fit,alpha_min,alpha_max,phase_angle_range,N_alpha_low,N_nights,N_iter,nfev,ier,N_mag_err)

                            # print("push fit to db")
                            #
                            # HG_params_str=""
                            # for p in range(len(self.phase_curve[i])):
                            #     HG_params_str+="phase_curve_{}%(ms)s_%(filt)s={},phase_curve_{}_err%(ms)s_%(filt)s={},".format(
                            #     self.phase_curve[i][p],params[p],
                            #     self.phase_curve[i][p],param_err_x[p]
                            #     ) % locals()
                            #
                            # tab_name=self.tab_name
                            # utc_date_now=self.utc_date_now
                            # mpc_number=self.mpc_number
                            #
                            # qry = u"""UPDATE %(tab_name)s SET
                            # detection_count_%(filt)s=%(detection_count_filt)s,
                            # %(HG_params_str)s
                            # phase_curve_N_fit%(ms)s_%(filt)s=%(N_data_fit)s,
                            # phase_curve_alpha_min%(ms)s_%(filt)s=%(alpha_min)s,
                            # phase_curve_alpha_max%(ms)s_%(filt)s=%(alpha_max)s,
                            # phase_angle_range_%(filt)s=%(phase_angle_range)s,
                            # phase_curve_N_alpha_low%(ms)s_%(filt)s=%(N_alpha_low)s,
                            # phase_curve_N_nights%(ms)s_%(filt)s=%(N_nights)s,
                            # phase_curve_N_iter%(ms)s_%(filt)s=%(N_iter)s,
                            # phase_curve_refresh_date_%(filt)s='%(utc_date_now)s',
                            # phase_curve_nfev%(ms)s_%(filt)s=%(nfev)s,
                            # phase_curve_ier%(ms)s_%(filt)s=%(ier)s,
                            # phase_curve_N_mag_err%(ms)s_%(filt)s=%(N_mag_err)s
                            # WHERE mpc_number=%(mpc_number)s;""" % locals()
                            #
                            # qry=qry.replace('=nan', '=NULL') # mysql doesn't understand nan, needs NULL
                            #
                            # print(qry)
                            # # exit()
                            #
                            # cursor2.execute(qry)
                            # cnx2.commit()

                            if self.plot_fig:

                                # self.plot_phase_fit(model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
                                # data_filt,data_zero_err,data_small_err,data_gal)
                                self.plot_phase_fit_iteration(model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
                                data_filt,data_zero_err,data_small_err,data_gal,data_diff)

                            # exit()

                            # # save data that was used to fit to file
                            # data_clip_file="results_analysis/fit_data/df_data_{}{}_{}.csv".format(self.mpc_number,ms,filt)
                            # print(data_clip_file)
                            # data.to_csv(data_clip_file)

                            break

                    # record stuff for plotting
                    label_iter_list.append(label)
                    model_iter_list.append(model)

                    # cut the outlying data
                    mask=self.data_clip_func(mag, model(alpha),self.std)

                    # record the cut data points for plotting
                    mag_cut_iter_list.append(mag[mask])
                    alpha_cut_iter_list.append(alpha[mask])

                    # define the data to keep
                    data=data[~mask]

                    # store these params as old params, to check for convergence of the next fit
                    if k>0:
                        old_params=params

                    k+=1

        print("N_mag_err before push = {}".format(df_obj.iloc[0]["phase_curve_N_mag_err_B89_o"]))
        print("N_fit before push = {}".format(df_obj.iloc[0]["phase_curve_N_fit_B89_o"]))

        # push ALL the data
        if self.push_fit==True:

            # # DO SOMETHING HERE TO STORE THE LINEAR PHASE FIT
            # if self.model_names_str[i]=="LinearPhaseFunc":
            #     print("store/push linear fit data")

            if df_obj.iloc[0]["phase_curve_N_mag_err_B89_o"]>df_obj.iloc[0]["phase_curve_N_fit_B89_o"]:
                # ERROR
                logging.warning("{} - {} - {} - {} - N_mag_err>N_data_fit".format(self.mpc_number,self.name,model_name,filt))

            self.push_obj_db(df_obj)

        # close all database connections
        # CLOSE CURSORS TOO?
        self.cnx1.disconnect()
        self.cnx2.disconnect()

        return df_obj # return the fit df

# if __name__ == "__main__":
#     # mpc_number=4986
#     mpc_number=1785
#     fit = phase_fit(mpc_number,push_fit_flag=False,plot_fig_flag=True,save_fig_flag=True,save_path="figs")#,show_fig_flag=True)
#     df=fit.calculate()
#     # print(df.iloc[0].to_string())
#     # df.to_csv("df_fit_{}.csv".format(mpc_number))
