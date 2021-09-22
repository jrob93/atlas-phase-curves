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
import calculate_phase.solar_apparitions as sa

class phase_fit():

    # default fixed values
    low_alpha_cut=5.0*u.deg # we want to quantify how many data points are fit at low phase angles, alpha < low_alpha_cut
    param_converge_check=0.01 # the model is fit until the change in parameters (e.g. H and G) is less than param_converge_check (or max iterations is reached)
    max_iters=30 # maximum number of attempts at fitting and cutting
    std=2 # standard deviation of the sigma data clip
    mag_err_threshold = 0.1 # limit for the error of "good" data, we record N_mag_err number of data points with error < mag_err_threshold
    # mag_err_small = 0.01 # we discount observations with error less than this
    mag_err_small = 0.005 # we discount observations with error less than this
    gal_lat_cut=10 # galatic latitude cut in degrees
    mag_med_cut=2 # initial magnitude difference cut on initial HG model
    phase_lin_min=5 # minimum phase angle for linear fit - MAKE OPTIONAL?
    phase_lin_max=25 # maximum phase angle for linear fit

    # set the clipping method for the phase_fit class
    def data_clip_sigma(self,data,data_predict,low=std,high=std):
        """ function to cut outliers by standard deviation, i.e. sigma clipping
        returns the mask of data points to cut"""
        std=np.std(data)
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
    data_clip_func=data_clip_sigma
    clip_label="{}-sigma_clip".format(std)

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
    fname_path = "{}/{}".format(Path(__file__).parent,"../create_table/atlas_objects_fields.txt")
    with open(fname_path,"r") as f:
        db_columns=f.readlines()
    db_columns=[d.rstrip() for d in db_columns]

    def __init__(self,
        mpc_number=False,
        name=False,
        save_path=".",
        save_file_suffix="",
        save_file_type="png",
        push_fit_flag=False,plot_fig_flag=False,show_fig_flag=False,save_fig_flag=False,hide_warning_flag=False,
        start_date=False,end_date=False,
        mag_diff_flag=False, # DEFAULT THIS TO BE TRUE?
        H_abs_mag_o=False,H_abs_mag_c=False,
        model_list=["HG", "HG1G2", "HG12", "HG12_Pen16"], # ADD LinearPhaseFunc here as default?
        filter_list=["o","c"],
        tab_name="atlas_phase_fits"):

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
        self.df_obj_datafit=pd.DataFrame(data=[],columns=self.db_columns) # create an empty series to store all values for the object: metadata and fit data
        print(self.df_obj_datafit)

        # all objects in the db have a unique orbital_elements_id - retrieve this from db using the mpc_number or name
        try:
            self.orbital_elements_id=get_orb_elements_id(self.cnx1,self.mpc_number,self.name)
        except:
            print("Cannot find object {}, {} in database".format(self.mpc_number,self.name))
            self.orbital_elements_id = None
            return

        # option to suppress warnings being shown to the screen
        if hide_warning_flag==1:
            print("hide warnings")
            warnings.filterwarnings('ignore')

        # set up a log file - This will append to an existing log file of the same name, delete the file if you need a new log
        logging.basicConfig(filename='{}.log'.format(os.path.basename(__file__).split('.')[0]), level=logging.INFO)

    def get_obj_data(self,cnx,orbid,t_start=False,t_end=False):
        """load data to be fitted, loads both filters (o & c)"""

        # query the database for photometry
        data_all_filt=atlas_SQL_query_orbid(cnx,orbid)

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

        # DROP ALL ROWS WITH A NAN
        N_data3 = len(data_all_filt)
        print("data before nan cut = {}".format(N_data3))
        data_all_filt=data_all_filt.dropna()
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

        for c,v in zip(cols,vals):
            print(c,v)

        N_cols=len(cols)
        N_vals=len(vals)

        # create a qry that sql can handle!
        col_vals_update=""
        for i in range(N_cols):
            print("update {} {}".format(cols[i],vals[i]))
            # get rid of if statement if possible? better coding? list comprehension?

            # catch any values that are nan
            if vals[i].lower()=="nan": # mysql doesn't understand nan, needs NULL
                vals[i]="NULL"

            # catch any names/dates that will need fixed
            if cols[i] in ["name","phase_curve_refresh_date_o","phase_curve_refresh_date_c"]: # these fields are strings and need quotation marks
                col_vals_update+="{}=\"{}\",".format(cols[i],vals[i])
                vals[i]="\"{}\"".format(vals[i])
            else:
                col_vals_update+="{}={},".format(cols[i],vals[i])

        col_vals_update=col_vals_update[:-1] # drop the last comma
        # print(col_vals_update)

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

    def plot_phase_fit(self,model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
    data_filt,data_zero_err,data_small_err,data_gal):
        """
        plot a figure.
        How best to pass large numbers of parameters?
        """

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
        gs = gridspec.GridSpec(2,1,height_ratios=[1,0.05])
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[1,0])

        # plot all the data from the SQL that goes into the fitting process
        ax1.errorbar(data_filt['phase_angle'],data_filt['reduced_mag'],data_filt['merr'], fmt='ko',label="data",zorder=0,markersize="2")
        s1=ax1.scatter(data_filt['phase_angle'],data_filt['reduced_mag'],c=data_filt['mjd'],s=10)
        cbar1=fig.colorbar(s1,ax2,use_gridspec=True, orientation='horizontal')
        cbar1.set_label("mjd")

        # highlight any measurements with zero uncertainty
        ax1.scatter(data_zero_err['phase_angle'],data_zero_err['reduced_mag'],edgecolor='r',facecolor="none",marker="^",s=50,label="{} mag_err = 0".format(len(data_zero_err)))
        ax1.scatter(data_small_err['phase_angle'],data_small_err['reduced_mag'],edgecolor='r',facecolor="none",marker="s",s=50,label="{} mag_err < {}".format(len(data_small_err),self.mag_err_small))
        # highlight low galactic latitude
        ax1.scatter(data_gal['phase_angle'],data_gal['reduced_mag'],edgecolor='r',facecolor="none",marker="o",s=50,label="{} galactic_latitude < {} deg".format(len(data_gal),self.gal_lat_cut))

        if self.mag_diff_flag:
            #plot objects dropped in initial cut
            ax1.scatter(data_diff['phase_angle'],data_diff['reduced_mag'],edgecolor='r',facecolor="none",marker="p",s=50,label="{} HG model diff > {}".format(len(data_diff),self.mag_med_cut))

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

        if self.save_fig:
            fname="{}/{}_{}_{}_{}_{}_iter{}.{}".format(self.save_path,os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt,self.save_file_suffix,self.save_file_type)
            print(fname)
            plt.savefig(fname, bbox_inches='tight')

        if self.show_fig:
            plt.show()
        else:
            plt.close()

        return

    def plot_phase_fit_iteration2(self,model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
    data_filt,data_zero_err,data_small_err,data_gal,data_diff):

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
        gs = gridspec.GridSpec(2,1,height_ratios=[1,0.05])
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[1,0])

        # plot all the data from the SQL that goes into the fitting process
        ax1.errorbar(data['phase_angle'],data['reduced_mag'],data['merr'], fmt='ko',label="data",zorder=0,markersize="2")
        s1=ax1.scatter(data['phase_angle'],data['reduced_mag'],c=data['mjd'],s=10)
        cbar1=fig.colorbar(s1,ax2,use_gridspec=True, orientation='horizontal')
        cbar1.set_label("mjd")

        # plot iterative fits and cuts
        alpha_fit=np.linspace(np.amin(alpha),np.amax(alpha),100)
        print(label_iter_list)
        print(model_iter_list)
        # print(k)
        for j in range(len(model_iter_list)):
            print(j,label_iter_list[j])
            ax1.plot(alpha_fit,model_iter_list[j](alpha_fit),label=label_iter_list[j])

        ax1.plot(alpha_fit,model(alpha_fit),label=label)

        ax1.set_xlabel('alpha(degrees)')
        ax1.set_ylabel('reduced mag')
        ax1.invert_yaxis()
        ax1.legend(prop={'size': 6})

        ax1.set_title("{}_{}_{}_{}_{}".format(os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt))
        plt.tight_layout()

        if self.save_fig:
            fname="{}/{}_{}_{}_{}_{}_iter{}.{}".format(self.save_path,os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt,self.save_file_suffix,self.save_file_type)
            print(fname)
            plt.savefig(fname, bbox_inches='tight')

        if self.show_fig:
            plt.show()
        else:
            plt.close()

        return

    def plot_phase_fit_iteration_2panel(self,model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
    data_filt,data_zero_err,data_small_err,data_gal,data_diff):

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
        gs = gridspec.GridSpec(3,1,height_ratios=[1,1,0.1])
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[1,0])
        ax3 = plt.subplot(gs[2,0])

        # plot just the data that was fit
        ax2.errorbar(np.array(alpha),np.array(mag),np.array(mag_err), fmt='ko',label="data fit",zorder=0,markersize="2")
        s2=ax2.scatter(np.array(alpha),np.array(mag),c=np.array(data["mjd"]),label="data mjd",s=10)
        cbar2=fig.colorbar(s2,ax3,use_gridspec=True, orientation='horizontal')

        # plot all the data from the SQL that goes into the fitting process
        ax1.errorbar(data_filt['phase_angle'],data_filt['reduced_mag'],data_filt['merr'], fmt='ko',label="all data",zorder=0,markersize="2")

        # highlight any measurements with zero uncertainty
        ax1.scatter(data_zero_err['phase_angle'],data_zero_err['reduced_mag'],edgecolor='r',facecolor="none",marker="^",s=50,label="{} mag_err = 0".format(len(data_zero_err)))
        ax1.scatter(data_small_err['phase_angle'],data_small_err['reduced_mag'],edgecolor='r',facecolor="none",marker="s",s=50,label="{} mag_err < {}".format(len(data_small_err),self.mag_err_small))
        # highlight low galactic latitude
        ax1.scatter(data_gal['phase_angle'],data_gal['reduced_mag'],edgecolor='r',facecolor="none",marker="o",s=50,label="{} galactic_latitude < {} deg".format(len(data_gal),self.gal_lat_cut))

        if self.mag_diff_flag:
            #plot objects dropped in initial cut
            ax1.scatter(data_diff['phase_angle'],data_diff['reduced_mag'],edgecolor='r',facecolor="none",marker="p",s=50,label="{} HG model diff > {}".format(len(data_diff),self.mag_med_cut))

        # plot iterative fits and cuts
        alpha_fit=np.linspace(np.amin(alpha),np.amax(alpha),100)
        print(label_iter_list)
        print(model_iter_list)
        for j in range(len(model_iter_list)):
            print(j,label_iter_list[j])
            ax1.plot(alpha_fit,model_iter_list[j](alpha_fit),label=label_iter_list[j])
            ax1.scatter(alpha_cut_iter_list[j],mag_cut_iter_list[j],marker="x",zorder=3)

        ax2.plot(alpha_fit,model(alpha_fit),label=label,c="C{}".format(j+1))

        ax1.plot(alpha_fit,model(alpha_fit),label=label)

        ax1.set_xlabel('alpha(degrees)')
        ax1.set_ylabel('reduced mag')
        ax1.invert_yaxis()
        ax1.legend(prop={'size': 6})

        ax2.set_xlabel('alpha(degrees)')
        ax2.set_ylabel('reduced mag')
        ax2.invert_yaxis()
        ax2.legend(prop={'size': 6})

        ax1.set_title("{}_{}_{}_{}_{}".format(os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt))
        plt.tight_layout()

        if self.save_fig:
            fname="{}/{}_{}_{}_{}_{}_iter{}.{}".format(self.save_path,os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt,self.save_file_suffix,self.save_file_type)
            print(fname)
            plt.savefig(fname, bbox_inches='tight')

        if self.show_fig:
            plt.show()
        else:
            plt.close()

        return

    def plot_phase_fit_fancy(self,model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
    data_filt,data_zero_err,data_small_err,data_gal,data_diff):

        if not self.show_fig:
            import matplotlib
            print("use agg")
            matplotlib.use('agg') # use agg backend to stop python stealing focus when plotting

        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        # plt.rcParams.update({'font.size': 16})

        # extract asteroid phase data from data, with units
        alpha = np.array(data['phase_angle']) * u.deg
        mag = np.array(data["reduced_mag"]) * u.mag
        mag_err = np.array(data["merr"]) * u.mag

        fig = plt.figure()
        gs = gridspec.GridSpec(1,1)
        ax1 = plt.subplot(gs[0,0])

        # plot all the data from the SQL that goes into the fitting process
        ax1.errorbar(data['phase_angle'],data['reduced_mag'],data['merr'], fmt='ko',label="fitted data",zorder=0,markersize="2")
        # s1=ax1.scatter(data_filt['phase_angle'],data_filt['reduced_mag'],c=data_filt['mjd'],s=10)

        # highlight rejected data
        alpha_reject = np.concatenate([np.array(data_zero_err['phase_angle']),np.array(data_small_err['phase_angle']),np.array(data_gal['phase_angle'])])
        print(alpha_reject)
        mag_reject = np.concatenate([np.array(data_zero_err['reduced_mag']),np.array(data_small_err['reduced_mag']),np.array(data_gal['reduced_mag'])])

        if self.mag_diff_flag:
            #plot objects dropped in initial cut
            alpha_reject = np.append(alpha_reject,np.array(data_diff['phase_angle']))
            mag_reject = np.append(mag_reject,np.array(data_diff['reduced_mag']))

        # plot iterative fits and cuts
        alpha_fit=np.linspace(np.amin(alpha),np.amax(alpha),100)
        print(label_iter_list)
        print(model_iter_list)
        # print(k)
        for j in range(len(label_iter_list)):
            print(j,label_iter_list[j])
            alpha_reject = np.append(alpha_reject,np.array(alpha_cut_iter_list[j]))
            mag_reject = np.append(mag_reject,np.array(mag_cut_iter_list[j]))
            if j==0 or j==len(label_iter_list)-1:
                ax1.plot(alpha_fit,model_iter_list[j](alpha_fit),label=label_iter_list[j].split(". ")[-1])

        # ax1.scatter(alpha_reject,mag_reject,edgecolor='r',facecolor="none",marker="o",s=50,label="rejected data".format(len(mag_reject)))
        ax1.scatter(alpha_reject,mag_reject,c='r',s=10,marker = "x",label="rejected data".format(len(mag_reject)))

        ax1.set_xlabel('alpha(degrees)')
        ax1.set_ylabel('reduced mag')
        ax1.invert_yaxis()
        ax1.legend()

        ax1.set_title("{}_{}_{}_{}_{}".format(os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt))
        plt.tight_layout()

        if self.save_fig:
            fname="{}/{}_{}_{}_{}_{}_fancy{}.{}".format(self.save_path,os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt,self.save_file_suffix,self.save_file_type)
            print(fname)
            plt.savefig(fname, bbox_inches='tight')

        if self.show_fig:
            plt.show()
        else:
            plt.close()

        plt.style.use('default')

        return

    def plot_epochs(self,model_func,model_name,model,data,data_all_filt,epochs,filt):
        """
        """

        if not self.show_fig:
            import matplotlib
            print("use agg")
            matplotlib.use('agg') # use agg backend to stop python stealing focus when plotting

        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        alpha = np.array(data['phase_angle']) * u.deg
        mag = np.array(data["reduced_mag"]) * u.mag
        mag_err = np.array(data["merr"]) * u.mag
        residuals = mag - model(alpha)

        fig = plt.figure()
        gs = gridspec.GridSpec(2,1)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[1,0])

        ax1.axhline(0,c="k")
        # ax1.scatter(data["mjd"],residuals, label = "fitted data")
        ax1.set_xlabel("mjd")
        ax1.set_ylabel("O-C")
        ax2.set_xlabel("phase_angle(degrees)")
        ax2.set_ylabel("reduced_mag")

        ax2.plot(alpha,model(alpha),c="k")

        for i in range(len(epochs)):
            ax1.axvline(epochs[i],color="r")

        # find time difference between each detection
        # ascending order sort on date
        sort_mask = np.argsort(data["mjd"])
        mjd = np.array(data["mjd"])[sort_mask]
        residuals = np.array(residuals)[sort_mask]
        alpha = np.array(data["phase_angle"])[sort_mask]
        mag = np.array(data["reduced_mag"])[sort_mask]
        mag_err = np.array(data["merr"])[sort_mask]

        for i in range(1,len(epochs)):

            m1 = epochs[i-1]
            m2 = epochs[i]
            N_days_epoch = m2-m1
            date_mask = ((mjd>=m1) & (mjd<m2))
            N_data_epoch = len(mjd[date_mask])
            # print(m1,m2,N_days_epoch,N_data_epoch)

            ax1.scatter(mjd[date_mask],residuals[date_mask],c="C{}".format(i), label = "epoch {}".format(i))

            # try fit model to epoch data, provided there is enough data
            # if N_data_epoch>0 & len(alpha)>len(model_func.param_names):
            if N_data_epoch>len(model_func.param_names):

                res_med = np.median(residuals[date_mask])
                ax1.hlines(res_med,m1,m2,color="r")

                _alpha = np.array(alpha[date_mask]) * u.deg
                _mag = np.array(mag[date_mask]) * u.mag
                _mag_err = np.array(mag_err[date_mask]) * u.mag

                ax2.scatter(_alpha,_mag,c="C{}".format(i), label = "epoch {}".format(i))
                _model = self.fitter(model_func, _alpha, _mag, weights=1.0/np.array(_mag_err))
                ax2.plot(_alpha,_model(_alpha),c="C{}".format(i))

        # plot all data points
        ax1.scatter(data_all_filt["mjd"],data_all_filt["reduced_mag"]- np.array(model(np.array(data_all_filt["phase_angle"]) * u.deg)),
        s=1,c="k", label = "all data")

        ax1.invert_yaxis()
        ax2.invert_yaxis()

        ax1.set_title("{}_{}_{}_{}_{}_epochs".format(os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt))
        plt.tight_layout()

        if self.save_fig:
            fname="{}/{}_{}_{}_{}_{}_epochs{}.{}".format(self.save_path,os.path.basename(__file__).split('.')[0],self.file_identifier,model_name,self.clip_label,filt,self.save_file_suffix,self.save_file_type)
            print(fname)
            plt.savefig(fname, bbox_inches='tight')

        if self.show_fig:
            plt.show()
        else:
            plt.close()

        return

    def calculate(self):
        """calculate the phase curves on the phase_fit object"""

        if self.orbital_elements_id is None:
            print("object not in db, nothing to fit")
            return

        # get the object observation data, cutting on date range and dropping rows with nan
        data_all_filt=self.get_obj_data(self.cnx1,self.orbital_elements_id,self.start_date,self.end_date)

        # get the object metadata and combine with the phase fit dataframe structure
        df_obj=self.get_obj_meta(self.cnx1,self.orbital_elements_id,self.mpc_number,self.name)
        print(df_obj)
        d1=df_obj
        d2=self.df_obj_datafit
        df_obj=d2.append(d1) # add the values to the df

        # update the detection_count, dropping nans and outside date range
        df_obj["detection_count"]=len(data_all_filt) # note that atlas_objects may not have up to date detection count, call update_atlas_objects

        # set the mpc_number and name from df_obj !!!
        self.mpc_number=df_obj.iloc[0]['mpc_number']
        self.name=df_obj.iloc[0]['name']

        # Find the solar apparitions from elongation
        # ADD CONDITIONS ON ORBIT ON WHETHER OR NOT TO USE JPL?
        # or just ignore results in the database for q < 1.3 AU
        print(df_obj[["a_semimajor_axis","e_eccentricity","i_inclination_deg"]])
        orbital_period_yrs = df_obj.iloc[0]["a_semimajor_axis"]**1.5
        sol = sa.solar_apparitions(mpc_number = self.mpc_number, name = self.name, df_data = data_all_filt)
        epochs = sol.solar_elongation(-1.0,period = orbital_period_yrs)
        # epochs = sol.solar_elongation_JPL(JPL_step="7d")

        print(epochs)
        N_app = len(epochs)-1 # number of apparitions detected in both filters
        df_obj["N_apparitions"]=N_app


        # do a seperate fit for data in each filter
        for filt in self.filters:

            # retrieve astorb H and G values for the predicted fit
            G_slope=float(df_obj.iloc[0]['G_slope'])
            H_abs_mag=float(df_obj.iloc[0]['H_abs_mag'])
            print("filt = {}\nG_slope = {}\nH_abs_mag = {}".format(filt,G_slope,H_abs_mag))

            # do filter correction from V band (Heinze et al. 2020) - see also Erasmus et al 2020 for the c-o colours of S and C types (0.388 and 0.249 respectively)
            # if H_abs_mag_o/c can be used to manually set the filter correction
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

            # We could fit an object multiple times, using the previous H as the new starting point?
            # use the updated field to track how many times this iteration has been done?

            # select all data from a certain filter
            data_filt=data_all_filt[data_all_filt['filter']==filt]

            # Record number of detections in the filter AFTER nan and date cuts but BEFORE other cuts are made
            # This updates the detection_count value from rockatlas
            N_data_start = len(data_filt)
            detection_count_filt=N_data_start
            df_obj["detection_count_{}".format(filt)]=detection_count_filt # update the number of detections in the df

            # Also update the rockatlas value of phase angle range, before fits are done
            df_obj["phase_angle_range_{}".format(filt)] = np.amax(data_filt["phase_angle"]) - np.amin(data_filt["phase_angle"])

            # cut starting data for this filter
            print("{} starting data".format(len(data_filt)))

            # drop any measurements with zero uncertainty
            mask_zero = data_filt['merr']==0
            data_zero_err=data_filt[mask_zero]
            data_filt=data_filt[~mask_zero]
            print("{} after zero error cut".format(len(data_filt)))
            # drop measurements with small (or zero) uncertainty
            mask_err = data_filt['merr']<self.mag_err_small
            data_small_err=data_filt[mask_err]
            data_filt=data_filt[~mask_err]
            print("{} after small error cut".format(len(data_filt)))
            # drop measurements near galactic plane
            mask_gal = np.absolute(data_filt["galactic_latitude"])<self.gal_lat_cut
            data_gal=data_filt[mask_gal]
            data_filt=data_filt[~mask_gal]
            print("{} after galactic plane cut".format(len(data_filt)))

            print("{} data after cuts".format(len(data_filt)))

            # RECORD THE NUMBER OF DATA POINTS THAT HAVE BEEN CUT
            N_data_zero_err = len(data_zero_err)
            N_data_small_err = len(data_small_err)
            N_data_gal = len(data_gal)
            N_data_cut = N_data_zero_err + N_data_small_err + N_data_gal
            print("CUT data_zero_err = {}\nCUT data_small_err = {}\nCUT data_gal = {}".format(
            N_data_zero_err,N_data_small_err,N_data_gal))
            print("TOTAL CUT N_data_cut = {}".format(N_data_cut))

            # if no data remains after cuts, then nothing can be fit
            if len(data_filt)==0:
                print("no data, cannot fit")
                break

            # iterate over all models
            for model_name,model_values in self.selected_models.items():

                # cut on phase angle range only if using the Linear phase function
                # MAKE THIS OPTIONAL?
                if model_name=="LinearPhaseFunc":
                    mask_alpha = ((data_filt["phase_angle"]>=self.phase_lin_min) & (data_filt["phase_angle"]<=self.phase_lin_max))
                    N_alpha_lin=len(data_filt[~mask_alpha])
                    data_filt=data_filt[mask_alpha]
                    print("CUT N_alpha_lin = {}".format(N_alpha_lin))
                else:
                    N_alpha_lin=0

                print(model_name,model_values)

                # retrieve model names etc
                ms=model_values["model_name_short"]
                pc=model_values["model_parameters"]
                model_func=model_values["model_function"]

                print(ms,pc,model_func)

                # Set the starting params with some dummy values, old_params will be used to test if the solution has converged
                old_params=[999]*len(model_func.parameters)

                # store models/labels/cut data for plotting
                label_iter_list=[]
                model_iter_list=[]
                mag_cut_iter_list=[]
                alpha_cut_iter_list=[]

                print("{}, {}: fit {}, filter {}".format(self.name,self.mpc_number,model_name,filt))

                # initialise the data that we will iteratively fit and cut
                data=data_filt
                data=data.sort_values("phase_angle") # ensure that the dataframe is in order for plotting
                # data=data.sort_values("mjd") # ensure that the dataframe is in date order for finding epochs

                # iteratively fit and cut data, for a maximum of max_iters times
                k=0
                while k<self.max_iters:

                    print("iteration: {}, N_data={}".format(k,len(data)))

                    # extract asteroid phase data from data, with units
                    alpha = np.array(data['phase_angle']) * u.deg
                    mag = np.array(data["reduced_mag"]) * u.mag
                    mag_err = np.array(data["merr"]) * u.mag

                    if k==0:
                        # for first iteration start with the predicted HG mag (Bowell 1989)
                        model=HG(H = H_abs_mag * u.mag, G = G_slope)
                        model_str="astorb HG"
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
                        # check that there is still enough data to fit (need more data points than parameters)
                        if len(alpha)<=len(pc):
                            print("less data to fit than parameters")
                            break

                        # !!! RECORD ANY WARNINGS ETC? see fitter.fit_info
                        # logging will record these warnings
                        with warnings.catch_warnings(record=True) as w:
                            model = self.fitter(model_func, alpha, mag, weights=1.0/np.array(mag_err)) # fit using weights by uncertainty
                            if len(w)>0:
                                warning_message="{} - {} - {} - {} - {}".format(self.mpc_number,self.name,model_name,filt,w[-1].message)
                                print(warning_message)
                                logging.warning(warning_message)

                        param_names=model.param_names
                        params=model.parameters

                        # label each fit iteration
                        labels=["{}. {}".format(k,model_name)]
                        for l in range(len(model.parameters)):
                            labels.append("{}={:.2f} ".format(param_names[l],params[l]))
                        label=", ".join(labels)

                        # test for convergence, find difference in params
                        delta_params = np.absolute(old_params-params)
                        # if difference is less than threshold, we deem the solution to have converged and stop
                        if np.sum(delta_params<self.param_converge_check)==len(delta_params):
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

                            # retrieve errors in the parameters: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
                            param_err_x = np.sqrt(np.diag(param_cov))
                            param_shape=param_cov.shape
                            correl_coef=np.zeros(param_shape)
                            for l in range(param_shape[0]):
                                for m in range(param_shape[1]):
                                    correl_coef[l,m]=param_cov[l,m]/np.sqrt(param_cov[l,l]*param_cov[m,m]) # Hughes AND Hase eqn. 7.3

                            # retrieve metrics for the data
                            N_data_fit=len(data) # number of data points fit after clipping
                            N_nights=len(np.unique(np.array(data["mjd"]).astype(int))) # number of unique nights in data set
                            alpha_min=np.amin(alpha).value # minimum phase angle in data that was fit
                            alpha_max=np.amax(alpha).value # max phase angle
                            N_alpha_low=len(alpha[alpha<self.low_alpha_cut]) # Number of data points at low phase angle
                            N_iter=k # number of fit and clip iterations
                            nfev=self.fitter.fit_info['nfev'] # number of scipy function evalutations
                            ier=self.fitter.fit_info['ierr'] # scipy fitter flag
                            N_mag_err=len(mag_err[np.array(mag_err)<self.mag_err_threshold]) # Use leq? Number of data points with error below some threshold

                            # Record the number of data points cut during the fitting process
                            # Note that data_filt will have been cut by phase angle when the Linear function is fit
                            N_data_cut = len(data_filt) - N_data_fit
                            print("N_data_fit = {}\nN_data_cut = {}".format(N_data_fit,N_data_cut))

                            # check number of data points cut and fit all add up to starting data
                            N_data_tot = N_data_zero_err + N_data_small_err + N_data_gal + N_alpha_lin + N_data_cut + N_data_fit
                            print(N_data_zero_err,N_data_small_err,N_data_gal,N_alpha_lin,N_data_cut,N_data_fit)
                            print("N_data_tot = {}\nN_data_start = {}".format(N_data_tot,N_data_start))
                            if N_data_tot!=N_data_start:
                                print("Error, number of data points doesn't add up")
                                logging.warning("{} - {} - {} - {} - N_data_tot!=N_data_start".format(self.mpc_number,self.name,model_name,filt))

                            # calculate the residual properties
                            residuals = mag - model(alpha)
                            OC_mean = np.mean(residuals)
                            OC_std = np.std(residuals)
                            OC_range = np.absolute(np.amax(residuals)-np.amin(residuals))

                            # residuals for each epoch
                            sort_mask = np.argsort(data["mjd"])
                            mjd = np.array(data["mjd"])[sort_mask]
                            residuals = np.array(residuals)[sort_mask]
                            res_med_list = []
                            for i in range(1,len(epochs)):
                                m1 = epochs[i-1]
                                m2 = epochs[i]
                                N_days_epoch = m2-m1
                                date_mask = ((mjd>=m1) & (mjd<m2))
                                N_data_epoch = len(mjd[date_mask])
                                if N_data_epoch>len(pc): # need at least the same of data points as number of parameters
                                    res_med = np.median(residuals[date_mask])
                                    res_med_list.append(res_med)
                                    # print(m1,m2,N_days_epoch,N_data_epoch,res_med)
                            res_med_list = np.array(res_med_list)
                            print(res_med_list)
                            app_res_med = np.median(res_med_list) # median of the median residual for all apparitions
                            app_res_std = np.std(res_med_list) # std of the median residual for all apparitions
                            app_res_mean = np.mean(res_med_list) # mean of the median residual for all apparitions
                            app_res_range = np.absolute(np.amax(res_med_list)-np.amin(res_med_list)) # maximum difference between apparition residuals

                            print("total number of epochs = {}".format(N_app))
                            print("median median epoch residual = {}".format(app_res_med))
                            print("mean median epoch residual = {}".format(app_res_mean))
                            print("std median epoch residual = {}".format(app_res_std))
                            print("range median epoch residual = {}".format(app_res_range))

                            # check for errors
                            if N_mag_err>N_data_fit:
                                # ERROR in calculating N-mag_err?
                                print("N_mag_err={}".format(N_mag_err))
                                logging.warning("{} - {} - {} - {} - N_mag_err>N_data_fit".format(self.mpc_number,self.name,model_name,filt))

                            # ALSO check for any nans in metrics? e.g. OC?

                            print(self.fitter.fit_info['message'])

                            print(df_obj["detection_count_{}".format(filt)])

                            # populate the df_obj dataframe: add all fit params/metrics to df_obj
                            df_obj["phase_curve_N_fit{}_{}".format(ms,filt)]=N_data_fit
                            df_obj["phase_curve_alpha_min{}_{}".format(ms,filt)]=alpha_min
                            df_obj["phase_curve_alpha_max{}_{}".format(ms,filt)]=alpha_max
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

                            df_obj["phase_curve_N_cut{}_{}".format(ms,filt)]=N_data_cut

                            df_obj["phase_curve_app_res_med{}_{}".format(ms,filt)]=app_res_med
                            df_obj["phase_curve_app_res_std{}_{}".format(ms,filt)]=app_res_std
                            df_obj["phase_curve_app_res_range{}_{}".format(ms,filt)]=app_res_range

                            for p in range(len(pc)):
                                df_obj["phase_curve_{}{}_{}".format(pc[p],ms,filt)]=params[p]
                                df_obj["phase_curve_{}_err{}_{}".format(pc[p],ms,filt)]=param_err_x[p]

                            # might want to record a different set of parameters in df obj for LinearPhaseFunc
                            # if model_name=="LinearPhaseFunc":
                            # if "LinearPhaseFunc" in model_list:
                                # Set only the LinearPhaseFunc parameters

                            if self.plot_fig:

                                # Plot this fit (different plotting functions available)

                                # self.plot_phase_fit(model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
                                # data_filt,data_zero_err,data_small_err,data_gal)

                                # self.plot_phase_fit_iteration(model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
                                # data_filt,data_zero_err,data_small_err,data_gal,data_diff)

                                # self.plot_phase_fit_iteration2(model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
                                # data_filt,data_zero_err,data_small_err,data_gal,data_diff)

                                self.plot_phase_fit_iteration_2panel(model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
                                data_filt,data_zero_err,data_small_err,data_gal,data_diff)

                                # self.plot_phase_fit_fancy(model,model_name,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
                                # data_filt,data_zero_err,data_small_err,data_gal,data_diff)

                                # plot epochs
                                self.plot_epochs(model_func,model_name,model,data,data_all_filt,epochs,filt)
                                # exit()

                            # # save data that was used to fit to file
                            # data_clip_file="results_analysis/fit_data/df_data_{}{}_{}.csv".format(self.mpc_number,ms,filt)
                            # print(data_clip_file)
                            # data.to_csv(data_clip_file)

                            break

                    # if params did not converge we clip outliers and repeat

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

        # push df_obj to the database
        if self.push_fit==True:

            # Will this fail if the names in df_obj do not match fields in mysql db?
            # push_obj_db SHOULD only push selected columns in db_columns
            # ADD LinearPhaseFunc fields to db!

            # catch a weird error
            if df_obj.iloc[0]["phase_curve_N_mag_err_B89_o"]>df_obj.iloc[0]["phase_curve_N_fit_B89_o"]:
                logging.warning("{} - {} - {} - {} - N_mag_err>N_data_fit".format(self.mpc_number,self.name,model_name,filt))

            self.push_obj_db(df_obj)

        # close all database connections
        # CLOSE CURSORS TOO?
        self.cnx1.disconnect()
        self.cnx2.disconnect()

        return df_obj # return the fit df
