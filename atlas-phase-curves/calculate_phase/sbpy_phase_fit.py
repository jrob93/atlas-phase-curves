# https://sbpy.readthedocs.io/en/latest/sbpy/photometry.html#disk-integrated-phase-function-models

import numpy as np
import astropy.units as u
# from astropy import units as u
from astropy.modeling.fitting import LevMarLSQFitter
from sbpy.photometry import HG, HG1G2, HG12, HG12_Pen16
import pandas as pd
import os
import mysql.connector
# from optparse import OptionParser
# import subprocess
# import sys
import datetime
# import time

# # import our classes
if __name__ == "__main__": # import other modules as a script
    from atlas_SQL_query_df import atlas_SQL_query
    from atlas_database_connection import database_connection
else: # import as a package
    from calculate_phase.atlas_SQL_query_df import atlas_SQL_query
    from calculate_phase.atlas_database_connection import database_connection

class phase_fit():

    # default values
    low_alpha_cut=5.0*u.deg # we want to quantify how many data points are fit at low phase angles, alpha < low_alpha_cut
    param_converge_check=0.01 # the model is fit until the change in parameters (e.g. H and G) is less than param_converge_check (or max iterations is reached)
    max_iters=30 # maximum number of attempts at fitting and cutting
    std=2 # standard deviation of the sigma data clip
    mag_err_threshold = 0.1 # limit for the error of "good" data, we record N_mag_err number of data points with error < mag_err_threshold
    mag_err_small = 0.01 # we discount observations with error less than this
    gal_lat_cut=10 # galatic latitude cut in degrees

    # set the clipping method for the phase_fit class
    def data_clip_sigma(self,data,data_predict,low=std,high=std):
        """ function to cut outliers by standard deviation, i.e. sigma clipping"""
        std=np.std(data)
        clip_mask=((data < (data_predict-(std*low))) | (data > (data_predict+(std*high))))
        return clip_mask

    def data_clip_diff(self,data,data_predict,diff=1):
        # cut outliers by diff (this function doesn't like astropy units, use np arrays)
        clip_mask=(np.absolute(data_predict-data)>diff)
        return clip_mask

    # define which clipping function to use
    data_clip_func=data_clip_sigma
    clip_label="{}-sigma_clip".format(std)

    # set up the sbpy fitter and models
    fitter = LevMarLSQFitter()
    model_names_str = ["HG", "HG1G2", "HG12", "HG12_Pen16"]
    model_short = ["_B89","_3M10","_2M10","_P16"] # use shorthand for all models
    phase_curve = [["H","G"],["H","G1","G2"],["H","G12"],["H","G12"]]
    model_names = [HG(), HG1G2(), HG12(), HG12_Pen16()]

    def __init__(self,
        mpc_number,
        save_path=".",
        push_fit_flag=False,plot_fig_flag=False,show_fig_flag=False,save_fig_flag=False,
        start_date=False,end_date=False,
        filter_list=["o","c"]):

        # set up the class
        # define object and some flags

        self.mpc_number=mpc_number
        self.push_fit=push_fit_flag # flag to push fit to database
        self.plot_fig=plot_fig_flag # flag to generate plot for each object
        self.show_fig=show_fig_flag # flag to display interactive plot
        self.save_fig=save_fig_flag # flag to save the figure
        self.utc_date_now=datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S") # time at which the script is run (UTC)
        self.save_path=save_path # where to save figures

        self.start_date=start_date
        self.end_date=end_date
        # self.date_hard_cap=? # set a fixed date beyond which we do not consider data
        self.filters=filter_list # list of filters to use: ["o"], ["c"], ["o","c"]

        # cnx1 is the connection where we retrieve the data from
        self.cnx1=database_connection().connect()
        self.cursor1 = self.cnx1.cursor()

        # cnx2 is the connection where we push the fit data to
        self.cnx2 = self.cnx1
        self.cursor2 = self.cursor1

        self.tab_name="atlas_phase_fits" # table name within the database

    def get_obj_data(self,cnx,mpc_num,t_start=False,t_end=False):
        # load data to be fitted, loads both filters (o & c)
        data_all_filt=atlas_SQL_query(cnx=cnx,mpc_number=mpc_num)

        # cut by date range. Note that different filters may have different dates/epochs!
        if t_start:
            data_all_filt=data_all_filt[data_all_filt['mjd']>float(t_start)]
        if t_end:
            data_all_filt=data_all_filt[data_all_filt['mjd']<float(t_end)]
        # else: # if t_end is not defined default to the hard cap for date
        #     data_all_filt=data_all_filt[data_all_filt['mjd']<self.date_hard_cap]

        return data_all_filt

    def get_obj_metadata(self,cnx,mpc_num):

        # get the last bits of data. might be out of date if rockatlas isn't running...?
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
        FROM atlas_objects WHERE mpc_number=%(mpc_num)s
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

    def get_obj_HG(self,cnx,mpc_num):
        # get the astorb H and G (B89)
        # qry_HG="select G_slope, H_abs_mag from orbital_elements where primaryId='{}';".format(self.mpc_number)
        qry_HG="select G_slope, H_abs_mag from orbital_elements where mpc_number='{}';".format(mpc_num)
        # print(qry_HG)
        df_HG=pd.read_sql_query(qry_HG,cnx)
        return df_HG

    def get_obj_meta(self,cnx,mpc_num):

        # get the last bits of data. might be out of date if rockatlas isn't running...?
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
        FROM atlas_objects a, orbital_elements o WHERE a.mpc_number=%(mpc_num)s AND o.mpc_number=%(mpc_num)s
        """ % locals()
        # print(qry_obj1)

        df_obj=pd.read_sql_query(qry_obj1,cnx)
        # print(df_obj)
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

    def database_obj_entry(self,mpc_num,df_obj,mpc_check):

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

        print(qry)
        # exit()

        self.cursor2.execute(qry)
        self.cnx2.commit()

    def plot_phase_fit(self,model,model_str,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
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

        ax1.set_title("{}_{}_{}_{}_{}".format(os.path.basename(__file__).split('.')[0],self.mpc_number,model_str,self.clip_label,filt))
        plt.tight_layout()

        if self.save_fig:
            fname="{}/{}_{}_{}_{}_{}.png".format(self.save_path,os.path.basename(__file__).split('.')[0],self.mpc_number,model_str,self.clip_label,filt)
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
        data_all_filt=self.get_obj_data(self.cnx1,self.mpc_number,self.start_date,self.end_date)
        print(data_all_filt)

        # get the object metadata
        # df_obj=self.get_obj_metadata(self.cnx1,self.mpc_number)
        df_obj=self.get_obj_meta(self.cnx1,self.mpc_number)
        print(df_obj)

        # update the detection_count
        df_obj["detection_count"]=len(data_all_filt) # note that atlas_objects may not have up to date detection count...

        # # get the astorb H and G (B89)
        # # qry_HG="select G_slope, H_abs_mag from orbital_elements where primaryId='{}';".format(self.mpc_number)
        # qry_HG="select G_slope, H_abs_mag from orbital_elements where mpc_number='{}';".format(self.mpc_number)
        # print(qry_HG)
        # df_HG=pd.read_sql_query(qry_HG,self.cnx1)
        # df_HG=self.get_obj_HG(self.cnx1,self.mpc_number)
        # print(df_HG)
        # exit()

        # is the object already in the tab_name at the storage database?
        self.cursor2.execute("SELECT COUNT(1) FROM {} where mpc_number={}".format(self.tab_name,self.mpc_number))
        mpc_check=self.cursor2.fetchone()[0]
        print(mpc_check)

        if self.push_fit==True: # only update the entry if we are going to push results
            # MOVE THIS TO PUSH FIT VALUES
            qry_obj = self.database_obj_entry(self.mpc_number,df_obj,mpc_check)
            # qry_obj = self.database_obj_entry2(self.mpc_number,df_obj,mpc_check)
            print(qry_obj)
            self.cursor2.execute(qry_obj)
            self.cnx2.commit()
            # exit()

        for filt in self.filters:

            # H and G values for the predicted fit
            # G_slope=float(df_HG.iloc[0]['G_slope'])
            # H_abs_mag=float(df_HG.iloc[0]['H_abs_mag'])
            G_slope=float(df_obj.iloc[0]['G_slope'])
            H_abs_mag=float(df_obj.iloc[0]['H_abs_mag'])
            # print("G_slope = {}\nH_abs_mag = {}".format(G_slope,H_abs_mag))

            # do filter correction from V band (Heinze et al. 2020)
            if filter=="o":
                H_abs_mag+=-0.332
            if filter=="c":
                H_abs_mag+=0.054

            # select all data from a certain filter
            data_filt=data_all_filt[data_all_filt['filter']==filt]
            # print(data_filt)
            detection_count_filt=len(data_filt) # !!! SHOULD WE RECORD NUMBER OF DETECTIONS BEFORE/AFTER CUTS?
            # print(detection_count_filt)
            # continue

            # iterate over all models
            for i,model_name in enumerate(self.model_names):

                ms=self.model_short[i]
                pc=self.phase_curve[i]

                # if ms!=model_short[1]:
                #     continue

                old_params=[999]*len(model_name.parameters)

                # store models/labels/cut data for plotting
                label_iter_list=[]
                model_iter_list=[]
                mag_cut_iter_list=[]
                alpha_cut_iter_list=[]

                print("{}: fit {}, filter {}".format(self.mpc_number,self.model_names_str[i],filt))

                # initialise the data that we will iteratively fit and cut
                data=data_filt
                data=data.sort_values("phase_angle") # ensure that the dataframe is in order for plotting

                # drop any measurements with zero uncertainty
                data_zero_err=data[data['merr']==0]
                data=data[data['merr']!=0]
                # drop measurements with small (or zero) uncertainty MOVE THIS TO A VARIABLE UP TOP!
                data_small_err=data[data['merr']<self.mag_err_small]
                data=data[~(data['merr']<self.mag_err_small)]

                # drop measurements near galactic plane
                data_gal=data[np.absolute(data["galactic_latitude"])<self.gal_lat_cut]
                data=data[~(np.absolute(data["galactic_latitude"])<self.gal_lat_cut)]

                if len(data)==0:
                    print("no data, cannot fit")
                    break

                # iteratively fit and cut data
                k=0
                while k<self.max_iters:

                    # print("iteration: {}".format(k))
                    print("iteration: {}, N_data={}".format(k,len(data)))

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

                    else:
                        # apply the fitter to the data
                        if len(alpha)<=len(pc): # check that there is enough data to fit
                            print("less data to fit than parameters")
                            break
                        model = self.fitter(model_name, alpha, mag, weights=1.0/np.array(mag_err)) # fit using weights by uncertainty
                        # if k==1:
                        #     model = fitter(model_name, alpha, mag) # drop the weights for the first fit
                        # else:
                        #     model = fitter(model_name, alpha, mag, weights=1.0/np.array(mag_err))

                        model_str=self.model_names_str[i]

                        # !!! RECORD ANY WARNINGS ETC? see fitter.fit_info

                        param_names=model.param_names
                        params=model.parameters
                        # print(model.__dict__)
                        # print()
                        # print(fitter.fit_info['cov_x']) # The scipy.optimize.leastsq result for the most recent fit https://docs.astropy.org/en/stable/api/astropy.modeling.fitting.LevMarLSQFitter.html
                        # print(err_x)

                        # label each fit
                        labels=["{}. {}".format(k,model_str)]
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
                            N_mag_err=len(mag_err[np.array(mag_err)<self.mag_err_threshold])

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

                            if self.push_fit==True:

                                self.push_fit_db(params,param_err_x,pc,ms,filt,detection_count_filt,N_data_fit,alpha_min,alpha_max,phase_angle_range,N_alpha_low,N_nights,N_iter,nfev,ier,N_mag_err)

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

                                self.plot_phase_fit(model,model_str,filt,label,data,label_iter_list,model_iter_list,alpha_cut_iter_list,mag_cut_iter_list,
                                data_filt,data_zero_err,data_small_err,data_gal)

                            # exit()
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

        # close all database connections
        # CLOSE CURSORS TOO?
        self.cnx1.disconnect()
        self.cnx2.disconnect()

        return

if __name__ == "__main__":
    mpc_number=4985
    # fit = phase_fit(mpc_number,push_fit_flag=True,plot_fig_flag=True,save_fig_flag=True)
    fit = phase_fit(mpc_number,push_fit_flag=True)#,plot_fig_flag=True,save_fig_flag=True)
    # print(fit.mpc_number)
    # print(fit.low_alpha_cut)
    # print(fit.__dict__)
    # # data=np.random.rand(10)
    # # data2=np.zeros(10)+np.median(data)
    # # print(data,np.median(data))
    # # print(data2)
    # # print(fit.data_clip_func(data,data2,fit.std,fit.std))
    # print(fit.clip_label)
    # print(fit.cnx1)
    # print(fit.tab_name)

    fit.calculate()

    exit()

"""parser = OptionParser()
parser.add_option( "-n", "--mpc-number", dest="mpc_number", help="mpc_number", metavar="MPC_NUMBER" ) # mpc number of object to fit
parser.add_option( "-s", "--save-path", dest="save_path", help="save_path", metavar="SAVE_PATH" ) # path to save figures
parser.add_option( "-f", "--filter", dest="filter", help="filter", metavar="FILTER" ) # filters to use, declare o or c, default is both
parser.add_option( "-w", "--warnings", dest="warnings", help="warnings", metavar="WARNINGS" ) # Suppress warnings?
parser.add_option( "-i", "--start-date", dest="start_date", help="start_date", metavar="START_DATE" ) # starts date for data
parser.add_option( "-j", "--end-date", dest="end_date", help="end_date", metavar="END_DATE" ) # starts date for data

(options,args)=parser.parse_args()

if options.mpc_number:
    mpc_number=int(options.mpc_number)
else:
    mpc_number="4986"
if options.save_path:
    save_path=options.save_path
else:
    save_path="."
if options.filter:
    filters=list(options.filter)
else:
    filters=["o","c"]
if options.warnings:
    warning_flag=int(options.warnings)
else:
    warning_flag=0
# if options.start_date:
#     start_date=int(start_date)
# if options.end_date:
#     end_date=int(end_date)

# Suppress warnings if the output is annoying. Be very careful suppressing warnings...
if warning_flag==1:
    import warnings
    warnings.filterwarnings('ignore')
"""
# print(filters)

# # define object and some flags
# obj_number=mpc_number
# print(mpc_number)
# low_alpha_cut=5.0*u.deg # we want to quantify how many data points are fit at low phase angles, alpha < low_alpha_cut
# param_converge_check=0.01 # the model is fit until the change in parameters (e.g. H and G) is less than param_converge_check (or max iterations is reached)
# max_iters=30 # maximum number of attempts at fitting and cutting
# std=2 # standard deviation of the sigma data clip
# mag_err_threshold = 0.1 # limit for the error of "good" data, we record N_mag_err number of data points with error < mag_err_threshold
# mag_err_small = 0.01 # we discount observations with error less than this
# gal_lat_cut=10 # galatic latitude cut in degrees
# push_fit=False # flag to push fit to database
# plot_fig=False # flag to generate plot for each object
# show_fig=False # flag to display interactive plot
# save_fig=False # flag to save the figure
# # push_fit=False # flag to push fit to database
# # plot_fig=True # flag to generate plot for each object
# # show_fig=False # flag to display interactive plot
# # save_fig=True # flag to save the figure

# if not show_fig:
#     import matplotlib
#     matplotlib.use('agg') # use agg backend to stop python stealing focus when plotting

# utc_date_now=datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S") # time at which the script is run (UTC)

# def data_clip_sigma(data,data_predict,low=2,high=2):
#     # cut outliers by sigma
#     std=np.std(data)
#     clip_mask=((data < (data_predict-(std*low))) | (data > (data_predict+(std*high))))
#     return clip_mask
#
# # def data_clip_diff(data,data_predict,diff=1):
# #     # cut outliers by diff (this function doesn't like astropy units, use np arrays)
# #     clip_mask=(np.absolute(data_predict-data)>diff)
# #     return clip_mask

# data_clip_func=data_clip_sigma
# clip_label="{}-sigma_clip".format(std)
# # clip_label="{}-mag_diff_clip".format(diff)

# # fit new models to model data
# fitter = LevMarLSQFitter()
# model_names_str = ["HG", "HG1G2", "HG12", "HG12_Pen16"]
# model_short = ["_B89","_3M10","_2M10","_P16"] # use shorthand for all models
# phase_curve = [["H","G"],["H","G1","G2"],["H","G12"],["H","G12"]]
#
# model_names = [HG(), HG1G2(), HG12(), HG12_Pen16()]

# # Predicted mag using parameters from orbital_elements table (USE ONLY FOR THE 1-MAG DIFF CUT?)
# # connect to the database (must be connected to dormammu: sshpass -e ssh dormammu) to read from
# config1 = {
#   'user': 'af',
#   'password': 'afPass',
#   'host': 'localhost',
#   'port': '3308',
#   'database': 'atlas_moving_objects',
#   'raise_on_warnings': True
# }
# cnx1 = mysql.connector.connect(**config1)
# cursor1 = cnx1.cursor()
#
# # connect to the database to write to
# # !!! too many connections! use the same connection but switch db in qry?
# # config2 = {
# #   'user': 'af',
# #   'password': 'afPass',
# #   'host': 'localhost',
# #   'port': '3308',
# #   'database': 'jamie_test_db',
# #   'raise_on_warnings': True
# # }
# # cnx2 = mysql.connector.connect(**config2)
# # cursor2 =cnx2.cursor()
# # tab_name="test_table2"
#
# # !!! merge conections when we start writing to the same db, try:
# cnx2 = cnx1
# cursor2 = cursor1
# tab_name="atlas_phase_fits"

# # get the last bits of data. might be out of date if rockatlas isn't running...?
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
# FROM atlas_objects WHERE mpc_number=%(mpc_number)s
# """ % locals()
# # print(qry_obj1)
#
# df_obj=pd.read_sql_query(qry_obj1,cnx1)
# # print(pd.Series(df_obj.iloc[0]))
# # print(df_obj.to_string())
# # print(df_obj[['detection_count',  'detection_count_c',  'detection_count_o']])
#
# # fix the date strings
# dates=["dateLastModified","last_photometry_update_date_c","last_photometry_update_date_o"]
# for d in dates:
#     if df_obj[d].iloc[0] is not None:
#         df_obj.loc[0,d]="'{}'".format(str(df_obj[d].iloc[0]))
# # for some reason None needs changed to Null?
# df_obj=df_obj.fillna(value="NULL")
# # print(df_obj.to_string())
# # print(df_obj[['detection_count',  'detection_count_c',  'detection_count_o']])

# qry_HG="select G_slope, H_abs_mag from orbital_elements where primaryId='{}';".format(obj_number)
# df_HG=pd.read_sql_query(qry_HG,cnx1)

# cursor2.execute("SELECT COUNT(1) FROM {} where mpc_number={}".format(tab_name,mpc_number))
# mpc_check=cursor2.fetchone()[0]
# print("mpc_check = {}".format(mpc_check))

# # load data to be fitted, loads both filters (o & c)
# data_all_filt=atlas_SQL_query(cnx=cnx1,mpc_number=obj_number)
# detection_count=len(data_all_filt) # note that atlas_objects may not have up to date detection count...
#
# # start = time.process_time()
# # cut by date range. Note that different filters may have different dates/epochs!
# if options.start_date:
#     data_all_filt=data_all_filt[data_all_filt['mjd']>float(options.start_date)]
# if options.end_date:
#     data_all_filt=data_all_filt[data_all_filt['mjd']<float(options.end_date)]
# # print(time.process_time() - start)

# print(data_all_filt)
# print(list(data_all_filt))
# exit()

# if mpc_check==0:
#
#     qry_obj = u"""INSERT INTO {}
#     (mpc_number,
#     dateLastModified,
#     detection_count,
#     last_detection_mjd,
#     last_photometry_update_date_c,
#     last_photometry_update_date_o,
#     name,
#     orbital_elements_id,
#     primaryId)
#     VALUES
#     (%(mpc_number)s,{},{},{},{},{},"{}",{},{});""".format(tab_name,
#     str(df_obj['dateLastModified'].iloc[0]),
#     detection_count,
#     float(df_obj['last_detection_mjd']),
#     str(df_obj['last_photometry_update_date_c'].iloc[0]),
#     str(df_obj['last_photometry_update_date_o'].iloc[0]),
#     str(df_obj['name'].iloc[0]),
#     int(df_obj['orbital_elements_id']),
#     int(df_obj['primaryId'])) % locals()
#
#     # !!! insert detection count etc here?
#
# else:
#     qry_obj = u"""UPDATE {} SET
#     dateLastModified={},
#     detection_count={},
#     last_detection_mjd={},
#     last_photometry_update_date_c={},
#     last_photometry_update_date_o={},
#     name="{}",
#     orbital_elements_id={},
#     primaryId={}
#     WHERE mpc_number=%(mpc_number)s;""".format(tab_name,
#     str(df_obj['dateLastModified'].iloc[0]),
#     detection_count,
#     float(df_obj['last_detection_mjd']),
#     str(df_obj['last_photometry_update_date_c'].iloc[0]),
#     str(df_obj['last_photometry_update_date_o'].iloc[0]),
#     str(df_obj['name'].iloc[0]),
#     int(df_obj['orbital_elements_id']),
#     int(df_obj['primaryId'])) % locals()
#
# # !!! detection count is NaN a lot...? CHECK METRICS
# # !!! use detection_count to track if object should be refit? no new data, no refit...
# # !!! or use last photometry update date
#
# # print(qry_obj)
#
# # !!! WHAT TO DO WITH updated FLAG? use it to track number of times the fit has been done? Change to a date?
#
# if push_fit==True: # only update the entry if we are going to push results
#     cursor2.execute(qry_obj)
#     cnx2.commit()

# for filt in filters:
#
#     # H and G values for the predicted fit
#     G_slope=float(df_HG.iloc[0]['G_slope'])
#     H_abs_mag=float(df_HG.iloc[0]['H_abs_mag'])
#     # print("G_slope = {}\nH_abs_mag = {}".format(G_slope,H_abs_mag))
#
#     # do filter correction from V band (Heinze et al. 2020)
#     if filter=="o":
#         H_abs_mag+=-0.332
#     if filter=="c":
#         H_abs_mag+=0.054
#
#     # select all data from a certain filter
#     data_filt=data_all_filt[data_all_filt['filter']==filt]
#     # print(data_filt)
#     detection_count_filt=len(data_filt) # !!! SHOULD WE RECORD NUMBER OF DETECTIONS BEFORE/AFTER CUTS?
#     # print(detection_count_filt)
#     # continue
#
#     # iterate over all models
#     for i,model_name in enumerate(model_names):
#
#         ms=model_short[i]
#
#         # if ms!=model_short[1]:
#         #     continue
#
#         old_params=[999]*len(model_name.parameters)
#
#         # store models/labels/cut data for plotting
#         label_iter_list=[]
#         model_iter_list=[]
#         mag_cut_iter_list=[]
#         alpha_cut_iter_list=[]
#
#         print("{}: fit {}, filter {}".format(mpc_number,model_names_str[i],filt))
#
#         # initialise the data that we will iteratively fit and cut
#         data=data_filt
#         data=data.sort_values("phase_angle") # ensure that the dataframe is in order for plotting
#
#         # make a cut by date to match Dave's fits. Add apparition effect stuff here?
#         # data=data[data['mjd']<58500]
#         # data=data[data['mjd']<57850]
#         # fig = plt.figure()
#         # gs = gridspec.GridSpec(1,1)
#         # ax1 = plt.subplot(gs[0,0])
#         # ax1.hist(data['mjd'],bins="auto")
#         # plt.show()
#
#         # drop any measurements with zero uncertainty
#         data_zero_err=data[data['merr']==0]
#         data=data[data['merr']!=0]
#         # drop measurements with small (or zero) uncertainty MOVE THIS TO A VARIABLE UP TOP!
#         data_small_err=data[data['merr']<mag_err_small]
#         data=data[~(data['merr']<mag_err_small)]
#
#         # drop measurements near galactic plane
#         data_gal=data[np.absolute(data["galactic_latitude"])<gal_lat_cut]
#         data=data[~(np.absolute(data["galactic_latitude"])<gal_lat_cut)]
#
#         if len(data)==0:
#             print("no data, cannot fit")
#             break
#
#         # iteratively fit and cut data
#         k=0
#         while k<max_iters:
#
#             # print("iteration: {}".format(k))
#             print("iteration: {}, N_data={}".format(k,len(data)))
#
#             # extract asteroid phase data from data
#             alpha = np.array(data['phase_angle']) * u.deg
#             mag = np.array(data["reduced_mag"]) * u.mag
#             mag_err = np.array(data["merr"]) * u.mag
#
#             if k==0:
#                 # for first iteration start with the predicted HG mag
#                 model=HG(H = H_abs_mag * u.mag, G = G_slope)
#                 model_str="predicted mag HG"
#
#                 param_names=model.param_names
#                 params=model.parameters
#
#                 # label the first fit
#                 labels=[model_str]
#                 for l in range(len(model.parameters)):
#                     labels.append("{}={:.2f} ".format(param_names[l],params[l]))
#                 label=", ".join(labels)
#
#             else:
#                 # apply the fitter to the data
#                 if len(alpha)<=len(phase_curve[i]): # check that there is enough data to fit
#                     print("less data to fit than parameters")
#                     break
#                 model = fitter(model_name, alpha, mag, weights=1.0/np.array(mag_err)) # fit using weights by uncertainty
#                 # if k==1:
#                 #     model = fitter(model_name, alpha, mag) # drop the weights for the first fit
#                 # else:
#                 #     model = fitter(model_name, alpha, mag, weights=1.0/np.array(mag_err))
#
#                 model_str=model_names_str[i]
#
#                 # !!! RECORD ANY WARNINGS ETC? see fitter.fit_info
#
#                 param_names=model.param_names
#                 params=model.parameters
#                 # print(model.__dict__)
#                 # print()
#                 # print(fitter.fit_info['cov_x']) # The scipy.optimize.leastsq result for the most recent fit https://docs.astropy.org/en/stable/api/astropy.modeling.fitting.LevMarLSQFitter.html
#                 # print(err_x)
#
#                 # label each fit
#                 labels=["{}. {}".format(k,model_str)]
#                 for l in range(len(model.parameters)):
#                     labels.append("{}={:.2f} ".format(param_names[l],params[l]))
#                 label=", ".join(labels)
#
#                 # test for convergence
#                 # find difference in params
#                 delta_params = np.absolute(old_params-params)
#                 # print(old_params)
#                 # print(fitter.fit_info)
#                 # print("nfev:{}".format(fitter.fit_info['nfev']))
#                 # print("ierr:{}".format(fitter.fit_info['ierr']))
#                 # print("message:{}".format(fitter.fit_info['message']))
#                 # print("N_data:{}".format(len(data)))
#                 # exit()
#                 # print(params)
#                 # print("delta params = {}".format(delta_params))
#
#                 if np.sum(delta_params<param_converge_check)==len(delta_params):
#                     # print("converged")
#                     print(params)
#
#                     # retrieve the fit metrics
#                     x_vals=params
#                     param_cov=fitter.fit_info['param_cov'] # see notes of https://docs.astropy.org/en/stable/api/astropy.modeling.fitting.LevMarLSQFitter.html for difference between param_cov and cov_x
#                     if param_cov is None:
#                         print("A value of None indicates a singular matrix, which means the curvature in parameters x is numerically flat")
#                         # What should I do here?
#                         break
#
#                     # cov_x=fitter.fit_info['cov_x']
#                     # print("cov_x:\n{}".format(cov_x))
#                     # err_x = np.sqrt(np.diag(cov_x))
#                     # print("err_x: {}".format(err_x))
#                     # print(param_cov/cov_x)
#
#                     param_err_x = np.sqrt(np.diag(param_cov)) # retrieve errors in the parameters: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
#                     param_shape=param_cov.shape
#                     correl_coef=np.zeros(param_shape)
#                     for l in range(param_shape[0]):
#                         for m in range(param_shape[1]):
#                             correl_coef[l,m]=param_cov[l,m]/np.sqrt(param_cov[l,l]*param_cov[m,m]) # Hughes AND Hase eqn. 7.3
#                     N_data_fit=len(data)
#                     # nfev=fitter.fit_info['nfev']
#                     # print(np.array(data["mjd"]).astype(int))
#                     N_nights=len(np.unique(np.array(data["mjd"]).astype(int)))
#                     alpha_min=np.amin(alpha).value
#                     alpha_max=np.amax(alpha).value
#                     phase_angle_range=alpha_max-alpha_min
#                     # N_alpha_low=sum(alpha<low_alpha_cut)
#                     N_alpha_low=len(alpha[alpha<low_alpha_cut])
#                     N_iter=k
#                     nfev=fitter.fit_info['nfev']
#                     ier=fitter.fit_info['ierr']
#                     N_mag_err=len(mag_err[np.array(mag_err)<mag_err_threshold])
#
#                     # # time len() vs sum()
#                     # import time
#                     # start = time.process_time()
#                     # sum(alpha<low_alpha_cut)
#                     # sum(np.array(mag_err)<mag_err_threshold)
#                     # print(time.process_time() - start)
#                     #
#                     # start = time.process_time()
#                     # len(alpha[alpha<low_alpha_cut])
#                     # len(mag_err[np.array(mag_err)<mag_err_threshold])
#                     # print(time.process_time() - start)
#
#                     # exit()
#
#                     # print(fitter.fit_info)
#                     print(fitter.fit_info['message'])
#
#                     if push_fit==True:
#
#                         print("push fit to db")
#
#                         # print(params)
#                         # print(np.isnan(params))
#                         # print(param_err_x)
#                         # print(np.isnan(param_err_x))
#                         # for prms in [params, param_err_x]:
#                         #     prms=list(prms)
#                         # #     mask=np.isnan(prms)
#                         # #     prms[mask]="NaN"
#                         # #     print(prms)
#                         # param_err_x[0]=None
#                         # print(params)
#                         # print(param_err_x)
#
#
#                         HG_params_str=""
#                         for p in range(len(phase_curve[i])):
#                             HG_params_str+="phase_curve_{}%(ms)s_%(filt)s={},phase_curve_{}_err%(ms)s_%(filt)s={},".format(
#                             phase_curve[i][p],params[p],
#                             phase_curve[i][p],param_err_x[p]
#                             ) % locals()
#
#                         qry = u"""UPDATE %(tab_name)s SET
#                         detection_count_%(filt)s=%(detection_count_filt)s,
#                         %(HG_params_str)s
#                         phase_curve_N_fit%(ms)s_%(filt)s=%(N_data_fit)s,
#                         phase_curve_alpha_min%(ms)s_%(filt)s=%(alpha_min)s,
#                         phase_curve_alpha_max%(ms)s_%(filt)s=%(alpha_max)s,
#                         phase_angle_range_%(filt)s=%(phase_angle_range)s,
#                         phase_curve_N_alpha_low%(ms)s_%(filt)s=%(N_alpha_low)s,
#                         phase_curve_N_nights%(ms)s_%(filt)s=%(N_nights)s,
#                         phase_curve_N_iter%(ms)s_%(filt)s=%(N_iter)s,
#                         phase_curve_refresh_date_%(filt)s='%(utc_date_now)s',
#                         phase_curve_nfev%(ms)s_%(filt)s=%(nfev)s,
#                         phase_curve_ier%(ms)s_%(filt)s=%(ier)s,
#                         phase_curve_N_mag_err%(ms)s_%(filt)s=%(N_mag_err)s
#                         WHERE mpc_number=%(mpc_number)s;""" % locals()
#
#                         qry=qry.replace('=nan', '=NULL') # mysql doesn't understand nan, needs NULL
#
#                         print(qry)
#                         # exit()
#
#                         cursor2.execute(qry)
#                         cnx2.commit()
#
#                     if plot_fig:
#
#                         import matplotlib.pyplot as plt
#                         import matplotlib.gridspec as gridspec
#
#                         fig = plt.figure()
#                         # gs = gridspec.GridSpec(1,1)
#                         # ax1 = plt.subplot(gs[0,0])
#                         gs = gridspec.GridSpec(3,1,height_ratios=[1,1,0.1])
#                         ax1 = plt.subplot(gs[0,0])
#                         ax2 = plt.subplot(gs[1,0])
#                         ax3 = plt.subplot(gs[2,0])
#
#                         # plot all the data from the SQL that goes into the fitting process
#                         ax1.errorbar(data_filt['phase_angle'],data_filt['reduced_mag'],data_filt['merr'], fmt='ko',label="data",zorder=0,markersize="2")
#
#                         # highlight any measurements with zero uncertainty
#                         ax1.scatter(data_zero_err['phase_angle'],data_zero_err['reduced_mag'],c='r',marker="+",s=50)
#                         ax1.scatter(data_small_err['phase_angle'],data_small_err['reduced_mag'],edgecolor='r',facecolor="none",marker="o",s=50)
#                         # highlight low galactic latitude
#                         ax1.scatter(data_gal['phase_angle'],data_gal['reduced_mag'],c='r',marker="x",s=50)
#
#                         # plot iterative fits and cuts
#                         alpha_fit=np.linspace(np.amin(alpha),np.amax(alpha),100)
#                         print(label_iter_list)
#                         print(model_iter_list)
#                         print(k)
#                         for j in range(len(model_iter_list)):
#                             print(j,label_iter_list[j])
#                             ax1.plot(alpha_fit,model_iter_list[j](alpha_fit),label=label_iter_list[j])
#                             ax1.scatter(alpha_cut_iter_list[j],mag_cut_iter_list[j],marker="x",zorder=3)
#
#                         # plot the phase data after the cuts
#                         ax2.errorbar(np.array(alpha),np.array(mag),np.array(mag_err),fmt='k.',alpha=0.2,label="data error",zorder=0,markersize="2")
#                         s2=ax2.scatter(np.array(alpha),np.array(mag),c=np.array(data["mjd"]),label="data time",s=10)
#                         cbar2=fig.colorbar(s2,ax3,use_gridspec=True, orientation='horizontal')
#
#                         # plot the final fit
#                         ax2.plot(alpha_fit,model(alpha_fit),label=label)
#
#                         ax1.set_xlabel('alpha(degrees)')
#                         ax1.set_ylabel('mag')
#                         ax1.invert_yaxis()
#                         ax1.legend(prop={'size': 6})
#
#                         ax2.set_xlabel('alpha(degrees)')
#                         ax2.set_ylabel('mag')
#                         ax2.invert_yaxis()
#                         ax2.legend(prop={'size': 6})
#
#                         # fig.set_size_inches(4,10)
#
#                         ax3.set_xlabel("MJD")
#
#                         ax1.set_title("{}_{}_{}_{}_{}".format(os.path.basename(__file__).split('.')[0],obj_number,model_names_str[i],clip_label,filt))
#                         plt.tight_layout()
#
#                         if save_fig:
#                             fname="{}/{}_{}_{}_{}_{}.png".format(save_path,os.path.basename(__file__).split('.')[0],obj_number,model_names_str[i],clip_label,filt)
#                             print(fname)
#                             plt.savefig(fname, bbox_inches='tight')
#
#                         if show_fig:
#                             plt.show()
#                         else:
#                             plt.close()
#
#                     # exit()
#                     break
#
#             # record stuff for plotting
#             label_iter_list.append(label)
#             model_iter_list.append(model)
#
#             # cut the outlying data
#             mask=data_clip_func(mag, model(alpha),std)
#
#             # record the cut data points for plotting
#             mag_cut_iter_list.append(mag[mask])
#             alpha_cut_iter_list.append(alpha[mask])
#
#             # define the data to keep
#             data=data[~mask]
#
#             # store these params as old params, to check for convergence of the next fit
#             if k>0:
#                 old_params=params
#
#             k+=1
#
# cnx1.close()
# cnx2.close()
