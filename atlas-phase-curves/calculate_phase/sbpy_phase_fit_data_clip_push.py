# https://sbpy.readthedocs.io/en/latest/sbpy/photometry.html#disk-integrated-phase-function-models

import numpy as np
import astropy.units as u
from astropy.modeling.fitting import LevMarLSQFitter
from sbpy.photometry import HG, HG1G2, HG12, HG12_Pen16
import pandas as pd
import os
import mysql.connector
from optparse import OptionParser
import subprocess
import sys
from atlas_SQL_query_df import atlas_SQL_query
import datetime

# !!! FUNCTIONALISE THIS?
parser = OptionParser()
parser.add_option( "-n", "--mpc-number", dest="mpc_number", help="mpc_number", metavar="MPC_NUMBER" )
# parser.add_option( "-o", "--obj-name", dest="obj_name", help="obj_name", metavar="OBJ_NUMBER" )
# parser.add_option( "-s", "--save-path", dest="save_path", help="save_path", metavar="SAVE_PATH" )

(options,args)=parser.parse_args()

if options.mpc_number:
    mpc_number=options.mpc_number
    # print("mpc_number:",mpc_number)
else:
    mpc_number="4986"

# define object
# LOAD OR RUN QUERY? ADD FILTER INFO (o or c?)
obj_number=mpc_number
print(mpc_number)
low_alpha_cut=5.0*u.deg # we want to quantify how many data points are fit at low phase angles, alpha < low_alpha_cut
param_converge_check=0.01 # the model is fit until the change in parameters (e.g. H and G) is less than param_converge_check (or max iterations is reached)
max_iters=30 # maximum number of attempts at fitting
std=2 # standard deviation of the sigma data clip
mag_err_threshold = 0.1 # limit for the error of "good" data, we record N_mag_err number of data points with error < mag_err_threshold

utc_date_now=datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S") # time at which the script is run (UTC)

def data_clip_sigma(data,data_predict,low=2,high=2):
    # cut outliers by sigma
    std=np.std(data)
    clip_mask=((data < (data_predict-(std*low))) | (data > (data_predict+(std*high))))
    return clip_mask

# def data_clip_diff(data,data_predict,diff=1):
#     # cut outliers by diff (this function doesn't like astropy units, use np arrays)
#     diff=1.0
#     clip_mask=(np.absolute(data_predict-data)>diff)
#     return clip_mask

data_clip_func=data_clip_sigma

# fit new models to model data
fitter = LevMarLSQFitter()
model_names_str = ["HG", "HG1G2", "HG12", "HG12_Pen16"]
filters = ["o","c"]
model_short = ["_B89","_3M10","_2M10","_P16"] # use shorthand for all models
phase_curve = [["H","G"],["H","G1","G2"],["H","G12"],["H","G12"]]

model_names = [HG(), HG1G2(), HG12(), HG12_Pen16()]

# Predicted mag using parameters from orbital_elements table (USE ONLY FOR THE 1-MAG DIFF CUT?)
# connect to the database (must be connected to dormammu: sshpass -e ssh dormammu) to read from
config1 = {
  'user': 'af',
  'password': 'afPass',
  'host': 'localhost',
  'port': '3308',
  'database': 'atlas_moving_objects',
  'raise_on_warnings': True
}
cnx1 = mysql.connector.connect(**config1)
cursor1 = cnx1.cursor()

# connect to the database to write to
# !!! too many connections! use the same connection but switch db in qry?
config2 = {
  'user': 'af',
  'password': 'afPass',
  'host': 'localhost',
  'port': '3308',
  'database': 'jamie_test_db',
  'raise_on_warnings': True
}
cnx2 = mysql.connector.connect(**config2)
cursor2 =cnx2.cursor()
# !!! merge conections when we start writing to the same db, try:
# cnx2 = cnx1
# cursor2 = cursor1

tab_name="test_table"

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
FROM atlas_objects WHERE mpc_number=%(mpc_number)s
""" % locals()
# print(qry_obj1)

df_obj=pd.read_sql_query(qry_obj1,cnx1)
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

qry_HG="select G_slope, H_abs_mag from orbital_elements where primaryId='{}';".format(obj_number)
df_HG=pd.read_sql_query(qry_HG,cnx1)

cursor2.execute("SELECT COUNT(1) FROM {} where mpc_number={}".format(tab_name,mpc_number))
mpc_check=cursor2.fetchone()[0]
# print("mpc_check = {}".format(mpc_check))

# load data to be fitted, loads both filters (o & c)
data_all_filt=atlas_SQL_query(cnx=cnx1,mpc_number=obj_number)
detection_count=len(data_all_filt) # note that atlas_objects may not have up to date detection count...
# print(data_all_filt)

if mpc_check==0:
    # qry_obj = u"""INSERT INTO {}
    # (mpc_number,
    # last_detection_mjd,
    # name,
    # orbital_elements_id,
    # primaryId)
    # VALUES
    # (%(mpc_number)s,{},'{}',{},{});""".format(tab_name,
    # float(df_obj['last_detection_mjd']),
    # str(df_obj['name'].iloc[0]),
    # int(df_obj['orbital_elements_id']),
    # int(df_obj['primaryId'])) % locals()

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
    (%(mpc_number)s,{},{},{},{},{},'{}',{},{});""".format(tab_name,
    str(df_obj['dateLastModified'].iloc[0]),
    detection_count,
    float(df_obj['last_detection_mjd']),
    str(df_obj['last_photometry_update_date_c'].iloc[0]),
    str(df_obj['last_photometry_update_date_o'].iloc[0]),
    str(df_obj['name'].iloc[0]),
    int(df_obj['orbital_elements_id']),
    int(df_obj['primaryId'])) % locals()

    # !!! insert detection count etc here?

else:
    qry_obj = u"""UPDATE {} SET
    dateLastModified={},
    detection_count={},
    last_detection_mjd={},
    last_photometry_update_date_c={},
    last_photometry_update_date_o={},
    name='{}',
    orbital_elements_id={},
    primaryId={}
    WHERE mpc_number=%(mpc_number)s;""".format(tab_name,
    str(df_obj['dateLastModified'].iloc[0]),
    detection_count,
    float(df_obj['last_detection_mjd']),
    str(df_obj['last_photometry_update_date_c'].iloc[0]),
    str(df_obj['last_photometry_update_date_o'].iloc[0]),
    str(df_obj['name'].iloc[0]),
    int(df_obj['orbital_elements_id']),
    int(df_obj['primaryId'])) % locals()

# !!! detection count is NaN a lot...? CHECK METRICS
# !!! use detection_count to track if object should be refit? no new data, no refit...
# !!! or use last photometry update date

# print(qry_obj)

# !!! WHAT TO DO WITH updated FLAG? use it to track number of times the fit has been done? Change to a date?

cursor2.execute(qry_obj)
cnx2.commit()

for filt in filters:

    # H and G values for the predicted fit
    G_slope=float(df_HG.iloc[0]['G_slope'])
    H_abs_mag=float(df_HG.iloc[0]['H_abs_mag'])
    # print("G_slope = {}\nH_abs_mag = {}".format(G_slope,H_abs_mag))

    # !!! SHIFT HG TO o AND c FILTERS!

    # select all data from a certain filter
    data_filt=data_all_filt[data_all_filt['filter']==filt]
    # print(data_filt)
    detection_count_filt=len(data_filt) # !!! SHOULD WE RECORD NUMBER OF DETECTIONS BEFORE/AFTER CUTS?
    # print(detection_count_filt)
    # continue

    # iterate over all models
    for i,model_name in enumerate(model_names):

        ms=model_short[i]

        # if ms!=model_short[1]:
        #     continue

        old_params=[999]*len(model_name.parameters)

        # iterate over data cuts?
        # !!! SELECT ONE AND SET IT: 2-sigma

        print("fit {}, filter {}".format(model_names_str[i],filt))
        # print("fit {}, {}".format(model_names_str[i],j))

        # initialise the data that we will iteratively fit and cut
        data=data_filt
        data=data.sort_values("phase_angle") # ensure that the dataframe is in order for plotting
        data_zero_err=data[data['merr']==0]
        data=data[data['merr']!=0] # drop any measurements with zero uncertainty

        # print(data)

        if len(data)==0:
            print("no data, cannot fit")
            break

        # plot the asteroid phase data
        # extract asteroid phase data from data
        alpha = np.array(data['phase_angle']) * u.deg
        mag = np.array(data["reduced_mag"]) * u.mag
        mag_err = np.array(data["merr"]) * u.mag

        # iteratively fit and cut data
        k=0
        while k<max_iters:

            print("iteration: {}".format(k))

            # extract asteroid phase data from data
            alpha = np.array(data['phase_angle']) * u.deg
            mag = np.array(data["reduced_mag"]) * u.mag
            mag_err = np.array(data["merr"]) * u.mag
            # try:
            #     alpha_fit=np.linspace(np.amin(alpha),np.amax(alpha),100)
            # except:
            #     break
            # print(mag_err)

            if k==0:
                # for first iteration start with the predicted HG mag
                model=HG(H = H_abs_mag * u.mag, G = G_slope)
                model_str="predicted mag HG"

                param_names=model.param_names
                params=model.parameters
                # labels=[model_str]
                # for l in range(len(model.parameters)):
                #     labels.append("{}={:.2f} ".format(param_names[l],params[l]))
                # label=", ".join(labels)

            else:
                # apply the fitter to the data
                # print(mag_err)
                # print(1.0/np.array(mag_err))
                # exit()

                if len(alpha)<=len(phase_curve[i]):
                    print("less data to fit than parameters")
                    break

                model = fitter(model_name, alpha, mag, weights=1.0/np.array(mag_err)) # add weights/uncertainties here?
                model_str=model_names_str[i]

                # !!! RECORD ANY WARNINGS ETC? see fitter.fit_info

                param_names=model.param_names
                params=model.parameters
                # print(model.__dict__)
                # print()
                # print(fitter.fit_info['cov_x']) # The scipy.optimize.leastsq result for the most recent fit https://docs.astropy.org/en/stable/api/astropy.modeling.fitting.LevMarLSQFitter.html
                # print(err_x)

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
                if np.sum(delta_params<param_converge_check)==len(delta_params):
                    # print("converged")
                    print(params)

                    # retrieve the fit metrics
                    x_vals=params
                    param_cov=fitter.fit_info['param_cov'] # see notes of https://docs.astropy.org/en/stable/api/astropy.modeling.fitting.LevMarLSQFitter.html for difference between param_cov and cov_x
                    if param_cov is None:
                        print("A value of None indicates a singular matrix, which means the curvature in parameters x is numerically flat")
                        # What should I do here?
                        break

                    # cov_x=fitter.fit_info['cov_x']
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
                    # nfev=fitter.fit_info['nfev']
                    # print(np.array(data["mjd"]).astype(int))
                    N_nights=len(np.unique(np.array(data["mjd"]).astype(int)))
                    alpha_min=np.amin(alpha).value
                    alpha_max=np.amax(alpha).value
                    phase_angle_range=alpha_max-alpha_min
                    # N_alpha_low=sum(alpha<low_alpha_cut)
                    N_alpha_low=len(alpha[alpha<low_alpha_cut])
                    N_iter=k
                    nfev=fitter.fit_info['nfev']
                    ier=fitter.fit_info['ierr']
                    N_mag_err=len(mag_err[np.array(mag_err)<mag_err_threshold])

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

                    # print(fitter.fit_info)
                    print(fitter.fit_info['message'])

                    HG_params_str=""
                    for k in range(len(phase_curve[i])):
                        HG_params_str+="phase_curve_{}%(ms)s_%(filt)s={},phase_curve_{}_err%(ms)s_%(filt)s={},".format(
                        phase_curve[i][k],params[k],
                        phase_curve[i][k],param_err_x[k]
                        ) % locals()

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

                    # print(qry)
                    # exit()

                    cursor2.execute(qry)
                    cnx2.commit()

                    break

            # cut and plot the outlying data
            # !!! remove options for speed
            # if j==0:
            # std=2
            mask=data_clip_func(mag, model(alpha),std)
            # clip_label="{}-sigma_clip".format(std)
            # elif j==1:
            #     std=3
            #     mask=data_clip_sigma(mag, model(alpha),std)
            #     clip_label="{}-sigma_clip".format(std)
            # else:
            #     diff=1
            #     mask=data_clip_diff(np.array(mag), np.array(model(alpha)),diff)
            #     clip_label="{}-mag_diff_clip".format(diff)

            mag_cut=mag[mask]
            alpha_cut=alpha[mask]

            # define the data to keep
            data=data[~mask]

            if k>0:
                old_params=params

            k+=1

cnx1.close()
cnx2.close()
