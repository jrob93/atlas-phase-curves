import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u
from calculate_phase.atlas_SQL_query_df import get_orb_elements_id,atlas_SQL_query_orbid,atlas_SQL_query_orbid_expname
from calculate_phase import atlas_database_connection
from astropy.time import Time
from datetime import datetime

plot_option=1
orbfit_separation_arcsec_cut =0.5

mjd_lim = 59259 # 2021-02-14 00:00:00.000
today = Time(datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")).mjd

mpc_number = False
# name = "2010 ET56"
name = "Osipovia"

# connect to database
cnx=atlas_database_connection.database_connection().connect()

# get data
orbid=get_orb_elements_id(cnx,mpc_number,name)
data_all_filt=atlas_SQL_query_orbid_expname(cnx,orbid)
data = data_all_filt
print(data)

time=np.array(data["mjd"])
alpha = np.array(data['phase_angle'])# * u.deg
mag = np.array(data["reduced_mag"])# * u.mag
mag_err = np.array(data["merr"])# * u.mag

if plot_option==0:
    fig = plt.figure()
    gs = gridspec.GridSpec(2,1)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0])

    # fig.set_size_inches(20,12)

    s1=ax1.scatter(alpha,mag,c=time,s=10,label="data time")
    ax1.errorbar(np.array(alpha),np.array(mag),np.array(mag_err),fmt='k.',alpha=0.2,label="data error",zorder=0,markersize="2")
    fig.colorbar(s1, ax=ax1,label="MJD")

    s2=ax2.scatter(time,mag,c=data["orbfit_separation_arcsec"],s=10,cmap="plasma")
    fig.colorbar(s2, ax=ax2,label="orbfit_separation_arcsec")
    ax2.axvline(mjd_lim)
    ax2.axvline(today)

    # # cut mag
    # mag_med_cut=5
    # mag_med = np.nanmedian(mag)
    # mask=np.absolute(np.array(mag-mag_med))>mag_med_cut
    # ax1.axhline(mag_med)
    # ax1.axhline(mag_med-mag_med_cut)
    # ax1.scatter(alpha[mask],mag[mask],c='r',marker="x")

    ax1.set_xlabel('alpha(degrees)')
    ax1.set_ylabel('mag')
    ax1.invert_yaxis()

    ax1.legend()#prop={'size': 6})

    plt.show()

if plot_option==1:
    fig = plt.figure()
    gs = gridspec.GridSpec(2,1)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0])

    ax1.hist(data["orbfit_separation_arcsec"],bins="auto")
    ax1.axvline(orbfit_separation_arcsec_cut,c="r")

    s2=ax2.scatter(alpha,mag,c=data["orbfit_separation_arcsec"],s=10,label="data time")
    ax2.errorbar(np.array(alpha),np.array(mag),np.array(mag_err),fmt='k.',alpha=0.2,label="data error",zorder=0,markersize="2")
    fig.colorbar(s2, ax=ax2,label="MJD")

    mask = data["orbfit_separation_arcsec"]<orbfit_separation_arcsec_cut
    print(len(alpha[mask])/len(alpha))
    ax2.scatter(alpha[mask],mag[mask],edgecolor="r",facecolor="none",s=10)
    ax2.invert_yaxis()

    plt.show()
