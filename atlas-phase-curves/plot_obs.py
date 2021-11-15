import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u
from calculate_phase.atlas_SQL_query_df import get_orb_elements_id,atlas_SQL_query_orbid,atlas_SQL_query_orbid_expname
from calculate_phase import atlas_database_connection
from astropy.time import Time
from datetime import datetime
import os

plot_option=1
orbfit_separation_arcsec_cut=1.0
size = 30
save_path = "/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/results_analysis/orbfit_figs"
save_file_type="png"

mjd_lim = 59259 # 2021-02-14 00:00:00.000
today = Time(datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")).mjd

mpc_number = False
# name = "2010 ET56"
name = "Osipovia"
# name = "Fitzsimmons"
# name = "Bettina"
# name = "Moguntia"

if name:
    obj_title=name
else:
    obj_title=mpc_number

# connect to database
cnx=atlas_database_connection.database_connection().connect()

# get data
orbid=get_orb_elements_id(cnx,mpc_number,name)
data_all_filt=atlas_SQL_query_orbid_expname(cnx,orbid)
data = data_all_filt
print(data)

data = data.sort_values("orbfit_separation_arcsec")

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

    s1=ax1.scatter(alpha,mag,c=time,s=size,label="data time")
    ax1.errorbar(np.array(alpha),np.array(mag),np.array(mag_err),fmt='k.',alpha=0.2,label="data error",zorder=0,markersize="2")
    fig.colorbar(s1, ax=ax1,label="MJD")

    s2=ax2.scatter(time,mag,c=data["orbfit_separation_arcsec"],s=size,cmap="plasma")
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

    ax1.set_xlabel('phase angle')
    ax1.set_ylabel('reduced mag')
    ax2.set_xlabel('time(mjd)')
    ax2.set_ylabel('reduced mag')

    ax1.invert_yaxis()
    ax2.invert_yaxis()

    # ax1.legend()#prop={'size': 6})

    fig.suptitle(obj_title)
    plt.tight_layout()
    plt.show()

if plot_option==1:
    fig = plt.figure()
    gs = gridspec.GridSpec(2,2)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[0,1])
    ax3 = plt.subplot(gs[1,1])
    ax4 = plt.subplot(gs[1,0])

    ax1.hist(data["orbfit_separation_arcsec"],bins="auto")
    ax1.axvline(orbfit_separation_arcsec_cut,c="r")

    s2=ax2.scatter(alpha,mag,c=data["orbfit_separation_arcsec"],s=size,label="data time")
    ax2.errorbar(np.array(alpha),np.array(mag),np.array(mag_err),fmt='k.',alpha=0.2,label="data error",zorder=0,markersize="2")
    fig.colorbar(s2, ax=ax2,label="orbfit_separation_arcsec")

    mask = data["orbfit_separation_arcsec"]<orbfit_separation_arcsec_cut
    print(len(alpha[mask])/len(alpha))
    # ax2.scatter(alpha[mask],mag[mask],edgecolor="r",facecolor="none",s=size)
    ax2.scatter(alpha[~mask],mag[~mask],edgecolor="r",facecolor="none",s=size)

    s3=ax3.scatter(time,mag,c=data["orbfit_separation_arcsec"],s=size,label="data time")
    ax3.errorbar(np.array(time),np.array(mag),np.array(mag_err),fmt='k.',alpha=0.2,label="data error",zorder=0,markersize="2")
    fig.colorbar(s3, ax=ax3,label="orbfit_separation_arcsec")
    # ax3.scatter(time[mask],mag[mask],edgecolor="r",facecolor="none",s=size)
    ax3.scatter(time[~mask],mag[~mask],edgecolor="r",facecolor="none",s=size)

    ax4.scatter(time,data["orbfit_separation_arcsec"],s=size)
    ax4.scatter(time[~mask],data["orbfit_separation_arcsec"][~mask],edgecolor="r",facecolor="none",s=size)

    ax2.invert_yaxis()
    ax3.invert_yaxis()
    ax1.set_xlabel("orbfit_separation_arcsec")
    ax2.set_xlabel("phase angle")
    ax3.set_xlabel("time(mjd)")
    ax4.set_xlabel("time(mjd)")
    ax1.set_ylabel("number of detections")
    ax2.set_ylabel("reduced mag")
    ax3.set_ylabel("reduced mag")
    ax4.set_ylabel("orbfit_separation_arcsec")

    fig.suptitle(obj_title)

    plt.tight_layout()

    # fname="{}/{}_{}.{}".format(save_path,os.path.basename(__file__).split('.')[0],"_".join(obj_title.split(" ")),save_file_type)
    # print(fname)
    # plt.savefig(fname, bbox_inches='tight')

    plt.show()
