import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u
from calculate_phase.atlas_SQL_query_df import get_orb_elements_id,atlas_SQL_query_orbid
from calculate_phase import atlas_database_connection
from astropy.time import Time
from datetime import datetime

mjd_lim = 59259 # 2021-02-14 00:00:00.000
today = Time(datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")).mjd

mpc_number = False
name = "Pandion"

# connect to database
cnx=atlas_database_connection.database_connection().connect()

# get data
orbid=get_orb_elements_id(cnx,mpc_number,name)
data_all_filt=atlas_SQL_query_orbid(cnx,orbid)
data = data_all_filt
print(data)

time=np.array(data["mjd"])
alpha = np.array(data['phase_angle'])# * u.deg
mag = np.array(data["reduced_mag"])# * u.mag
mag_err = np.array(data["merr"])# * u.mag

fig = plt.figure()
gs = gridspec.GridSpec(2,1)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])

# fig.set_size_inches(20,12)

s1=ax1.scatter(alpha,mag,c=time,s=10,label="data time")
ax1.errorbar(np.array(alpha),np.array(mag),np.array(mag_err),fmt='k.',alpha=0.2,label="data error",zorder=0,markersize="2")
fig.colorbar(s1, ax=ax1,label="MJD")

ax2.scatter(time,mag,s=10)
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
