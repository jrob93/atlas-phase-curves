import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u
from calculate_phase.atlas_SQL_query_df import atlas_SQL_query,get_orb_elements_id,get_unique_ids,atlas_SQL_query_orbid

# data = pd.read_csv("results_analysis/obs/df_data_168576.csv",index_col=0)
# data = pd.read_csv("results_analysis/obs/df_data_150621.csv",index_col=0)
# data = pd.read_csv("results_analysis/obs/df_data_141570.csv",index_col=0)
# data = pd.read_csv("results_analysis/obs/df_data_134649.csv",index_col=0)
data = pd.read_csv("results_analysis/obs/df_data_17707.csv",index_col=0)
print(data)
# exit()

time=np.array(data["mjd"])
alpha = np.array(data['phase_angle'])# * u.deg
mag = np.array(data["reduced_mag"])# * u.mag
mag_err = np.array(data["merr"])# * u.mag

mag_med_cut=5
mag_med = np.nanmedian(mag)
# mag=mag-mag_med
mask=np.absolute(np.array(mag-mag_med))>mag_med_cut
# print(mask)

fig = plt.figure()
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])

fig.set_size_inches(20,12)

s1=ax1.scatter(alpha,mag,c=time,s=10,label="data time")
ax1.errorbar(np.array(alpha),np.array(mag),np.array(mag_err),fmt='k.',alpha=0.2,label="data error",zorder=0,markersize="2")
fig.colorbar(s1, ax=ax1,label="MJD")

ax1.axhline(mag_med)
ax1.axhline(mag_med-mag_med_cut)
ax1.scatter(alpha[mask],mag[mask],c='r',marker="x")

ax1.set_xlabel('alpha(degrees)')
ax1.set_ylabel('mag')
ax1.invert_yaxis()

ax1.legend()#prop={'size': 6})

plt.show()
