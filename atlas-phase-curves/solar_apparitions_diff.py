"""
Calls the atlas_SQL_query function to get an object's observations from ATLAS
Uses the solar elongation to analyse the apparitions of an object.
Fits a saw tooth function to the solar elongation, when the wave drops at 360 deg an apparition ends
"""

import pandas as pd
from calculate_phase import atlas_SQL_query_df
from calculate_phase import atlas_database_connection
from optparse import OptionParser
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.time import Time
from astroquery.jplhorizons import Horizons

parser = OptionParser()
parser.add_option( "-n", "--mpc-number", dest="mpc_number", help="mpc_number", metavar="MPC_NUMBER" ) # mpc number of object to fit
parser.add_option( "-N", "--name", dest="name", help="name", metavar="NAME" ) # NAME of object to fit
parser.add_option( "-f", "--filter", dest="filter", help="filter", metavar="FILTER" ) # colour filter of data
(options,args)=parser.parse_args()

if options.mpc_number:
    mpc_number=int(options.mpc_number)
else:
    mpc_number=False
if options.name:
    name=options.name
else:
    name=False
if options.filter:
    filter=options.filter
else:
    filter=False

elong_mod = 360
sol_elong_diff = -10 # difference required for an epoch, should be minimised! (e.g. 76820)
mjd_jd = 2400000.5 # conversion from MJD to JD

# connect to database
cnx=atlas_database_connection.database_connection().connect()

# get orbit_id
orbid=atlas_SQL_query_df.get_orb_elements_id(cnx,mpc_number,name)

# do the query
df_data=atlas_SQL_query_df.atlas_SQL_query_orbid_expname(cnx,orbid)

print(df_data)
print(list(df_data))

# drop any nans
df_data = df_data.dropna()

# select only a particular filter
if filter:
    df_data=df_data[df_data["filter"]==filter]

# save df
# df_data.to_csv("results_analysis/obs/df_data_{}.csv".format(mpc_number))

# count extreme low error measurements
err_cut=0.005
N_merr = len(df_data[df_data["merr"]<err_cut])
print("N merr<{} : {}".format(err_cut,N_merr))

fig = plt.figure()
gs = gridspec.GridSpec(2, 1)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])

# plot the phase angle/reduced mag data
ax1.invert_yaxis()
ax1.errorbar(df_data['phase_angle'],df_data['reduced_mag'],df_data['merr'], fmt='ko',label="mag err",zorder=0,markersize=0)
ax1.set_xlabel("phase_angle")
ax1.set_ylabel("reduced_mag")

# plot reduced mag as a function of mjd date
# ax2.scatter(df_data['mjd'],df_data['reduced_mag'])
ax2.set_xlabel("mjd")
ax2.set_ylabel("reduced_mag")
ax2.invert_yaxis()

# plot the solar elongation retrieved from rockAtlas
# N.B. rockAtlas angle cycles from -180 to 180 deg, use mod to get 0 to 360
# or abs to force 0 -> 180

df_data=df_data.sort_values(['mjd'])
xdata=np.array(df_data["mjd"])

# ydata=np.array(df_data["sun_obs_target_angle"])
# lab="solar elongation"

ydata=np.array(df_data["sun_obs_target_angle"])%elong_mod
lab="solar elongation % 360"

# ydata=np.absolute(np.array(df_data["sun_obs_target_angle"]))
# lab="abs(solar elongation)"

ax3 = ax2.twinx()
s=3
ax3.scatter(xdata,ydata,s=s,c="k",label=lab)
ax3.scatter(df_data['mjd'],df_data['phase_angle'],s=s,c="gray",label="phase angle")
ax3.set_ylabel("angle")

# Get the JPL Horizons elongation
if name:
    obj_id = name
else:
    obj_id = mpc_number
loc = "T05"
t1=Time(int(np.amin(df_data["mjd"])),format="mjd")
t2=Time(int(np.amax(df_data["mjd"])),format="mjd")
epoch_list = {'start':t1.iso, 'stop':t2.iso, 'step':'1d'} # a range of epochs in Horizons format is FAST!
obj = Horizons(id=obj_id, location=loc, epochs=epoch_list)
eph = obj.ephemerides()
df_eph = eph.to_pandas()
df_eph["mjd"] = df_eph["datetime_jd"] - mjd_jd
# plot the exact elongation from JPL
ax3.plot(df_eph["mjd"],df_eph["elong"],c="r",alpha=0.2,label = "JPL elong ({})".format(loc),zorder=0)
mask = (df_eph["elong"]<10.0)
ax3.scatter(df_eph[mask]["mjd"],df_eph[mask]["elong"],c="r",s=1,label = "JPL elong<10 deg",zorder=0)


# find the turning points of the solar elongation data
ydiff = np.diff(ydata)

# find any -ve drop
# turning_points = xdata[1:][ydiff<0]

# select only drops of at least a certain size
turning_points = xdata[1:][ydiff<sol_elong_diff]

# add the first and last mjds in the dataset to fully bound all apparitions
turning_points = np.insert(turning_points,0,xdata[0])
turning_points = np.insert(turning_points,len(turning_points),xdata[-1])
print(turning_points)
for i,t in enumerate(turning_points):
    if i==0 or i==len(turning_points)-1:
        ax3.axvline(t,c="r",ls=":",alpha=0.3)
    else:
        ax3.axvline(t,c="r")

# use these turning points to divide into epochs
for i in range(1,len(turning_points)):
    t1=turning_points[i-1]
    t2=turning_points[i]

    df = df_data[(df_data["mjd"]>=t1) & (df_data["mjd"]<t2)]

    # ax2.scatter(df['mjd'],df['reduced_mag'],s=1,c="C{}".format(i))
    ax2.scatter(df['mjd'],df['reduced_mag'],s=2)
    ax1.scatter(df['phase_angle'],df['reduced_mag'],s=2,label="epoch:{}".format(i))

ax1.legend()
ax3.legend()

plt.tight_layout()

plt.show()
