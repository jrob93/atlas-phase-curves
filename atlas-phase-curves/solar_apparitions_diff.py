"""
Calls the atlas_SQL_query function to get an object's observations from ATLAS
Uses the solar elongation to analyse the apparitions of an object.
Looks for drops in solar elongation angle to define epochs.

N.B. will save queried data and try reload, these files may need updated over time!
Delete the loaded files to get fresh ones
"""

# non-interactive
import matplotlib
matplotlib.use('Agg')

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
import os

parser = OptionParser()
parser.add_option( "-n", "--mpc-number", dest="mpc_number", help="mpc_number", metavar="MPC_NUMBER" ) # mpc number of object to fit
parser.add_option( "-N", "--name", dest="name", help="name", metavar="NAME" ) # NAME of object to fit
parser.add_option( "-f", "--filter", dest="filter", help="filter", metavar="FILTER" ) # colour filter of data
parser.add_option( "-R", "--reload", dest="reload", help="reload", metavar="RELOAD" ) # delete existing files and reload
parser.add_option( "-s", "--save", dest="save", help="save", metavar="SAVE" ) # delete existing files and reload
parser.add_option( "-d", "--diff", dest="diff", help="diff", metavar="diff" ) # delete existing files and reload
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
if options.reload:
    reload=options.reload
else:
    reload=False
if options.save:
    save_path=options.save
else:
    save_path="solar_apparitions_figs"
if options.diff:
    sol_elong_diff = -float(options.diff)
else:
    sol_elong_diff = -10.0 # difference required for an epoch, should be minimised! (e.g. 76820)

elong_mod = 360
date_diff = 365.0 # search also for gaps of at least 100 days
P_est = 365.0 # estimated period for apparitions
mjd_jd = 2400000.5 # conversion from MJD to JD
loc = "T05" # which telescope to use for Horizons query, T05 or T08?

data_load_path = "/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/results_analysis/obs"
eph_load_path = "solar_apparitions_files"

if name:
    obj_id = name
    obj_id_save = "_".join(name.split())
else:
    obj_id = mpc_number
    obj_id_save = mpc_number

data_fname = "{}/df_data_{}.csv".format(data_load_path,obj_id_save)
eph_fname = "{}/df_eph_{}.csv".format(eph_load_path,obj_id_save)

print(data_fname)
print(eph_fname)

# delete existing data to force refresh
if reload:
    for f in [data_fname,eph_fname]:
        try:
            os.remove(f)
            print("delete {}".format(f))
        except:
            continue

# TRY LOAD AN EXISTING FILE
if os.path.isfile(data_fname):
    print("load data")
    df_data = pd.read_csv(data_fname,index_col=0)
else:
    print("query data")
    # connect to database
    cnx=atlas_database_connection.database_connection().connect()

    # get orbit_id
    orbid=atlas_SQL_query_df.get_orb_elements_id(cnx,mpc_number,name)

    # do the query
    df_data=atlas_SQL_query_df.atlas_SQL_query_orbid_expname(cnx,orbid)

    # save df
    df_data.to_csv(data_fname)

print(df_data)
print(list(df_data))

# drop any nans
df_data = df_data.dropna()

# select only a particular filter
if filter:
    df_data=df_data[df_data["filter"]==filter]

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

# plot galactic latitude also

# Get the JPL Horizons elongation

# TRY LOAD AN EXISTING FILE
if os.path.isfile(eph_fname):
    print("load Horizons")
    df_eph = pd.read_csv(eph_fname,index_col=0)
else:
    print("query Horizons")
    t1=Time(int(np.amin(df_data["mjd"])),format="mjd")
    t2=Time(int(np.amax(df_data["mjd"])),format="mjd")
    epoch_list = {'start':t1.iso, 'stop':t2.iso, 'step':'1d'} # a range of epochs in Horizons format is FAST!
    obj = Horizons(id=obj_id, location=loc, epochs=epoch_list)
    eph = obj.ephemerides()
    df_eph = eph.to_pandas()
    df_eph["mjd"] = df_eph["datetime_jd"] - mjd_jd
    # save df
    df_eph.to_csv(eph_fname)

# plot the exact elongation from JPL
ax3.plot(df_eph["mjd"],df_eph["elong"],c="r",alpha=0.2,label = "JPL elong ({})".format(loc),zorder=0)
# mask = (df_eph["elong"]<10.0)
# ax3.scatter(df_eph[mask]["mjd"],df_eph[mask]["elong"],c="r",s=1,label = "JPL elong<10 deg",zorder=0)


# find the turning points of the solar elongation data
ydiff = np.diff(ydata)

print("ydiff: {}".format(np.sort(ydiff)[:10]))

# find any -ve drop
# turning_points = xdata[1:][ydiff<0]

# select only drops of at least a certain size
turning_points = xdata[1:][ydiff<sol_elong_diff]

# # select also on gaps in time data
# xdiff = np.diff(xdata)
# # search for drops in elongation and large gaps
# turning_points = xdata[1:][(ydiff<sol_elong_diff) | (xdiff>date_diff)]

# add the first and last mjds in the dataset to fully bound all apparitions
turning_points = np.insert(turning_points,0,xdata[0])
turning_points = np.insert(turning_points,len(turning_points),xdata[-1])
print(turning_points)
for i,t in enumerate(turning_points):
    if i==0 or i==len(turning_points)-1:
        ax3.axvline(t,c="r",ls=":",alpha=0.3)
    else:
        ax3.axvline(t,c="r")

app_ranges = []
# use these turning points to divide into epochs
for i in range(1,len(turning_points)):
    t1=turning_points[i-1]
    t2=turning_points[i]

    df = df_data[(df_data["mjd"]>=t1) & (df_data["mjd"]<t2)]
    app_ranges.append(t2-t1)

    # ax2.scatter(df['mjd'],df['reduced_mag'],s=1,c="C{}".format(i))
    ax2.scatter(df['mjd'],df['reduced_mag'],s=2)
    ax1.scatter(df['phase_angle'],df['reduced_mag'],s=2,label="epoch:{}".format(i))

N_app = len(turning_points)-1
mjd_range = xdata[-1]-xdata[0]
print(app_ranges)
app_range_max = np.amax(app_ranges)
print("estimated apparitions = {}".format(mjd_range/P_est))
print("measured apparitions = {}".format(N_app))
print("max epoch range = {}".format(app_range_max))
print("max range / N_app = {}".format(mjd_range/app_range_max))

ax1.legend(fontsize=6)
ax3.legend(fontsize=6)

title = "{}_{}".format(os.path.basename(__file__).split('.')[0],obj_id_save)
fig.suptitle(title)
plt.tight_layout()

fname="{}/{}.png".format(save_path,title)
print(fname)
plt.savefig(fname, bbox_inches='tight')

plt.close()
# plt.show()
