"""
Calls the atlas_SQL_query function to get an object's observations from ATLAS
Uses the solar elongation to analyse the apparitions of an object.
"""

import pandas as pd
from calculate_phase import atlas_SQL_query_df
from calculate_phase import atlas_database_connection
from optparse import OptionParser
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import scipy.signal
from astropy.timeseries import LombScargle

parser = OptionParser()
parser.add_option( "-n", "--mpc-number", dest="mpc_number", help="mpc_number", metavar="MPC_NUMBER" ) # mpc number of object to fit
parser.add_option( "-N", "--name", dest="name", help="name", metavar="NAME" ) # NAME of object to fit
parser.add_option( "-f", "--filter", dest="filter", help="filter", metavar="FILTER" ) # colour filter of data
(options,args)=parser.parse_args()

if options.mpc_number:
    mpc_number=int(options.mpc_number)
else:
    # mpc_number="4986"
    mpc_number=False
if options.name:
    name=options.name
else:
    name=False
if options.filter:
    filter=options.filter
else:
    filter=False

# def tan_func(x, A, P, phi, B):
#     return (A*np.tan((x/P)+phi))+B
def tan_func(x, A, P, phi):
    return A*np.tan((x/P)+phi) + 180.0
def sin_func(x, A, P, phi, B):
    return (A*np.sin((x/P)+phi))+B
def saw_func(x, A, fi, offset, width, phi):
    # width is from scipy docs at https://www.pydoc.io/pypi/scipy-1.0.1/autoapi/signal/waveforms/index.html#signal.waveforms.sawtooth
    return (A * scipy.signal.sawtooth((x / fi)+phi, width)) + offset
# def saw_func(x,fi, phi):
#     # width is from scipy docs at https://www.pydoc.io/pypi/scipy-1.0.1/autoapi/signal/waveforms/index.html#signal.waveforms.sawtooth
#     return (180.0 * scipy.signal.sawtooth((x / fi)+phi, 1)) + 180.0

# connect to database
cnx=atlas_database_connection.database_connection().connect()

# primid,orbid=atlas_SQL_query_df.get_unique_ids(cnx,mpc_number,name)
# print(primid,orbid)
ids=atlas_SQL_query_df.get_unique_ids(cnx,mpc_number,name)
print(ids)
print(ids["primaryId"])
print(ids["orbital_elements_id"])
# exit()

# # run the SQL query using mpc number
# start = time.process_time()
# df_data=atlas_SQL_query_df.atlas_SQL_query(mpc_number=mpc_number,cnx=cnx)
# time1=time.process_time() - start

# # time single function for orbid and data
# start = time.process_time()
# df_data=atlas_SQL_query_df.atlas_SQL_query_test(mpc_number=mpc_number,cnx=cnx,name=name)
# time2=time.process_time() - start

# time separate funcs
start = time.process_time()
orbid=atlas_SQL_query_df.get_orb_elements_id(cnx,mpc_number,name)

# do the regular query
# df_data=atlas_SQL_query_df.atlas_SQL_query_orbid(cnx,orbid)
# get the query with exposure name
df_data=atlas_SQL_query_df.atlas_SQL_query_orbid_expname(cnx,orbid)
time3=time.process_time() - start

print(df_data[np.absolute(df_data["galactic_latitude"])<10][["expname",'ra_deg','dec_deg']])
# save dataframe
print(df_data)
print(list(df_data))
print(df_data[["expname",'ra_deg','dec_deg']])

print(df_data[["m",'merr','galactic_latitude']])

df_data = df_data.dropna()

if filter:
    df_data=df_data[df_data["filter"]==filter]

# df_data.to_csv("results_analysis/obs/df_data_{}.csv".format(mpc_number))

# print(time1/time2,time2/time2, time3/time2)

# err_cut=0.01
err_cut=0.005
N_merr = len(df_data[df_data["merr"]<err_cut])
print("N merr<{} : {}".format(err_cut,N_merr))


# fig = plt.figure()
# gs = gridspec.GridSpec(1, 1)
# ax1 = plt.subplot(gs[0,0])
#
# # ax1.hist(df_data["merr"],bins="auto")
# # ax1.axvline(err_cut,c="r")
# # ax1.set_xlabel("N")
# # ax1.set_ylabel("merr")
#
# # ax1.scatter(df_data["m"],df_data["merr"],c=df_data["mjd"])
# # ax1.axhline(0.01,c="r",label="merr=0.01")
# # ax1.axvline(np.median(df_data["m"]),label="median m")
# # ax1.legend()
# # ax1.set_xlabel("m")
# # ax1.set_ylabel("merr")
#
# ax1.errorbar(df_data['phase_angle'],df_data['reduced_mag'],df_data['merr'], fmt='ko',label="data",zorder=0,markersize="2")
# ax1.legend()
# ax1.set_xlabel("phase_angle")
# ax1.set_ylabel("reduced_mag")
# ax1.invert_yaxis()
#
# plt.show()

fig = plt.figure()
# gs = gridspec.GridSpec(4, 1)
# ax1 = plt.subplot(gs[0,0])
# ax2 = plt.subplot(gs[1,0])
# ax4 = plt.subplot(gs[2,0])
# ax5 = plt.subplot(gs[3,0])
gs = gridspec.GridSpec(2, 2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax4 = plt.subplot(gs[1,1])
ax5 = plt.subplot(gs[1,0])

ax1.invert_yaxis()
ax1.errorbar(df_data['phase_angle'],df_data['reduced_mag'],df_data['merr'], fmt='ko',label="data",zorder=0,markersize="2")
ax1.legend()
ax1.set_xlabel("phase_angle")
ax1.set_ylabel("reduced_mag")

ax5.scatter(df_data['mjd'],df_data['reduced_mag'])
ax5.set_xlabel("mjd")
ax5.set_ylabel("reduced_mag")
ax5.invert_yaxis()

ax3 = ax5.twinx()
s=3
# ax3.scatter(df_data['mjd'],df_data['sun_obs_target_angle'],c="C1",label="solar elongation",s=s)
ax3.scatter(df_data['mjd'],df_data['sun_obs_target_angle']%360,c="C1",label="solar elongation",s=s)
# ax3.plot(df_data['mjd'],df_data['sun_obs_target_angle']%360,c="C1",label="solar elongation")
ax3.scatter(df_data['mjd'],df_data['phase_angle'],c="C2",label="phase angle",s=s)
ax3.set_ylabel("angle")

# df_data=df_data[(df_data["mjd"]>58000) & (df_data["mjd"]<58400)]
xdata=df_data["mjd"]
# ydata=df_data["sun_obs_target_angle"]
ydata=df_data["sun_obs_target_angle"]%360

# period=450.0
# do Lomb Scargle Periodogram, see https://docs.astropy.org/en/stable/timeseries/lombscargle.html
# get the starting guess of the period
# frequency, power = LombScargle(xdata, ydata).autopower()
frequency = np.linspace(1.0/1e3, 1.0/1e2, 1000)
power = LombScargle(xdata,ydata).power(frequency)
print(np.amin(frequency),np.amax(frequency))
best_frequency = frequency[np.argmax(power)]
period_days = 1. / frequency
print(np.amin(period_days),np.amax(period_days))
best_period = period_days[np.argmax(power)]
print("best frequency = {}, P=1/f = {}, P/2pi={}".format(best_frequency,best_period,best_period/(2.0*np.pi)))
phase = (np.array(xdata) / best_period) % 1
period=best_period

ax2.scatter(period_days,power)
ax2.axvline(best_period,c="r")
ax2.set_xlabel("period(d)")
ax2.set_ylabel("power")

ax4.scatter(phase,ydata)
ax4.set_xlabel("phase")
ax4.set_ylabel("angle")

x=np.linspace(np.amin(df_data["mjd"]),np.amax(df_data["mjd"]),1000)

# popt_guess = np.array([period/(2.0*np.pi), 0])
# print(popt_guess)
# y = saw_func(x, *popt_guess)
# popt, pcov = curve_fit(saw_func, xdata, ydata, p0=popt_guess)
# print(popt)
# y_fit=saw_func(x, *popt)
# ax3.plot(x, y_fit, c='k')

# pass a lambda to curve_fit to fit only the phase of the saw function
A=180.0
offset=180.0
width=1.0
fi=period/(2.0*np.pi)
phi=0
popt_guess = np.array([A, fi, offset, width, phi])
print(popt_guess)
y = saw_func(x, *popt_guess)
# popt, pcov = curve_fit(saw_func, xdata, ydata)#, p0=popt_guess)
popt,pcov = curve_fit(lambda x, phi: saw_func(x, A, fi, offset, width, phi), xdata, ydata) # https://stackoverflow.com/questions/12208634/fitting-only-one-parameter-of-a-function-with-many-parameters-in-python
print(popt)
popt = np.array([A, fi, offset, width, popt[0]]) # only phi is returned, recombine the parameters
print(popt)
y_fit=saw_func(x, *popt)
ax3.plot(x, y_fit, c='k')


# ax3.set_ylim(np.amin(ydata),np.amax(ydata))
# ax3.set_ylim(0,360)

y_fit_diff = np.diff(y_fit)
turning_points = x[1:][y_fit_diff<0]
# ax3.vlines(turning_points,-10,370,color="r",zorder=5)
for t in turning_points:
    ax3.axvline(t,c="r")
# add the first and last mjds in the dataset to fully bound all apparitions

# ax3.legend()
plt.tight_layout()

plt.show()

print()
