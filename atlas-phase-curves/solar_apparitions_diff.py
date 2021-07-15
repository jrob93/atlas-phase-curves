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
    mpc_number=False
if options.name:
    name=options.name
else:
    name=False
if options.filter:
    filter=options.filter
else:
    filter=False

def tan_func(x, A, P, phi):
    return A*np.tan((x/P)+phi) + 180.0
def sin_func(x, A, P, phi, B):
    return (A*np.sin((x/P)+phi))+B
def saw_func(x, A, fi, offset, width, phi):
    # width is from scipy docs at https://www.pydoc.io/pypi/scipy-1.0.1/autoapi/signal/waveforms/index.html#signal.waveforms.sawtooth
    return (A * scipy.signal.sawtooth((x / fi)+phi, width)) + offset

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
gs = gridspec.GridSpec(2, 2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax4 = plt.subplot(gs[1,1])
ax5 = plt.subplot(gs[1,0])

# plot the phase angle/reduced mag data
ax1.invert_yaxis()
ax1.errorbar(df_data['phase_angle'],df_data['reduced_mag'],df_data['merr'], fmt='ko',label="data",zorder=0,markersize="2")
ax1.legend()
ax1.set_xlabel("phase_angle")
ax1.set_ylabel("reduced_mag")

# plot reduced mag as a function of mjd date
ax5.scatter(df_data['mjd'],df_data['reduced_mag'])
ax5.set_xlabel("mjd")
ax5.set_ylabel("reduced_mag")
ax5.invert_yaxis()

# plot the solar elongation retrieved from rockAtlas
# N.B. rockAtlas angle cycles from -180 to 180 deg, use mod to get 0 to 360
ax3 = ax5.twinx()
s=3
# ax3.scatter(df_data['mjd'],df_data['sun_obs_target_angle'],c="C1",label="solar elongation",s=s)
ax3.scatter(df_data['mjd'],df_data['sun_obs_target_angle']%360,c="C1",label="solar elongation",s=s)
# ax3.plot(df_data['mjd'],df_data['sun_obs_target_angle']%360,c="C1",label="solar elongation")
ax3.scatter(df_data['mjd'],df_data['phase_angle'],c="C2",label="phase angle",s=s)
ax3.set_ylabel("angle")

# select data to fit
# df_data=df_data[(df_data["mjd"]>58000) & (df_data["mjd"]<58400)]
xdata=df_data["mjd"]
# ydata=df_data["sun_obs_target_angle"]
ydata=df_data["sun_obs_target_angle"]%360

# do Lomb Scargle Periodogram, see https://docs.astropy.org/en/stable/timeseries/lombscargle.html
# to get the period of the solar elongation data
# frequency, power = LombScargle(xdata, ydata).autopower() # autopower may not get long periods (typically 400-500 days)
frequency = np.linspace(1.0/1e3, 1.0/1e2, 1000) # manually set the period (frequency) search array
power = LombScargle(xdata,ydata).power(frequency)
print(np.amin(frequency),np.amax(frequency))
best_frequency = frequency[np.argmax(power)]
period_days = 1. / frequency
print(np.amin(period_days),np.amax(period_days))
best_period = period_days[np.argmax(power)]
print("best frequency = {}, P=1/f = {}, P/2pi={}".format(best_frequency,best_period,best_period/(2.0*np.pi)))
phase = (np.array(xdata) / best_period) % 1
period=best_period

# plot power spectrum
ax2.scatter(period_days,power)
ax2.axvline(best_period,c="r")
ax2.set_xlabel("period(d)")
ax2.set_ylabel("power")

# plot phase folded data
ax4.scatter(phase,ydata)
ax4.set_xlabel("phase")
ax4.set_ylabel("angle")

# define x axis grid to plot saw tooth function
x=np.linspace(np.amin(df_data["mjd"]),np.amax(df_data["mjd"]),1000)

# pass a lambda to curve_fit to fit only the phase of the saw function
# fix the following parameters, including the period found by LS
# We define a saw tooth that is always rising and drops instantly (width=1)
# offset and amplitude (A) are set such that saw rises from 0 -> 360 deg
A=180.0
offset=180.0
width=1.0
fi=period/(2.0*np.pi)
# we use scipy curve fit to change only the phase of the saw tooth function
phi=0
# define the initial parameters, popt_guess[-1] (phi) is the only one we will actually fit
popt_guess = np.array([A, fi, offset, width, phi])
print(popt_guess)
y = saw_func(x, *popt_guess) # initial guess
# popt, pcov = curve_fit(saw_func, xdata, ydata)#, p0=popt_guess)
# use python lambda such that curve_fit only changes phi (NB we do not pass the starting guess for phi)
popt,pcov = curve_fit(lambda x, phi: saw_func(x, A, fi, offset, width, phi), xdata, ydata) # https://stackoverflow.com/questions/12208634/fitting-only-one-parameter-of-a-function-with-many-parameters-in-python
print(popt)
# only phi is returned from curve_fit, recombine with the fixed parameters IN CORRECT ORDER
popt = np.array([A, fi, offset, width, popt[0]])
print(popt)
y_fit=saw_func(x, *popt)
ax3.plot(x, y_fit, c='k')

# find the turning points of the fitted saw tooth function
# when diff is -ve the saw tooth has dropped from 360 to 0 deg and this is the epoch bounds
y_fit_diff = np.diff(y_fit)
turning_points = x[1:][y_fit_diff<0]
# add the first and last mjds in the dataset to fully bound all apparitions
turning_points = np.insert(turning_points,0,x[0])
turning_points = np.insert(turning_points,len(turning_points),x[-1])
print(turning_points)
for i,t in enumerate(turning_points):
    if i==0 or i==len(turning_points)-1:
        ax3.axvline(t,c="r",ls=":",alpha=0.3)
    else:
        ax3.axvline(t,c="r")

# ax3.legend()
plt.tight_layout()

plt.show()
