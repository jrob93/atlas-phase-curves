""" This script calls the phase fitting functions, and accepts command line options """

from calculate_phase import sbpy_phase_fit
# import datetime
from astropy.time import Time
import numpy as np

# set an end date so that no observations after this time are selected
# utc_date_now=datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
# print(utc_date_now)
utc_now=Time(Time.now(),scale='utc')
print(utc_now)
print(utc_now.jd)
print(utc_now.mjd)

end_date=utc_now.mjd
# end_date=Time("2021-02-01",scale='utc').mjd

print(end_date)
# exit()

mpc_numbers_14_1_21 = np.array([    42,     68,     80,    115,    196,    214,    245,   3451,
         3995,   4833,   8491,   8734,  10424,  11129,  17309,  17707,
        18973,  22758,  24007,  25090,  27035,  28452,  29421,  30070,
        31912,  33149,  33492,  36369,  38818,  48769,  49637,  50768,
        50801,  51509,  54298,  55068,  55709,  59534,  62132,  64033,
        66378,  71618,  72096,  75739,  80968,  86578,  86608,  87089,
        87855,  88308,  88748,  91585,  92551,  93481,  93649,  95267,
        97181,  97192,  99572, 100150, 101080, 101981, 101986, 106227,
       111093, 111818, 112556, 112727, 113743, 118832, 122594, 123731,
       125737, 130688, 130707, 134649, 136629, 150621, 155506, 168576,
       170233, 171220, 194549, 199566, 211157,  70716, 141570, 157745,
       200054])

mpc_numbers_29_1_21=np.array([    68,     80,    115,    196,    214,    245,   3451,   3995,
         4833,   8491,   8734,  10424,  11129,  17707,  18973,  22758,
        24007,  25090,  27035,  28452,  29421,  30070,  31912,  33149,
        33492,  36369,  38818,  48769,  49637,  50768,  50801,  51509,
        54298,  55068,  55709,  59534,  62132,  64033,  66378,  71618,
        72096,  75739,  80968,  86578,  86608,  87089,  87855,  88308,
        88748,  91585,  92551,  93481,  93649,  95267,  97181,  97192,
        99348,  99572, 100150, 101080, 101981, 101986, 106227, 111093,
       111818, 112556, 112727, 113743, 118832, 122594, 123731, 125737,
       130688, 130707, 134649, 136629, 150621, 155506, 168576, 170233,
       171220, 191393, 194549, 199566, 211157,  70716, 141570, 157745,
       200054])

mpc_numbers_1_2_21=np.array([    80,    115,    196,    214,    245,   3451,   3995,   4833,
         8491,   8734,  10424,  11129,  17707,  18973,  22758,  24007,
        25090,  27035,  28452,  29421,  30070,  31912,  33149,  33492,
        36369,  38818,  48769,  49637,  50768,  50801,  51509,  54298,
        55068,  55709,  59534,  62132,  64033,  66378,  71618,  72096,
        75739,  80968,  86578,  86608,  87089,  87855,  88308,  88748,
        91585,  92551,  93481,  93649,  95267,  97181,  97192,  99348,
        99572, 100150, 101080, 101981, 101986, 106227, 111093, 111818,
       112556, 112727, 113743, 118832, 122594, 123731, 125737, 130688,
       130707, 134649, 136629, 150621, 155506, 168576, 170233, 171220,
       191393, 194549, 199566, 211157,  70716, 141570, 157745, 200054])

# print(len(mpc_numbers_14_1_21),len(mpc_numbers_14_1_21))
# mask_14=np.isin(mpc_numbers_14_1_21,mpc_numbers_29_1_21)
# print(mask_14)
# print(mpc_numbers_14_1_21[~mask_14])
#
# mask_29=np.isin(mpc_numbers_29_1_21,mpc_numbers_14_1_21)
# print(mask_29)
# print(mpc_numbers_29_1_21[~mask_29])

mpc_numbers_5_2_21=[    1291,
    2478,
   33066,
  112357]

mpc_numbers=mpc_numbers_5_2_21
# exit()

ratio_list=[]

# mpc_numbers=[1291]

# wipe the log file
with open("sbpy_phase_fit.log","w") as f:
    f.close()

# mpc_number = 4986
# name="osipovia"

# mpc_number = False
# name="2017 QN14"

for mpc_number in mpc_numbers:
    # mpc_number=42
    name=False

    fit = sbpy_phase_fit.phase_fit(mpc_number=mpc_number,name=name,
    push_fit_flag=False,plot_fig_flag=False,save_fig_flag=False,
    # push_fit_flag=True,plot_fig_flag=True,save_fig_flag=True,
    # filter_list=["o"],
    # save_path="results_analysis/calculate_phase_figs",
    # save_path="50801_no_end_date",
    # save_path="50801_end_date",
    end_date=end_date,
    hide_warning_flag=False)

    df=fit.calculate()

    print(df.iloc[0].to_string())
    ratio = df.iloc[0]["phase_curve_N_mag_err_B89_o"]/df.iloc[0]["phase_curve_N_fit_B89_o"]
    print(ratio)

    ratio_list.append(ratio)

    if ratio>1:
        print(mpc_number)
        break

for m,r in zip(mpc_numbers,ratio_list):
    print(m,r)
