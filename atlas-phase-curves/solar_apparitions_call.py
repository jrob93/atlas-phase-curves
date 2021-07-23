import pandas as pd
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import tools.database_tools as dbt
import os

plot=False
# plot=True

# sol_elong_diff = 10.0 # difference required for an epoch, should be minimised! (e.g. 76820)
# sol_elong_diff = 5.0
sol_elong_diff = 1.0

# sample full distribution
fpath = "/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/results_analysis/fit_db_analysis"
fname = "atlas_phase_fits_orbs_2_7_2021.csv"
sample_file = "solar_apparitions_files/df_atlas_objects_sample.csv"
# save_path = "solar_apparitions_call_figs"
# save_path = "solar_apparitions_call_figs_5"
save_path = "solar_apparitions_call_figs_1"
N_samp = int(1e3)

# # sample apollo asteroids
# fpath = "/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/results_analysis/fit_db_analysis"
# fname = "df_atlas_apollos.csv"
# sample_file = "solar_apparitions_files/df_atlas_objects_apollos_sample.csv"
# save_path = "solar_apparitions_call_figs_apollos"
# N_samp = int(1e2)

# # sample eccentric asteroids
# fpath = "/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/results_analysis/fit_db_analysis"
# fname = "df_atlas_eccentrics.csv"
# sample_file = "solar_apparitions_files/df_atlas_objects_eccentrics_sample.csv"
# save_path = "solar_apparitions_call_figs_eccentrics"
# N_samp = int(1e2)

# sample brightness limited objects

# Load or generate a subsample of objects
if os.path.isfile(sample_file):
    print("load sample")
    df_samp = pd.read_csv(sample_file,index_col=0)
else:
    print("gen sample")

    df=dbt.load_atlas_phase_fits_orbs("{}/{}".format(fpath,fname))

    # select only inner solar system
    J_aphelion=5.46 # Jupiter aphelion
    mask=(df["a_semimajor_axis"]>J_aphelion)
    df = df[~mask]

    df_samp = df.sample(N_samp)
    df_samp.to_csv(sample_file)
    print("save {}".format(sample_file))

if plot:
    # plot orbital distribution
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0])

    ax1.scatter(df["a_semimajor_axis"],df["e_eccentricity"],s=1)
    ax2.scatter(df["a_semimajor_axis"],df["i_inclination_deg"],s=1)

    ax1.scatter(df_samp["a_semimajor_axis"],df_samp["e_eccentricity"],c="r",marker="x")
    ax2.scatter(df_samp["a_semimajor_axis"],df_samp["i_inclination_deg"],c="r",marker="x")

    ax1.set_xlabel("a_semimajor_axis")
    ax1.set_ylabel("e_eccentricity")

    ax2.set_xlabel("a_semimajor_axis")
    ax2.set_ylabel("i_inclination_deg")

    plt.tight_layout()
    plt.show()

    exit()

# find the apparitions of sample
mask = np.isnan(df_samp["mpc_number"])
mpc_list = np.array(df_samp[~mask]["mpc_number"]).astype(int)
name_list = df_samp[mask]["name"]
# name_list = np.array(['"{}"'.format(n) for n in name_list])

# mpc_list = mpc_list[mpc_list >= 495927]
# mpc_list = []
# mpc_list = mpc_list[0]
# name_list = []
print("objects with mpc_number: {}".format(len(mpc_list)))
print("objects with name: {}".format(len(name_list)))
# exit()

for cmd_flag,obj_id_list in zip(["n","N"],[mpc_list,name_list]):
    for obj_id in obj_id_list:
        cmd = "python solar_apparitions_diff.py -{} \"{}\" -s {} -d {}".format(cmd_flag, obj_id, save_path, sol_elong_diff)
        print(cmd)
        p=subprocess.Popen(cmd,shell=True).wait()
