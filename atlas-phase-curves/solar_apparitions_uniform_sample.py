# see also solar_apparitions_sample.ipynb
import pandas as pd
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import tools.database_tools as dbt
import calculate_phase.solar_apparitions as sa
import os

fpath = "/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/results_analysis/fit_db_analysis"
fname = "atlas_phase_fits_orbs_2_7_2021.csv"
fname_results = "/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/solar_apparitions_files/df_atlas_objects_uniform_sample_results.csv"
sample_file = "/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/solar_apparitions_files/df_atlas_objects_uniform_sample.csv"
data_load_path="/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/results_analysis/obs"
eph_load_path="/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/solar_apparitions_files"

print("load sample")
df_samp = pd.read_csv(sample_file,index_col=0)
sol_elong_diff = -1.0

result_list = []

for i in range(len(df_samp)):
    name = df_samp.iloc[i]["name"]
    mpc_number = df_samp.iloc[i]["mpc_number"]

    # sol = sa.solar_apparitions(name=name)
    sol = sa.solar_apparitions(name=name,data_load_path=data_load_path,eph_load_path=eph_load_path)
    tp1 = sol.solar_elongation_JPL(JPL_step="7d")
    tp2 = sol.solar_elongation(sol_elong_diff)

    if len(tp1)==0:
        N_app_JPL = np.nan
    else:
        N_app_JPL = len(tp1)-1
    N_app_diff = len(tp2)-1

    orb=list(df_samp.iloc[i][["a_semimajor_axis","e_eccentricity","i_inclination_deg","detection_count"]])
    results = [mpc_number,name,N_app_diff,N_app_JPL]+orb
    print(results)

    result_list.append(results)

    # break

columns = ["mpc_number","name",
             "N_app_diff","N_app_JPL",
                "a_semimajor_axis","e_eccentricity","i_inclination_deg",
                "detection_count"]

df_results = pd.DataFrame(result_list, columns = columns)

print("save {}".format(fname_results))
df_results.to_csv(fname_results)
