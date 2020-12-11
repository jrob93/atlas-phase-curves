from calculate_phase import atlas_SQL_query_df
from calculate_phase import atlas_database_connection

# cnx=atlas_database_connection.database_connection().connect()
# print(cnx)
# print(atlas_SQL_query_df.atlas_SQL_query(mpc_number=4985,cnx=cnx))

from calculate_phase import sbpy_phase_fit

import cProfile
import time

# mpc_number_list=[4985]#,4986,4987,4988]
# mpc_number_list=[1]
mpc_number_list=[2,24,62,213,222,229,261,316,335,379,383,419,426,431,505,554,704,762,1508,2100,3200,6411] # B-types Clark et al 2010

start = time.process_time()
for m in mpc_number_list:
    fit = sbpy_phase_fit.phase_fit(m,push_fit_flag=False,plot_fig_flag=True,save_fig_flag=True,save_path="calculate_phase/figs")
    fit.calculate()
end = time.process_time()


# start = time.process_time()
# fits=[]
# for m in mpc_number_list:
#     fits.append(sbpy_phase_fit.phase_fit(m,push_fit_flag=False))
# for f in fits:
#     f.calculate()
#
# end = time.process_time()

print((end-start)/len(mpc_number_list))
