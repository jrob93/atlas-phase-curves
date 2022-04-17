from calculate_phase import sbpy_phase_fit
import pandas as pd
import mysql.connector
import numpy as np
import time
from datetime import datetime
import multiprocessing
from multiprocessing import Pool
import subprocess
import os
import cProfile
import sys
from contextlib import redirect_stdout
from calculate_phase import atlas_database_connection
from calculate_phase.atlas_SQL_query_df import update_atlas_objects
from astropy.time import Time
import platform

# Set the required flags for the phase fit function
function_flags = {"push_fit_flag":True,"hide_warning_flag":False,"mag_diff_flag":True}
# function_flags = {"push_fit_flag":False,"hide_warning_flag":False,"mag_diff_flag":True}
print(function_flags)
# print(*function_flags)

# set the date beyond which no obs are counted
# end_date=Time(Time.now(),scale='utc',format="iso")
end_date=Time("2021-02-14",scale='utc',format="iso")
end_date_mjd=round(end_date.mjd)
print(end_date,end_date_mjd)

# define the fitting functions
def phase_fit_func_mpc(mpc_number,end_date):
    fit = sbpy_phase_fit.phase_fit(mpc_number=mpc_number,end_date=end_date,**function_flags)
    check=fit.calculate()
    return check

def phase_fit_func_name(name,end_date):
    fit = sbpy_phase_fit.phase_fit(name=name,end_date=end_date,**function_flags)
    check=fit.calculate()
    return check

# connect to the database
cnx=atlas_database_connection.database_connection().connect()

# perform the query and store results as a dataframe
qry="select mpc_number,name from atlas_objects;"
# qry="select mpc_number,name from atlas_objects limit 10;"
df=pd.read_sql_query(qry,cnx)
print(df)
cnx.close()

# make the mpc number list
df_mpc=df[~np.isnan(df["mpc_number"])]
mpc_number_list=df_mpc['mpc_number'].astype(int)
mpc_number_list=list(mpc_number_list)
print(mpc_number_list[:10])
print(len(mpc_number_list))

# make the name list
df_name=df[np.isnan(df["mpc_number"])]
name_list=list(df_name["name"])
print(name_list[:10])
print(len(name_list))

mpc_number_list = []
name_check = ['1998 QQ', '1999 FP19', '1999 LW1', '1999 TW16', '2000 EM26',
       '2000 SF8', '2003 JN14', '2004 TD10', '2005 GO59', '2006 KA40',
       '2006 NL', '2007 KD', '2008 CA6', '2008 PR9', '2009 CR4',
       '2009 EM1', '2009 HG', '2010 SA12', '2011 ED78', '2011 YB40',
       '2012 BB124', '2014 MH68', '2014 VU1', '2015 FO124', '2015 MW53',
       '2015 OL35', '2015 PR228', '2015 XL128', '2015 XJ261',
       '2015 XB379', '2016 AO165', '2016 AK193', '2016 BT', '2016 HC3',
       '2016 HJ19', '2016 JG18', '2016 LM8', '2016 LS9', '2016 LZ10',
       '2016 LY47', '2016 NV', '2016 NV15', '2016 OJ', '2016 PR8',
       '2016 PR38', '2016 UU80', '2016 VT2', '2017 FJ90', '2017 GM7',
       '2017 HH', '2017 HW48', '2017 JU2', '2017 KR34', '2017 LD',
       '2017 LW', '2017 MB1', '2017 MY2', '2017 MB3', '2017 NM6',
       '2017 OC', '2017 OD1', '2017 OP68', '2017 PS25', '2017 PJ26',
       '2017 PL26', '2017 QS17', '2017 QE18', '2017 QL18', '2017 QQ35',
       '2017 QR35', '2017 QT35', '2017 RU', '2017 RV1', '2017 RX1',
       '2017 RG2', '2017 RK2', '2017 RQ2', '2017 RS2', '2017 RU2',
       '2017 RV2', '2017 SN2', '2017 SP10', '2012 TT5', '2004 XP14',
       '2006 AL3', '2006 KZ86', '2009 DH39', '2013 KJ6', '2016 EO84',
       '2017 NK', '1996 AE2', '1998 HN3', '1999 DB2', '2000 EU70',
       '2000 WO148', '2002 JQ100', '2004 MD', '2005 MR5', '2006 DY',
       '2006 FH36', '2006 GX2', '2006 MA', '2006 OZ', '2006 PY17',
       '2008 KV2', '2008 LW16', '2008 YR27', '2008 YX32', '2010 WF3',
       '2010 XD11', '2013 CK89', '2013 PV2', '2016 AX147', '2000 TU28',
       '2002 JE9', '2003 YM1', '2003 YP17', '2004 MW2', '2004 QD3',
       '2004 QD17', '2004 RF84', '2007 VX137', '2009 BB', '2000 PE3',
       '2011 SO32']

print(len(mpc_number_list))
print(len(name_list))

for name in name_list:
    if name not in name_check:
        continue
    print(name)
    df = phase_fit_func_name(name,end_date_mjd)
    print(df)
    print("\n")
    break
