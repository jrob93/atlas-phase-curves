""" This script calls the phase fitting functions, and accepts command line options """

from calculate_phase import sbpy_phase_fit
from optparse import OptionParser

parser = OptionParser()
parser.add_option( "-n", "--mpc-number", dest="mpc_number", help="mpc_number", metavar="MPC_NUMBER" ) # mpc number of object to fit
parser.add_option( "-N", "--name", dest="name", help="name", metavar="NAME" ) # name of object to fit
parser.add_option( "-s", "--save-path", dest="save_path", help="save_path", metavar="SAVE_PATH" ) # path to save figures
parser.add_option( "-f", "--filter", dest="filter", help="filter", metavar="FILTER" ) # filters to use, declare o or c, default is both
parser.add_option( "-w", "--warnings", dest="warnings", help="warnings", metavar="WARNINGS" ) # Suppress warnings?
parser.add_option( "-i", "--start-date", dest="start_date", help="start_date", metavar="START_DATE" ) # starts date for data
parser.add_option( "-j", "--end-date", dest="end_date", help="end_date", metavar="END_DATE" ) # starts date for data

(options,args)=parser.parse_args()

if options.mpc_number:
    mpc_number=int(options.mpc_number)
else:
    mpc_number=False
if options.name:
    name=str(options.name)
else:
    name=False
if options.save_path:
    save_path=options.save_path
else:
    save_path="."
if options.filter:
    filters=list(options.filter)
else:
    filters=["o","c"]
if options.warnings:
    warning_flag=int(options.warnings)
else:
    warning_flag=False
if options.start_date:
    start_date=float(options.start_date)
else:
    start_date=False
if options.end_date:
    end_date=float(options.end_date)
else:
    # end_date=False
    end_date=59259 # 2021-02-14 00:00:00.000

# add options for push and plot flags

fit = sbpy_phase_fit.phase_fit(mpc_number,name,
push_fit_flag=False,plot_fig_flag=False,save_fig_flag=False,
# show_fig_flag=True,
hide_warning_flag=warning_flag,
save_path="results_analysis/calculate_phase_figs",
# save_path="results_analysis/calculate_phase_figs_0-01",
# save_path="results_analysis/calculate_phase_figs_0-005",
start_date=start_date, end_date=end_date,
model_list=["HG"],filter_list=["o"],
# save_file_type="pdf",
mag_diff_flag=True,save_file_suffix="_mag_diff"
)
df=fit.calculate()
print(df.iloc[0].to_string())
