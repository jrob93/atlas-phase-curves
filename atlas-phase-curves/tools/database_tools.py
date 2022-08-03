import numpy as np
import pandas as pd
import glob
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sbpy.data import Names
import os
from pathlib import Path

# define the path to this file so that database_meta_files can be loaded
# figure out how to set this automatically when importing
# abs_path="/Users/jrobinson/atlas-phase-curves/atlas-phase-curves/tools"
abs_path = "{}".format(Path(__file__).parent)
# print(abs_path)

#-------------------------------------------------------------------------------
# database loading functions
#-------------------------------------------------------------------------------

def load_atlas_phase_fits(fname,nrows=None):
    """
    Load our phase fit database.
    Uses the dtypes file to deal with a pandas low memory error

    fname = file name of the database csv, including the path
    nrows = number of rows to read in, leave as None to read whole file
    """

    # load the dtypes for the atlas_phase_fits table to avoid low memory error from read_csv
    # NB pandas has trouble with nans in int columns - use float - https://stackoverflow.com/questions/11548005/numpy-or-pandas-keeping-array-type-as-integer-while-having-a-nan-value
    dtypes_file="{}/database_meta_files/dtypes_atlas_phase_fits.txt".format(abs_path)
    with open(dtypes_file,"r") as f:
        dat=f.readlines()
    # parse the dtypes
    dat=[d.rstrip() for d in dat]
    keys=[k.split(" ")[0] for k in dat]
    vals=[k.split(" ")[-1].replace(",", "") for k in dat]
    # Create a zip object from two lists
    zipbObj = zip(keys,vals)
    # Create a dictionary from zip object
    dtype_dict = dict(zipbObj)

    # load as dataframe
    df=pd.read_csv(fname,index_col=0,dtype=dtype_dict,nrows=nrows)

    return df

def load_atlas_phase_fits_orbs(fname,nrows=None):
    """
    Load our phase fit database, with additional columns from the orbital_elements table (astorb).
    Uses the dtypes file to deal with a pandas low memory error

    fname = file name of the database csv, including the path
    nrows = number of rows to read in, leave as None to read whole file
    """

    # load the dtypes for the atlas_phase_fits table to avoid low memory error from read_csv
    # NB pandas has trouble with nans in int columns - use float - https://stackoverflow.com/questions/11548005/numpy-or-pandas-keeping-array-type-as-integer-while-having-a-nan-value
    dtypes_file="{}/database_meta_files/dtypes_atlas_phase_fits_orbs.txt".format(abs_path)
    with open(dtypes_file,"r") as f:
        dat=f.readlines()
    # parse the dtypes
    dat=[d.rstrip() for d in dat]
    keys=[k.split(" ")[0] for k in dat]
    vals=[k.split(" ")[-1].replace(",", "") for k in dat]
    # Create a zip object from two lists
    zipbObj = zip(keys,vals)
    # Create a dictionary from zip object
    dtype_dict = dict(zipbObj)

    # load as dataframe
    df=pd.read_csv(fname,index_col=0,dtype=dtype_dict,nrows=nrows)

    return df

def load_orbital_elements(fname,nrows=None):
    """
    Load our the rockAtlas orbital_elements table
    Uses the dtypes file to deal with a pandas low memory error

    fname = file name of the database csv, including the path
    nrows = number of rows to read in, leave as None to read whole file
    """

    # load the dtypes for the atlas_phase_fits table to avoid low memory error from read_csv
    # NB pandas has trouble with nans in int columns - use float - https://stackoverflow.com/questions/11548005/numpy-or-pandas-keeping-array-type-as-integer-while-having-a-nan-value
    dtypes_file="{}/database_meta_files/dtypes_orbital_elements.txt".format(abs_path)
    with open(dtypes_file,"r") as f:
        dat=f.readlines()
    # parse the dtypes
    dat=[d.rstrip() for d in dat]
    keys=[k.split(" ")[0] for k in dat]
    vals=[k.split(" ")[-1].replace(",", "") for k in dat]
    # Create a zip object from two lists
    zipbObj = zip(keys,vals)
    # Create a dictionary from zip object
    dtype_dict = dict(zipbObj)

    # load as dataframe
    df=pd.read_csv(fname,index_col=0,dtype=dtype_dict,nrows=nrows)

    return df

def atlas_phase_fits_new_cols(df,model="B89",filt="o"):

    """
    Calculate some new columns for the atlas fits dataframe
    """

    model_filt="{}_{}".format(model,filt)

    # add a new column - phase angle range for model_filt
    df.loc[:,"phase_angle_range_{}".format(model_filt)]=df["phase_curve_alpha_max_{}".format(model_filt)]-df["phase_curve_alpha_min_{}".format(model_filt)]
    # add a new column for fraction of good fits
    df.loc[:,"phase_curve_frac_good_fit_{}".format(model_filt)]=df["phase_curve_N_mag_err_{}".format(model_filt)]/df["phase_curve_N_fit_{}".format(model_filt)]

    return df

def load_mahlke_db(
    fname="/Users/jrobinson/asteroid_databases/mahlke2021_asteroids/phase.dat",
    nrows=None):

    """
    Load the Mahlke et al 2020 asteroid phase curve database.
    Requires the label and byte lists for pandas read fixed width file function

    fname = file name of the database csv, including the path
    nrows = number of rows to read in, leave as None to read whole file
    """

    label_list=np.load("{}/database_meta_files/label_list_mahlke.npy".format(abs_path))
    byte_start_list=np.load("{}/database_meta_files/byte_start_list_mahlke.npy".format(abs_path))
    byte_end_list=np.load("{}/database_meta_files/byte_end_list_mahlke.npy".format(abs_path))

    colspecs=[(i-1,j+1) for i,j in zip(byte_start_list,byte_end_list)]
    df=pd.read_fwf(fname,colspecs=colspecs,names=label_list,nrows=nrows)

    return df

def load_schemel_db(
    fname="/Users/jrobinson/asteroid_databases/schemel2021_trojans/schemel2021_trojans.txt",
    nrows=None):

    """
    """

    label_list=np.load("{}/database_meta_files/label_list_schemel2021.npy".format(abs_path))
    byte_list=np.load("{}/database_meta_files/byte_list_schemel2021.npy".format(abs_path))
    byte_list=[tuple(x) for x in byte_list] # np save and load doesn't preserve the tuples?

    df=pd.read_fwf(fname,skiprows=42,colspecs=byte_list,names=label_list,index_col=False,nrows=nrows)

    return df

def load_LCDB(
    fname="/Users/jrobinson/asteroid_databases/warner_LCDB/ast-lightcurve-database_V4_0/data/lc_summary.csv",
    nrows=None):

    """
    Data from https://alcdef.org/index.php, https://alcdef.org/docs/ALCDEF_ALL.zip
    See ast-lightcurve-database_V3_0/document/pds_readme.pdf for documentation
    """

    N_skip_lines=21
    header_line_number=17
    skip_list=np.delete(np.arange(N_skip_lines+1),header_line_number-1)

    df=pd.read_csv(fname,header=0,skiprows=skip_list,nrows=nrows)

    # I assume these are all placeholders for nan
    df=df.replace("-",np.nan)
    df=df.replace(-9,np.nan)
    df=df.replace(-9.99,np.nan)

    return df

def load_mpcorb_db(
    fname="/Users/jrobinson/asteroid_databases/mpcorb/MPCORB.DAT",
    nrows=None):

    """
    Load the mpcorb asteroid database.
    Requires the label and byte lists for pandas read fixed width file function

    fname = file name of the database csv, including the path
    nrows = number of rows to read in, leave as None to read whole file

    N.B. that if nrows is short the automatic pandas dtypes may not read the flags column properly
    """

    label_list=np.load("{}/database_meta_files/label_list_mpcorb.npy".format(abs_path))
    byte_list=np.load("{}/database_meta_files/byte_list_mpcorb.npy".format(abs_path))
    byte_list=[tuple(x) for x in byte_list] # np save and load doesn't preserve the tuples?

    # df=pd.read_fwf(fname,skiprows=43,colspecs=byte_list,names=label_list,nrows=nrows)

    # load the dtypes to avoid low memory error from read_csv
    # NB pandas has trouble with nans in int columns - use float - https://stackoverflow.com/questions/11548005/numpy-or-pandas-keeping-array-type-as-integer-while-having-a-nan-value
    dtypes_file="{}/database_meta_files/dtypes_MPCORB.txt".format(abs_path)
    with open(dtypes_file,"r") as f:
        dat=f.readlines()
    # parse the dtypes
    dat=[d.rstrip() for d in dat]
    keys=[k.split(" ")[0] for k in dat]
    vals=[k.split(" ")[-1].replace(",", "") for k in dat]
    # Create a zip object from two lists
    zipbObj = zip(keys,vals)
    # Create a dictionary from zip object
    dtype_dict = dict(zipbObj)

    # load as dataframe
    df=pd.read_fwf(fname,skiprows=43,colspecs=byte_list,names=label_list,nrows=nrows,dtype=dtype_dict)

    return df

def load_mpc_comets_db(
    fname="/Users/jrobinson/asteroid_databases/mpcorb/CometEls.txt",
    nrows=None):

    """
    Load the mpcorb comet database.
    Requires the label and byte lists for pandas read fixed width file function

    fname = file name of the database csv, including the path
    nrows = number of rows to read in, leave as None to read whole file

    N.B. that if nrows is short the automatic pandas dtypes may not read the flags column properly
    """

    label_list=np.load("{}/database_meta_files/label_list_mpc_comets.npy".format(abs_path))
    byte_list=np.load("{}/database_meta_files/byte_list_mpc_comets.npy".format(abs_path))
    byte_list=[tuple(x) for x in byte_list] # np save and load doesn't preserve the tuples?

    # load as dataframe
    df=pd.read_fwf(fname,colspecs=byte_list,names=label_list,nrows=nrows)

    return df

def convert_hex_digit(hexdigit):
    """
    mpcorb uses a packed designation for their orbit type flags
    This function unpacks the 4 hexdigit according to:

    https://www.minorplanetcenter.net/iau/info/MPOrbitFormat.html
    """

    result = {}

    if not pd.isna(hexdigit):
        hd = int(hexdigit, 16)
        result['orbitType'] = int(hd & 0X003F)
        result['isNeo'] = True if hd & 0X0800 else False
        result['NeoOver1km'] = True if hd & 0X1000 else False
        result['opposition'] = True if hd & 0X2000 else False
        result['critListNumbered'] = True if hd & 0X4000 else False
        result['isPHA'] = True if hd & 0X8000 else False
    else:
        result['orbitType'] = np.nan
        result['isNeo'] = np.nan
        result['NeoOver1km'] = np.nan
        result['opposition'] = np.nan
        result['critListNumbered'] = np.nan
        result['isPHA'] = np.nan

    return result

def get_orb_type(hexdigit):
    """
    Function to accept a hexdigit, unpack it, and return only the orbitType flag.
    Defined to be easy to apply to the mpcorb dataframe
    """
    result=convert_hex_digit(hexdigit)
    return result['orbitType']

"""this dict shows the conversion between mpcorb orbitType flag and the corresponding object class"""
mpcorb_label_dict={"undefined":0,"Atira":1,"Aten":2,"Apollo":3,"Amor":4,"q<1.665AU":5,"Hungaria":6,"unused/mpc":7,"Hilda":8,"Jupiter Trojan":9,"Distant object":10}

def return_orbitType(i_list=[0,1,2,3,4,5,6,7,8,9]):
    """
    i_list = list of mpcorb orbitType integers

    returns the corresponding name of the orbitType
    """

    key_list=list(mpcorb_label_dict.keys())
    labelled_list=[key_list[i_list.index(i)] for i in i_list]
    return labelled_list

def unpack_mpcorb_flag(df):
    """
    Apply orbitType unpack function to the mpcorb dataframe
    """

    df["orb_flag"]=df["flags"].apply(get_orb_type)
    return df

def _unpack_mpc_date(packed_date):
    """
    copied from
    https://github.com/dirac-institute/sso_tools/blob/master/sso_tools/catalogs/catalog_utils.py

    See https://minorplanetcenter.net/iau/info/PackedDates.html
    for MPC documentation on packed dates.
    Examples:
    1998 Jan. 18.73     = J981I73
    2001 Oct. 22.138303 = K01AM138303
    """
    packed_date = str(packed_date)
    year = _mpc_lookup(packed_date[0]) * 100 + int(packed_date[1:3]) # this used to be [1:2] which was wrong!
    # Month is encoded in third column.
    month = _mpc_lookup(packed_date[3])
    day = _mpc_lookup(packed_date[4])
    isot_string = '%d-%02d-%02d' % (year, month, day)
    if len(packed_date) > 5:
        fractional_day = float(packed_date[5:])
        isot_string += '.%f' % fractional_day
    t = Time(isot_string, format='isot', scale='tt')
    return t.mjd

def _mpc_lookup(x):
    """
    copied from
    https://github.com/dirac-institute/sso_tools/blob/master/sso_tools/catalogs/catalog_utils.py

    Convert the single character dates into integers.
    """
    try:
        x = int(x)
    except ValueError:
        x = ord(x) - 55
    if x < 0 or x > 31:
        raise ValueError
    return x

def mpc_get_epoch(packed_date):
    try:
        packed_date=packed_date.split()[0]
        d=_unpack_mpc_date(packed_date)
    except:
        d=np.nan
    return d

def unpack_mpcorb_epoch(df):
    """
    Apply the above functions to the mpcorb DataFrame
    THIS IS FAIRLY SLOW, TRY SPEED UP?
    ALSO, in some cases it introduced extra lines with nan?
    e.g. check around sed -n 547967,547969p MPCORB_unpacked.dat
    """
    df["epoch_mjd"]=df["epoch"].apply(mpc_get_epoch)
    return df

# def split_mpcorb_name(mpc_name):
#     """
#     Split the mpcorb name field to retrieve just the designation.
#     see also sbpy Names parse_asteroid
#     """
#
#     # catch any stray nans
#     try:
#         name_split=mpc_name.split(")")
#     except:
#         # return [np.nan,np.nan]
#         return np.nan
#
#     # get the name
#     name=name_split[-1].strip()
#
#     # # Also grab the mpc number, if it exists
#     # if len(name_split[0])>0:
#     #     mpc_number=name_split[0].split("(")[-1]
#     # else:
#     #     mpc_number = np.nan
#
#     # return [name,mpc_number]
#     return name
#
# def apply_split_mpcorb_name(df):
#     """
#     """
#
#     df[["name","mpc_number"]]=[split_mpcorb_name(x) for x in df["num_des"]]
#
#     return df

# def unpack_mpc_num_des(num_des):
#     """grab the number from the num_des column"""
#     x=Names.from_packed(str(num_des))
#
#     try:
#         num = int(x)
#     except:
#         num = np.nan
#
#     return num

def unpack_mpc_name_num(mpc_name):
    """
    Return either the name or designation from mpcorb name column, and the number if it exists.
    NB that sbpy Names also has from_packed + to_packed to read just the mpc numbers
    https://sbpy.readthedocs.io/en/latest/sbpy/data/names.html
    """

    try:
        x=Names.parse_asteroid(mpc_name)
    except:
        x=Names.parse_comet(mpc_name)

    # get the name if it exists, if not get the designation
    if 'name' in x:
        name=x['name']
    elif 'desig' in x:
        name=x['desig']
    else:
        name=np.nan

    # get the number if it exists
    if 'number' in x:
        number=x['number']
    else:
        number=np.nan

    return [name,number]

def mpcorb_load_unpack_save(fname="/Users/jrobinson/asteroid_databases/mpcorb/MPCORB.DAT",
    nrows=None,
    save_csv=False,
    save_path=None):

    """Perform all the loading and unpacking functions to mpcorb.
    Option to save to csv in order to save time when reloading
    """

    # Add option to download new MPCORB.DAT?

    print("load file")
    df = load_mpcorb_db(fname,nrows=nrows) # load the df directly from the mpcorb dat file
    print("unpack flags")
    df = unpack_mpcorb_flag(df) # unpack orbit flags

    # print("unpack epochs")
    df = unpack_mpcorb_epoch(df) # unpack the epochs THIS IS VERY SLOW!

    print("unpack names")
    df[["name","mpc_number"]] = [unpack_mpc_name_num(x) for x in df["num_name"]]

    # save the file to csv if desired
    if save_csv==True:
        if save_path:
            spath="{}/{}".format(save_path,fname.split("/")[-1].split(".dat")[0]+"_unpacked.dat")
        else:
            spath=fname.split(".DAT")[0]+"_unpacked.dat"

        print("save unpacked mpcorb: {}".format(spath))

        df.to_csv(spath)

    return df


def load_mpcorb_db_unpacked(
    fname="/Users/jrobinson/asteroid_databases/mpcorb/MPCORB_unpacked.dat"):

    """
    Load an already unpacked version of mpcorb,
    Note this file must be updated when a new mpcorb is downloaded
    """

    # load the dtypes to avoid low memory error from read_csv
    # NB pandas has trouble with nans in int columns - use float - https://stackoverflow.com/questions/11548005/numpy-or-pandas-keeping-array-type-as-integer-while-having-a-nan-value
    dtypes_file="{}/database_meta_files/dtypes_MPCORB.txt".format(abs_path)
    with open(dtypes_file,"r") as f:
        dat=f.readlines()
    # parse the dtypes
    dat=[d.rstrip() for d in dat]
    keys=[k.split(" ")[0] for k in dat]
    vals=[k.split(" ")[-1].replace(",", "") for k in dat]
    # Create a zip object from two lists
    zipbObj = zip(keys,vals)
    # Create a dictionary from zip object
    dtype_dict = dict(zipbObj)

    # load as dataframe

    # df=pd.read_csv(fname,index_col=0)
    df=pd.read_csv(fname,index_col=0,dtype=dtype_dict)

    return df

def load_astdys(fname="/Users/jrobinson/asteroid_databases/astdys/all.syn"):
    col_list1 = ["Name","mag.","a(AU)","e","sin(I)","n(deg/yr)",'g("/yr)','s("/yr)',"LCEx1E6","My"]
    byte_list1 = [(0,10),(10,18),(18,29),(29,39),(39,50),(50,63),(63,75),(75,88),(88,96),(96,98)]
    df=pd.read_fwf(fname,skiprows=2,names=col_list1,colspecs=byte_list1)
    return df
#-------------------------------------------------------------------------------
# database plotting functions
#-------------------------------------------------------------------------------

def plot_hist(df,col_name,n_bins=100,plot_log=True):
    """generate a single histogram"""

    fig = plt.figure()
    # fig.set_size_inches(10,4)
    gs = gridspec.GridSpec(1,1)

    # make a histogram for the df column
    a = plt.subplot(gs[0,0])
    a.hist(df[col_name],bins=n_bins,log=plot_log)
    a.set_xlabel(col_name)

    plt.tight_layout()
    plt.show()

def plot_hist_grid(df,x_fig,y_fig,col_names,n_bins,plot_log):
    """
    generate a large grid of histograms

    df = the dataframe to plot
    x_fig, y_fig = the dimensions of the grid of histograms
    col_names = list of columns names in df to plot as a grid (requires len(col_names)<=(x_fig * y_fig)
    n_bins = list of number of bins for each histogram
    plot_log = list of True/False values to determine if each hist should be log yscale or not
    """

    fig = plt.figure()
    fig.set_size_inches(2*y_fig,2*x_fig)
    gs = gridspec.GridSpec(y_fig,x_fig)

    count=0
    for i in range(y_fig):

        for j in range(x_fig):

            if count>len(col_names)-1:
                break

            col=col_names[count] # get the column name

            print(count,i,j,col)

            # make a histogram for each column
            a = plt.subplot(gs[i,j])
            a.hist(df[col],bins=n_bins[count],log=plot_log[count])
            a.set_xlabel(col)
            count+=1

    plt.tight_layout()

    fig.set_size_inches(12,8)
    fname="database_grid.pdf"
    # print(fname)
    # plt.savefig(fname, bbox_inches='tight')
    plt.show()

def plot_xyc(df,x_plot,y_plot,c_plot):
    """makes a matplotlib fig plotting x against y and coloured by c"""

    # order by colour property
    df=df.sort_values(c_plot)

    fig = plt.figure()
    gs = gridspec.GridSpec(1,1)
    ax1 = plt.subplot(gs[0,0])

    s1=ax1.scatter(df[x_plot],df[y_plot],c=df[c_plot],s=1)
    cbar1=fig.colorbar(s1)
    ax1.set_xlabel(x_plot)
    ax1.set_ylabel(y_plot)
    cbar1.set_label(c_plot)

    plt.tight_layout()
    plt.show()

def plot_xyc_hist(df,x_plot,y_plot,c_plot,n_bins=100,plot_log=True):
    """makes a matplotlib fig plotting x against y and coloured by c
    plot a hist dist for each axis"""

    # set up fig and grid
    fig = plt.figure()
    gs = gridspec.GridSpec(3,2,height_ratios=[0.5,1,0.1],width_ratios=[1,0.5])
    ax_fig = plt.subplot(gs[1,0])
    ax_cbar = plt.subplot(gs[2,0])
    ax_histx = plt.subplot(gs[0,0],sharex=ax_fig)
    ax_histy = plt.subplot(gs[1,1],sharey=ax_fig)
    fig.set_size_inches(5,5)

    # plot the figure and colour bar
    s1=ax_fig.scatter(df[x_plot],df[y_plot],c=df[c_plot],s=1)
    cbar1=fig.colorbar(s1,ax_cbar,use_gridspec=True,orientation='horizontal')
    ax_fig.set_xlabel(x_plot)
    ax_fig.set_ylabel(y_plot)
    cbar1.set_label(c_plot)

    # plot the histograms
    ax_histx.hist(df[x_plot],bins=n_bins,log=plot_log)
    ax_histy.hist(df[y_plot],bins=n_bins,log=plot_log,orientation="horizontal")
    ax_histx.set_ylabel("n")
    ax_histy.set_xlabel("n")

    plt.tight_layout()
    plt.show()

#-------------------------------------------------------------------------------
# Other tools
#-------------------------------------------------------------------------------

def col_sigma_clip(df,col,sigma=3):
    """
    perform sigma clipping on a particular col of a df.

    df = dataframe containing the data
    col = name of column in df to clip
    sigma = keep data points <+sigma*std and >-sigma*std

    returns the mask of values to keep, this mask must be applied to df
    """
    dat=df[col]
    std=np.std(dat)
    med=np.median(dat)
    # print(med,std)
    mask1=(dat>(med-(sigma*std)))
    mask2=(dat<(med+(sigma*std)))
    # print(sum(mask1),sum(mask2),sum(mask1 & mask2),sum(~(mask1 & mask2)))
    return mask1 & mask2
