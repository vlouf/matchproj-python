#! /usr/bin/env python
"""
MSGR Matching Satellite and Ground Radar
========================================

@author: Valentin Louf
@date: 2016-12-06 (creation) 2017-10-05 (current version)
@email: valentin.louf@bom.gov.au
@company: Monash University/Bureau of Meteorology
"""
# Standard library import
import os
import re
import glob
import time
import argparse
import datetime
import warnings
import traceback
import configparser

from multiprocessing import Pool

# Others lib.
import numpy as np
import pandas as pd

# Custom lib.
from msgr import cross_validation
from msgr.utils.misc import *  # bunch of useful functions
from msgr.io.save_data import save_data


def get_orbit_number(infile):
    """
    Look for an orbit number in the given filename.

    Parameters:
    ===========
    infile: str
        Input file.

    Returns:
    ========
    orbit: str
        Supposed orbit number
    """
    orbit = re.findall("[0-9]{6}", infile)[-1]  # Get orbit number
    return orbit


def get_satfile_list(satdir, date, l_gpm):
    """
    Get a list of satellite files.

    Parameters:
    ===========
    satdir: str
        Path to satellite data directory.
    date: str
        Date, format YYYYMMDD
    l_gpm: bool
        Is this GPM or TRMM?

    Returns:
    ========
    satfiles: str
        List of GPM files or TRMM 2A23 files.
    satfiles2: str
        List of TRMM 2A23 files (None for GPM).
    """
    # Looking for satellite data files.
    if l_gpm:
        satfiles = glob.glob(satdir + '/*' + date + '*.HDF5')
        satfiles2 = None
    else:
        satfiles = glob.glob(satdir + '/*2A23*' + date + '*.HDF')
        satfiles2 = glob.glob(satdir + '/*2A25*' + date + '*.HDF')

    # Checking if found any satellite data file.
    if len(satfiles) == 0:
        satfiles, satfiles2 = None, None

    return satfiles, satfiles2


def production_line_manager(configuration_file, the_date, outdir, satdir, rid,
                            gr_offset, l_cband=True, l_dbz=True, l_gpm=True,
                            l_atten=True, l_write=True):
    """
    Here we locate the satellite files, call the comparison function
    match_volumes, and send the results for saving. The real deal is the
    match_volumes from msgr.core.msgr module.

    Parameters
    ==========
        the_date: datetime
            The day for comparison.
        parameters_dict: dict
            Dictionnary containing all parameters from the configuration file.
    """
    print("")
    date = the_date.strftime("%Y%m%d")

    # Looking for satellites.
    satfiles, satfiles2 = get_satfile_list(satdir, date, l_gpm)
    if satfiles is None:
        print_red("No satellite swaths for %s." % (date))
        return None

    if len(satfiles) > 5:
        print_red("There are more than 5 files for {}. Something probably wrong in the files name. Skipping this date.".format(date))
        return None

    # Looping over satellite file
    for one_sat_file in satfiles:
        # Start chrono
        tick = time.time()

        # Get orbit number
        orbit = get_orbit_number(one_sat_file)
        print_with_time("Orbit #{} -- {}.".format(orbit, date))

        if not l_gpm:
            try:
                # Trying to find corresponding 2A25 TRMM file based on the orbit
                fd_25 = find_file_with_string(satfiles2, orbit)
            except IndexError:
                print_red("Found no matching 2A25 file for TRMM.")
                continue
        else:
            fd_25 = None

        # Calling processing function for TRMM
        match_vol = cross_validation.match_volumes(configuration_file, None, one_sat_file,
                                                   sat_file_2A25_trmm=fd_25, dtime=the_date, l_cband=l_cband,
                                                   l_dbz=l_dbz, l_gpm=l_gpm, l_atten=l_atten)

    return None


def multiproc_manager(kwargs):
    """
    Buffer function that handles Exceptions while running the multiprocessing.
    All the arguments are identical to the
    """
    try:
        production_line_manager(*kwargs)
    except Exception:
        traceback.print_exc()
        pass

    return None


def main():
    """
    Reading general informations from configuration file like ncpu, start_date,
    end_date, all the switches. Loop over dates, spawn multiprocessing, and call the
    production_line_manager.
    """
    print_with_time("Loading configuration file.")
    #  Reading configuration file
    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)

    # General info.
    general = config['general']
    ncpu = general.getint('ncpu')
    date1 = general.get('start_date')
    date2 = general.get('end_date')

    # Data path
    path = config['path']
    raddir = path.get('ground_radar')
    satdir = path.get('satellite')
    outdir = path.get('output')

    # Check if dirs exist.
    if not os.path.isdir(raddir):
        print_red("Wrong radar directory in configuration file.")
        return None
    if not os.path.isdir(satdir):
        print_red("Wrong satellite directory in configuration file.")
        return None
    if not os.path.isdir(outdir):
        print_red("Wrong output directory in configuration file.")
        return None

    # Switch for writing out volume-matched data
    switch = config['switch']
    l_write = switch.getboolean('write')
    l_cband = switch.getboolean('cband')   # Switch for C-band GR
    l_dbz = switch.getboolean('dbz')       # Switch for averaging in dBZ
    l_gpm = switch.getboolean('gpm')       # Switch for GPM PR data
    l_atten = switch.getboolean('correct_gr_attenuation')

    # General info about the ground radar (ID and OFFSET to apply.)
    GR_param = config['radar']
    rid = GR_param.get('radar_id')
    try:
        gr_offset = GR_param.getfloat('offset')
    except KeyError:
        gr_offset = 0

    start_date = datetime.datetime.strptime(date1, '%Y%m%d')
    end_date = datetime.datetime.strptime(date2, '%Y%m%d')

    print_yellow("Generating ground radar file list.")
    total_radar_file_list = get_files(raddir)
    print_yellow("Found {} supported radar files in {}.".format(len(total_radar_file_list), raddir))

    date_list = pd.date_range(start_date, end_date)
    args_list = [None] * len(date_list)
    for cnt, onedate in enumerate(date_list):
        mydate = onedate.strftime("%Y%m%d")
        # radar_file_list = [f for f in total_radar_file_list if mydate in f]

        # Argument list for multiprocessing.
        args_list[cnt] = (CONFIG_FILE, onedate, outdir, satdir, rid, gr_offset,
                          l_cband, l_dbz, l_gpm, l_atten, l_write)

    # Start multiprocessing.
    with Pool(ncpu) as pool:
        pool.map(multiproc_manager, args_list)

    return None


if __name__ == '__main__':
    """
    Global variables definition.

    Parameters:
    ===========
        CONFIG_FILE: str
            Configuration file (.ini)
    """
    parser_description = "Start MSGR - volume Matching Satellite and Ground Radar."
    parser = argparse.ArgumentParser(description=parser_description)

    parser.add_argument('-c', '--config', type=str, dest='config_file', help='Path to configuration file.', default=None, required=True)

    args = parser.parse_args()
    CONFIG_FILE = args.config_file

    warnings.simplefilter("ignore")
    main()
