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
from msgr.utils.calc import compute_offset


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


def production_line_manager(configuration_file, the_date, outdir, radar_file_list, satdir, rid, gr_offset, pass_number=0):
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
    #  Reading configuration file
    config = configparser.ConfigParser()
    config.read(configuration_file)
    # Switch for writing out volume-matched data
    switch = config['switch']
    l_write = switch.getboolean('write')
    l_cband = switch.getboolean('cband')   # Switch for C-band GR
    l_dbz = switch.getboolean('dbz')       # Switch for averaging in dBZ
    l_gpm = switch.getboolean('gpm')       # Switch for GPM PR data
    l_atten = switch.getboolean('correct_gr_attenuation')

    outfilename = None

    print("")
    date = the_date.strftime("%Y%m%d")

    if len(radar_file_list) == 0:
        print_red("No radar file for this date %s." % (date))
        return None

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
        match_vol = cross_validation.match_volumes(configuration_file, radar_file_list, one_sat_file,
                                                   sat_file_2A25_trmm=fd_25, dtime=the_date, l_cband=l_cband,
                                                   l_dbz=l_dbz, l_gpm=l_gpm, l_atten=l_atten, gr_offset=gr_offset)

        if match_vol is None:
            continue

        print_with_time("Comparison done in %.2fs." % (time.time() - tick))

        # Saving data
        if l_write:
            # Output file name.
            outfilename = "RID_{}_ORBIT_{}_DATE_{}_PASS_{}.nc".format(rid, orbit, date, pass_number + 1)
            outfilename = os.path.join(outdir, outfilename)
            print_green("Saving data to {}.".format(outfilename), bold=True)
            save_data(outfilename, match_vol, the_date, offset=gr_offset, nb_pass=pass_number)

    return outfilename


def multiproc_manager(configuration_file, onedate, outdir, radar_file_list, satdir, rid, gr_offset):
    """
    Buffer function that handles Exceptions while running the multiprocessing.
    Automatically runs the comparison 2 times. First with the raw radar data,
    then it computes the offset between ground radars and satellites and runs
    a second time with that offset.
    """
    for c in range(2):
        try:
            outdata_file = production_line_manager(configuration_file, onedate, outdir, radar_file_list, satdir, rid, gr_offset, pass_number=c)
        except Exception:
            traceback.print_exc()
            return None

        if outdata_file is None:
            return None

        if not os.path.exists(outdata_file):
            return None

        gr_offset = compute_offset(outdata_file)
        if np.abs(gr_offset) < 1:
            print_green(f"No significant difference between ground radar and satellite found for {onedate}.")
            break
        elif np.isnan(gr_offset):
            print_red(f"Invalid offset found. Stopping comparison for this {onedate}.")
            return None
        elif c == 0:
            print_magenta(f"The difference between the ground radar data and the satellite data " +
                          f"for {onedate} is of {gr_offset:0.2} dB.")

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

    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

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
    args_list = []
    for onedate in date_list:
        mydate = onedate.strftime("%Y%m%d")
        radar_file_list = [f for f in total_radar_file_list if mydate in f]

        if len(radar_file_list) == 0:
            print_yellow(f"No ground radar file found for this date {mydate}")
            continue

        # Argument list for multiprocessing.
        args_list.append((CONFIG_FILE, onedate, outdir, radar_file_list, satdir, rid, gr_offset))

    if len(args_list) == 0:
        print_red("Nothing to do. Is configuration file correct?")
        return None

    # Start multiprocessing.
    with Pool(ncpu) as pool:
        pool.starmap(multiproc_manager, args_list)

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
