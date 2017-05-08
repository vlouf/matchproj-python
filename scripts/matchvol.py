#! /usr/bin/env python
"""
MSGR Matching Satellite and Ground Radar
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@date: 2016-12-06 (creation) 2017-3-21 (current version)
@email: valentin.louf@bom.gov.au
@company: Monash University/Bureau of Meteorology
"""

__title__ = 'matchvol'
__version__ = '0.6'
__author__ = 'Valentin Louf'
__license__ = 'MIT'
__copyright__ = 'Copyright 2017 Valentin Louf'

# Standard library import
import os
import sys
import glob
import time
import argparse
import datetime
import configparser
from multiprocessing import Pool

# Others lib.
import pyart
import pyproj  # For cartographic transformations and geodetic computations
import numpy as np
import pandas as pd

# Custom lib.
import msgr

from msgr import config_codes
from msgr import processing_codes
from msgr.util_fun import *  # bunch of useful functions
from msgr.io.save_data import save_data
from msgr.instruments.satellite import get_orbit_number, satellite_params
from msgr.instruments.ground_radar import radar_gaussian_curve


def chunks(l, n):
    """
    Yield successive n-sized chunks from l.
    Use it to cut a big list into smaller chunks. => memory efficient
    """
    for i in range(0, len(l), n):
        yield l[i:i + n] # type: Generator[list of strings]


def production_line_manager(the_date, parameters_dict):
    """
    Here we locate the satellite files, call the comparison function
    matchproj_fun, and send the results for saving. The real deal is the
    matchproj_fun from msgr.core.msgr module.

    Parameters
    ==========
        the_date: datetime
            The day for comparison.
        parameters_dict: dict
            Dictionnary containing all parameters from the configuration file.
    """

    # Unpack the parameter dictionnaries
    PATH_params = parameters_dict['PATH_params']
    PROJ_params = parameters_dict['PROJ_params']
    RADAR_params = parameters_dict['RADAR_params']
    SAT_params = parameters_dict['SAT_params']
    SWITCH_params = parameters_dict['SWITCH_params']
    THRESHOLDS_params = parameters_dict['THRESHOLDS_params']

    l_write = SWITCH_params['l_write']
    l_gpm = SWITCH_params['l_gpm']
    outdir = PATH_params['outdir']
    satdir = PATH_params['satdir']
    rid = RADAR_params['rid']
    gr_reflectivity_offset = RADAR_params['gr_reflectivity_offset']

    # Julian day corresponding to 00 UTC
    julday = datetime.datetime(the_date.year, the_date.month, the_date.day, 0, 0, 0)
    date = julday.strftime("%Y%m%d")

    # Looking for satellite data files.
    if l_gpm:
        satfiles = glob.glob(satdir + '/*' + date + '*.HDF5')
    else:
        satfiles = glob.glob(satdir + '/*2A23*' + date + '*.HDF')
        satfiles2 = glob.glob(satdir + '/*2A25*' + date + '*.HDF')

    # Checking if found any satellite data file.
    if len(satfiles) == 0:
        print_red("\nNo satellite swaths for {}.".format(julday.strftime("%d %b %Y")))
        return None

    # Looping over satellite file
    for one_sat_file in satfiles:
        # Get orbit number
        orbit = get_orbit_number(one_sat_file)
        print_with_time("\nOrbit {} -- {}.".format(orbit, julday.strftime("%d %B %Y")))

        # Start chrono
        st_time = time.time()

        if l_gpm:
            # Calling processing function for GPM
            match_vol = processing_codes.matchproj_fun(
                PATH_params, PROJ_params, RADAR_params, SAT_params,
                SWITCH_params, THRESHOLDS_params, one_sat_file, dtime=julday)
        else:
            # It's TRIMM, need to find corresponding pair file.
            try:
                # Trying to find corresponding 2A25 TRMM file based on the orbit
                fd_25 = find_file_with_string(satfiles2, orbit)
            except IndexError:
                print_red("Found no matching 2A25 file for TRMM.")
                continue

            # Calling processing function for TRMM
            match_vol = processing_codes.matchproj_fun(PATH_params,
                PROJ_params, RADAR_params, SAT_params, SWITCH_params,
                THRESHOLDS_params, one_sat_file, fd_25, dtime=julday)

        if match_vol is None:
            continue

        print_with_time("Comparison done in %.2fs." % (time.time() - st_time))

        # Saving data
        if not l_write:
            continue

        # Output file name.
        outfilename = "RID_{}_ORBIT_{}_DATE_{}_OFFSET_{:0.2f}dB".format(rid,
            orbit, julday.strftime("%Y%m%d"), gr_reflectivity_offset)
        outfilename = os.path.join(outdir, outfilename)

        print_green("Saving data to {}, for orbit {} on {}.".format(outfilename,
            orbit, julday.strftime("%d %B %Y")), bold=True)
        save_data(outfilename, match_vol)

    return None


def main(configuration_file):
    """
    MAIN
    Multiprocessing control room
    """

    # Reading the configuration file.
    start_date, end_date, ncpu, parameters_dict = config_codes.read_configuration_file(configuration_file)

    # Generating the date range.
    date_range = pd.date_range(start_date, end_date)

    # Cutting the file list into smaller chunks. (The multiprocessing.Pool instance
    # is freed from memory, at each iteration of the main for loop).
    date_range_chunk = chunks(date_range, ncpu*2)
    for date_list in date_range_chunk:
        # Because we use multiprocessing, we need to send a list of tuple as argument of Pool.starmap.
        args_list = [None]*len(date_list)  # yes, I like declaring empty array.
        for cnt, onedate in enumerate(date_list):
            args_list[cnt] = (onedate, parameters_dict)

        # Start multiprocessing.
        with Pool(ncpu) as pool:
            pool.starmap(production_line_manager, args_list)

    return None


if __name__ == '__main__':
    parser_description = "Start MSGR - volume Matching Satellite and Ground Radar."
    parser = argparse.ArgumentParser(description=parser_description)

    parser.add_argument('-c', '--config', type=str, dest='config_file',
        help='Path to configuration file.', default=None)

    args = parser.parse_args()
    config_file = args.config_file

    if config_file is None:
        parser.error("Configuration file required.")
        sys.exit()

    main(config_file)
