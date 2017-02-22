#!python

"""
MSGR Matching Satellite and Ground Radar
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@date: 2016-12-06 (creation) 2017-2-22 (current version)
@email: valentin.louf@bom.gov.au
@company: Monash University/Bureau of Meteorology
"""

__title__ = 'matchvol'
__version__ = '0.5'
__author__ = 'Valentin Louf'
__license__ = 'MIT'
__copyright__ = 'Copyright 2017 Valentin Louf'

import glob
import time
import pyart
import pyproj  # For cartographic transformations and geodetic computations
import datetime
import warnings
import configparser
import numpy as np
import pandas as pd
from numpy import sqrt, cos, sin, pi, exp
from multiprocessing import Pool

# Custom modules
from core.util_fun import * # bunch of useful functions
from core.msgr import matchproj_fun
from core.io.save_data import save_data
from core.instruments.ground_radar import radar_gaussian_curve # functions related to the ground radar data
from core.instruments.satellite import get_orbit_number, satellite_params # functions related to the satellite data


def MAIN_matchproj_fun(the_date):
    """
    MAIN_MATCHPROJ_FUN
    Here we locate the satellite files, call the comparison function
    matchproj_fun, and send the results for saving.

    Parameters
    ==========
        the_date: datetime
            The day for comparison.
    """

    # Note the Julian day corresponding to 00 UTC
    julday = datetime.datetime(the_date.year, the_date.month, the_date.day, 0, 0, 0)
    date = julday.strftime("%Y%m%d")

    # Note the number of satellite overpasses on this day
    if l_gpm:
        satfiles = glob.glob(satdir + '/*' + date + '*.HDF5')
    else:
        satfiles = glob.glob(satdir + '/*2A23*' + date + '*.HDF')
        satfiles2 = glob.glob(satdir + '/*2A25*' + date + '*.HDF')

    if len(satfiles) == 0:
        print('')  # line break
        txt = 'No satellite swaths for ' + julday.strftime("%d %b %Y")
        print_red(txt)
        nerr[0] += 1
        return None

    for the_file in satfiles:
        orbit = get_orbit_number(the_file)

        print('')  # line break
        print_with_time("Orbit " + orbit + " -- " + julday.strftime("%d %B %Y"))

        st_time = time.time()
        if l_gpm:
            match_vol = matchproj_fun(PATH_params, PROJ_params, RADAR_params,
                                      SAT_params, SWITCH_params,
                                      THRESHOLDS_params, the_file,
                                      dtime=julday)
        else:
            try:
                # Trying to find corresponding 2A25 TRMM file based on the orbit
                fd_25 = find_file_with_string(satfiles2, orbit)
            except IndexError:
                print_red("No matching 2A25 file for TRMM.")
                continue
            # match_vol = matchproj_fun(the_file, fd_25, dtime=julday)
            match_vol = matchproj_fun(PATH_params, PROJ_params, RADAR_params,
                                      SAT_params, SWITCH_params,
                                      THRESHOLDS_params, the_file, fd_25,
                                      dtime=julday)

        end_time = time.time()
        if match_vol is None:
            continue

        print_with_time("Comparison took %.2fs." % (end_time - st_time))

        # Saving data
        if l_write:
            outfilename = "RID_" + rid + "_ORBIT_" + orbit + "_DATE_" + \
                       julday.strftime("%Y%m%d") + \
                       "_OFFSET_%1.2fdB" % (gr_reflectivity_offset)

            out_name = os.path.join(outdir, outfilename)

            txt = "Saving data to " + out_name + \
                  ". For orbit " + orbit + " on " + julday.strftime("%d %B %Y")
            print_green(txt, bold=True)
            save_data(out_name, match_vol)

    return None


def main():
    """
    MAIN
    Multiprocessing control room
    """

    date_range = pd.date_range(start_date, end_date)

    # Chunking the date_range list in order to make it smaller to ingest in
    # multiprocessing. This allows to clear multiprocessing memory at every
    # chunks and not going to cray with memory eating. It's just a little trick.
    if len(date_range) > ncpu*2:

        date_range_chunk = chunks(date_range, ncpu*2)  # Type: Generator
        for date_range_slice in date_range_chunk:
            with Pool(ncpu) as pool:
                date_list = list(date_range_slice)
                pool.map(MAIN_matchproj_fun, date_list)

    else:
        with Pool(ncpu) as pool:
            pool.map(MAIN_matchproj_fun, date_range)

    return None


if __name__ == '__main__':
    """
    GLOBAL variables declaration
    Reading configuration file.
    """

    #  Reading configuration file
    config = configparser.ConfigParser()
    config.read('config.ini')

    general = config['general']
    ncpu = general.getint('ncpu')
    date1 = general.get('start_date')
    date2 = general.get('end_date')

    switch = config['switch']
    l_write = switch.getboolean('write')   # Switch for writing out volume-matched data
    l_cband = switch.getboolean('cband')   # Switch for C-band GR
    l_dbz = switch.getboolean('dbz')       # Switch for averaging in dBZ
    l_gpm = switch.getboolean('gpm')       # Switch for GPM PR data
    l_atten = switch.getboolean('correct_gr_attenuation')       # Switch for GPM PR data

    path = config['path']
    raddir = path.get('ground_radar')
    satdir = path.get('satellite')
    outdir = path.get('output')

    GR_param = config['radar']
    radstr = GR_param.get('radar_name')
    rmin = GR_param.getfloat('rmin')  # minimum GR range (m)
    rmax = GR_param.getfloat('rmax')  # maximum GR range (m)
    rid = GR_param.get('radar_id')
    lon0 = GR_param.getfloat('longitude')
    lat0 = GR_param.getfloat('latitude')
    z0 = GR_param.getfloat('altitude')
    bwr = GR_param.getfloat('beamwidth')
    gr_reflectivity_offset = GR_param.getfloat('offset')

    thresholds = config['thresholds']
    minprof = thresholds.getint('min_profiles')  # minimum number of PR profiles with precip
    maxdt = thresholds.getfloat('max_time_delta')   # maximum PR-GR time difference (s)
    minrefg = thresholds.getfloat('min_gr_reflec')  # minimum GR reflectivity
    minrefp = thresholds.getfloat('min_sat_reflec')  # minimum PR reflectivity
    minpair = thresholds.getint('min_pair')  # minimum number of paired samples
    """ End of the section for user-defined parameters """

    start_date = datetime.datetime.strptime(date1, '%Y%m%d')
    end_date = datetime.datetime.strptime(date2, '%Y%m%d')

    # Map Projection
    # Options: projection transverse mercator, lon and lat of radar, and
    # ellipsoid WGS84
    pyproj_config = "+proj=tmerc +lon_0=%f +lat_0=%f +ellps=WGS84" % (lon0, lat0)
    smap = pyproj.Proj(pyproj_config)

    # Gaussian radius of curvatur for the radar's position
    earth_gaussian_radius = radar_gaussian_curve(lat0)

    """Stocking parameters in dictionnaries"""
    if l_gpm:
        SAT_params = satellite_params('gpm')
    else:
        SAT_params = satellite_params('trmm')

    PATH_params = dict()
    PROJ_params = dict()
    RADAR_params = dict()
    SWITCH_params = dict()
    THRESHOLDS_params = dict()

    SWITCH_params['l_cband'] = l_cband
    SWITCH_params['l_dbz'] = l_dbz
    SWITCH_params['l_gpm'] = l_gpm
    SWITCH_params['l_atten'] = l_atten

    PATH_params['raddir'] = raddir
    PATH_params['outdir'] = outdir

    RADAR_params['xmin'] = -1*rmax
    RADAR_params['xmax'] = rmax
    RADAR_params['ymin'] = -1*rmax
    RADAR_params['ymax'] = rmax
    RADAR_params['rmax'] = rmax
    RADAR_params['rmin'] = rmin
    RADAR_params['rid'] = rid
    RADAR_params['z0'] = z0
    RADAR_params['bwr'] = bwr
    RADAR_params['gr_reflectivity_offset'] = gr_reflectivity_offset

    THRESHOLDS_params['minprof'] = minprof
    THRESHOLDS_params['maxdt'] = maxdt
    THRESHOLDS_params['minrefg'] = minrefg
    THRESHOLDS_params['minrefp'] = minrefp
    THRESHOLDS_params['minpair'] = minpair

    PROJ_params['earth_gaussian_radius'] = earth_gaussian_radius
    PROJ_params['smap'] = smap

    # Printing some information about the global variables and switches
    welcome_message(l_gpm, l_atten, l_dbz, l_write, outdir, satdir, raddir,
                    ncpu, start_date, end_date)

    # Serious business starting here.
    main()