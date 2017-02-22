#!python

"""
MSGR Matching Satellite and Ground Radar
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@date: 2016-12-06 (creation) 2017-2-22 (current version)
@email: valentin.louf@bom.gov.au
@company: Monash University/Bureau of Meteorology
"""

import os
import sys
import glob
import msgr
import time
# import pyart
import pyproj  # For cartographic transformations and geodetic computations
import datetime
import configparser
import pandas as pd
from multiprocessing import Pool

# Custom modules
from msgr.core.parser import parse
from msgr.core.util_fun import * # bunch of useful functions
from msgr.core.msgr import matchproj_fun
from msgr.core.io.save_data import save_data
from msgr.core.instruments.ground_radar import radar_gaussian_curve # functions related to the ground radar data
from msgr.core.instruments.satellite import get_orbit_number, satellite_params # functions related to the satellite data


def read_configuration_file(config_file):
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

    welcome_message(l_gpm, l_atten, l_dbz, l_write, outdir, satdir, raddir,
                    ncpu, start_date, end_date)

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
    SWITCH_params['l_write'] = l_write
    SWITCH_params['l_atten'] = l_atten

    PATH_params['raddir'] = raddir
    PATH_params['satdir'] = satdir
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

    PARAMETERS_dict = dict()

    PARAMETERS_dict['PATH_params'] = PATH_params
    PARAMETERS_dict['PROJ_params'] = PROJ_params
    PARAMETERS_dict['RADAR_params'] = RADAR_params
    PARAMETERS_dict['SAT_params'] = SAT_params
    PARAMETERS_dict['SWITCH_params'] = SWITCH_params
    PARAMETERS_dict['THRESHOLDS_params'] = THRESHOLDS_params

    return start_date, end_date, ncpu, PARAMETERS_dict


def MAIN_matchproj_fun(kwarg):
    """
    MAIN_MATCHPROJ_FUN
    Here we locate the satellite files, call the comparison function
    matchproj_fun, and send the results for saving. The real deal is the
    matchproj_fun from msgr.core.msgr module.

    Parameters
    ==========
        kwarg: tuple containing:
            the_date: datetime
                The day for comparison.
            PARAMETERS_dict: dict
                Dictionnary containing all parameters from the configuration
                file.
    """

    the_date, PARAMETERS_dict = kwarg

    # Unpack the parameter dictionnaries
    PATH_params = PARAMETERS_dict['PATH_params']
    PROJ_params = PARAMETERS_dict['PROJ_params']
    RADAR_params = PARAMETERS_dict['RADAR_params']
    SAT_params = PARAMETERS_dict['SAT_params']
    SWITCH_params = PARAMETERS_dict['SWITCH_params']
    THRESHOLDS_params = PARAMETERS_dict['THRESHOLDS_params']

    l_write = SWITCH_params['l_write']
    l_gpm = SWITCH_params['l_gpm']
    outdir = PATH_params['outdir']
    satdir = PATH_params['satdir']
    rid = RADAR_params['rid']
    gr_reflectivity_offset = RADAR_params['gr_reflectivity_offset']

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

        print_with_time("Comparison done in %.2fs." % (end_time - st_time))

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


def main(argv):
    """
    MAIN
    Multiprocessing control room
    """

    # Parsing argv to get the configuration file.
    configuration_file = parse(argv)

    if configuration_file is None:
        print_red("No configuration file given. Type `generate_config_matchvol" + \
                  " -e` for generating a configuration file example.")
        sys.exit()

    # Reading the configuration file.
    start_date, end_date, ncpu, PARAMETERS_dict = read_configuration_file(configuration_file)

    # Printing some information about the global variables and switches
    date_range = pd.date_range(start_date, end_date)

    # Chunking the date_range list in order to make it smaller to ingest in
    # multiprocessing. This allows to clear multiprocessing memory at every
    # chunks and not going to cray with memory eating. It's just a little trick.
    if len(date_range) > ncpu*2:

        date_range_chunk = chunks(date_range, ncpu*2)  # Type: Generator
        for date_range_slice in date_range_chunk:
            with Pool(ncpu) as pool:
                date_list = list(date_range_slice)
                kwarg = []
                for dl in date_list:
                    kwarg.append((dl, PARAMETERS_dict))

                pool.map(MAIN_matchproj_fun, kwarg)

    else:
        with Pool(ncpu) as pool:
            kwarg = []
            for dl in date_range:
                kwarg.append((dl, PARAMETERS_dict))

            pool.map(MAIN_matchproj_fun, kwarg)

    return None


if __name__ == '__main__':
    main(sys.argv)
