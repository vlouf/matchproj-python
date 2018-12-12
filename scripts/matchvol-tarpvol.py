"""
MSGR Matching Satellite and Ground Radar
========================================

@author: Valentin Louf
@date: 2016-12-06 (creation) 2018-04-18 (current version)
@email: valentin.louf@bom.gov.au
@company: Monash University/Bureau of Meteorology

Adapted by Joshua Soderholm for NCI odim pvol archive
"""
# Standard library import
import os
import glob
import time
import argparse
import datetime
import warnings
import traceback
import configparser
import tarfile
import shutil
import tempfile

from multiprocessing import Pool

# Others lib.
import numpy as np

# Custom lib
from msgr import cross_validation
from msgr.io.read_gpm import read_date_from_GPM
from msgr.io.read_trmm import read_date_from_TRMM
from msgr.io.save_data import save_data
from msgr.utils.misc import *


def compute_offset(matchvol_data,l_mw):
    z        = matchvol_data['z']
    zbb      = matchvol_data['zbb']
    zbbw     = matchvol_data['bbwidth']
    ref_GR   = matchvol_data['ref2']  #GR
    ref_GR_K = matchvol_data['ref5']  #GR (Kuband)
    ref_SR   = matchvol_data['ref1']  #SR
    ref_SR_S = matchvol_data['ref3']  #SR (Sband/Snow)
    std_SR   = matchvol_data['stdv1'] #SR std
    std_GR   = matchvol_data['stdv2'] #GR std
    nrej_GR  = matchvol_data['nrej1']
    nrej_SR  = matchvol_data['nrej2']
    ntot_GR  = matchvol_data['ntot1']
    ntot_SR  = matchvol_data['ntot2']
    
    if l_mw:
        #mask implemented by michael whimpey
        pos = (z > 4e3) | (ref_SR_K >= 36) | (std_SR > 4) | (std_GR > 4) | (ref_GR_K >= 36) | (ref_SR_K == 0) | (ref_GR_K < 21)
    else:
        pos   = np.logical_and(z > (zbb - zbbw/2), z < (zbb + zbbw/2)) | (ref_GR < 21) | (ref_SR_S < 21) | (nrej_GR/ntot_GR > 0.3) | (nrej_SR/ntot_SR > 0.3)
    
    ref_GR[pos] = np.NaN

    dref = (ref_GR - ref_SR_S)
    dref = dref[~np.isnan(dref)]

    if len(dref) <= 1:
        return np.NaN
    else:
        offset = np.median(dref)
        return offset


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
        sat_path = '/'.join([satdir, '**', f'*{date}*.HDF5'])
        satfiles = glob.glob(sat_path, recursive=True)
        satfiles2 = None
    else:
        satfiles = sorted(glob.glob(os.path.join(satdir, f'*2A23*{date}*.HDF')))
        satfiles2 = sorted(glob.glob(os.path.join(satdir, f'*2A25*{date}*.HDF')))

    # Checking if found any satellite data file.
    if len(satfiles) == 0:
        satfiles, satfiles2 = None, None

    return satfiles, satfiles2


def check_directory(radar_dir, satellite_dir, output_dir):
    # Check if dirs exist.
    if not os.path.isdir(radar_dir):
        raise FileNotFoundError('Ground radar data directory not found.')
    if not os.path.isdir(satellite_dir):
        raise FileNotFoundError('Satellite data directory not found.')

    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass

    return None


def multiprocessing_driver(CONFIG_FILE, ground_radar_file, one_sat_file, sat_file_2A25_trmm,
                           satellite_dtime, l_cband, l_dbz, l_gpm, l_atten, l_mw, gr_offset,
                           l_write, rid, orbit, outdir):
    """
    Buffer function that handles Exceptions while running the multiprocessing.
    Automatically runs the comparison 2 times. First with the raw radar data,
    then it computes the offset between ground radars and satellites and runs
    a second time with that offset.
    """
    datestr      = satellite_dtime.strftime('%Y%m%d')
    pass1_offset = 0
    
    for pass_number in range(2):
        try:
            # Calling processing function for TRMM
            tick = time.time()
            match_vol = cross_validation.match_volumes(CONFIG_FILE, ground_radar_file, one_sat_file,
                                                       sat_file_2A25_trmm, satellite_dtime, l_cband,
                                                       l_dbz, l_gpm, l_atten, gr_offset, rid)
        except Exception:
            traceback.print_exc()
            return None

        print_with_time("Comparison done in %.2fs." % (time.time() - tick))
        if match_vol is None:
            print_red('The comparison returned nothing.')
            return None

        delta_zh = compute_offset(match_vol, l_mw)
        # Saving data
        if l_write:
            # Output file name.
            out_path    = '/'.join([outdir, rid])
            if not os.path.exists(out_path):
                os.makedirs(out_path)
            outfilename = "RID_{}_ORBIT_{}_DATE_{}_PASS_{}.nc".format(rid, orbit, datestr, pass_number + 1)
            outfilename = os.path.join(out_path, outfilename)
            print_green("Saving data to {}.".format(outfilename), bold=True)
            if pass_number == 0:
                save_data(outfilename, match_vol, satellite_dtime, offset1=gr_offset, nb_pass=pass_number)
            else:
                save_data(outfilename, match_vol, satellite_dtime, offset1=pass1_offset, offset2=delta_zh, nb_pass=pass_number)

        #for pass 2, save the offset from pass 1
        pass1_offset = delta_zh
        #for reprocessing pass 2, use the negative of the offset as this value was inverted in compute_offset
        gr_offset    = -delta_zh
        
        if np.abs(gr_offset) < 1:
            print_green(f"No significant difference between ground radar and satellite found for {datestr}. Not doing anymore pass.")
            break
        elif np.isnan(gr_offset):
            print_red(f"Invalid offset found. Stopping comparison for {datestr}.")
            return None
        elif pass_number == 0:
            print_magenta(f"The difference between the ground radar data and the satellite data for {datestr} is of {gr_offset:0.2} dB.")

    return None


def main():
    print_with_time("Loading configuration file.")
    #  Reading configuration file
    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)

    # General info.
    general = config['general']
    date1 = general.get('start_date')
    date2 = general.get('end_date')

    # Data path
    path = config['path']
    radar_dir = path.get('ground_radar')
    satellite_dir = path.get('satellite')
    outdir = path.get('output')

    # Thresholds
    thresholds = config['thresholds']
    max_time_delta = thresholds.getfloat('max_time_delta')  # in second

    # General info about the ground radar (ID and OFFSET to apply.)
    # Radar location too.
    GR_param = config['radar']
    rid = GR_param.get('radar_id')
    radar_lat = GR_param.getfloat('latitude')
    radar_lon = GR_param.getfloat('longitude')
    rmax      = GR_param.getfloat('rmax')
    try:
        gr_offset = GR_param.getfloat('offset')
    except KeyError:
        gr_offset = 0

    # Switches.
    switch = config['switch']
    l_write = switch.getboolean('write')  # switch to save data.
    l_cband = switch.getboolean('cband')  # Switch for C-band GR
    l_dbz = switch.getboolean('dbz')  # Switch for averaging in dBZ
    l_gpm = switch.getboolean('gpm')  # Switch for GPM PR data
    l_atten = switch.getboolean('correct_gr_attenuation')
    l_mw    = switch.getboolean('mw_mask')
    # Finish reading configuration file.

    #create temp dir
    tempdir_path = tempfile.mkdtemp()
    #tempdir_path = '/home/548/jss548/dev/tmp/matchproj/g_in'
    
    
    check_directory(radar_dir, satellite_dir, outdir)

    start_date = datetime.datetime.strptime(date1, '%Y%m%d')
    end_date = datetime.datetime.strptime(date2, '%Y%m%d')
    nbdays = (end_date - start_date).days
    date_list = [start_date + datetime.timedelta(days=d) for d in range(nbdays)]

    print_yellow("Generating ground radar file list.")
    total_radar_file_list = get_files(radar_dir)
    print_yellow(f"Found {len(total_radar_file_list)} supported radar files in {radar_dir}.")

    args_list = []
    for date in date_list:
        datestr = date.strftime('%Y%m%d')

        # Extracting radar file list for this date from the total radar file list.

        radar_file_list = [f for f in total_radar_file_list if datestr in f]

        if len(radar_file_list) == 0:
            print_yellow(f"No ground radar file found for this date {datestr}")
            continue
            
        #Extract tarfile file list if requires
        if 'pvol.tar' in radar_file_list[0]:
            l_tar  = True
            tar_fn = radar_file_list[0]
            tar_id = tarfile.open(tar_fn)
            radar_file_list = tar_id.getnames()
            tar_member_list = tar_id.getmembers()
            tar_id.close()
        else:
            l_tar = False 
            
        # Looking for satellite data corresponding to this date.
        satfiles, satfiles2 = get_satfile_list(satellite_dir, datestr, l_gpm)
        if satfiles is None:
            print_red(f"No satellite data for {datestr}.")
            continue

        # if len(satfiles) > 5:
        #     print_red(f"There are more than 5 files for {datestr}. Something probably wrong in the files name. Skipping this date.")
        #     continue

        # Obtaining the satellite file(s) and reading its exact date and time.
        for idx, one_sat_file in enumerate(satfiles):
            if not l_gpm:
                try:
                    sat_file_2A25_trmm = satfiles2[idx]
                except:
                    print('no matching 2A25 file found for ', one_sat_file)
                    continue
                satellite_dtime, satellite_dist = read_date_from_TRMM(one_sat_file, radar_lat, radar_lon)
            else:
                sat_file_2A25_trmm = None
                satellite_dtime, satellite_dist = read_date_from_GPM(one_sat_file, radar_lat, radar_lon)

            orbit = get_orbit_number(one_sat_file)

            # check satellite dist
            if satellite_dist > rmax:
                continue

            # Get the datetime for each radar files
            radar_dtime = [get_time_from_filename(radfile, datestr) for radfile in radar_file_list]
            radar_dtime = list(filter(None, radar_dtime))  # Removing None values

            closest_dtime_rad = get_closest_date(radar_dtime, satellite_dtime)
            time_difference = np.abs(satellite_dtime - closest_dtime_rad)
            if time_difference.seconds > max_time_delta:
                print_red(f'Time difference is {time_difference.seconds}s while the maximum time difference allowed is {max_time_delta}s.', bold=True)
                continue

            # Radar file corresponding to the nearest scan time
            ground_radar_file = get_filename_from_date(radar_file_list, closest_dtime_rad)

            #transfer file from tar to temp dir
            if l_tar:
                #extract to temp dir
                member_idx = radar_file_list.index(ground_radar_file)
                tar_member = tar_member_list[member_idx]
                tar_id = tarfile.open(tar_fn)
                tar_id.extract(tar_member, path=tempdir_path)
                tar_id.close()
                #generate path to temp dir
                ground_radar_file = '/'.join([tempdir_path, ground_radar_file])
            
            # Argument list for multiprocessing.
            args_list.append((CONFIG_FILE, ground_radar_file, one_sat_file,
                              sat_file_2A25_trmm, satellite_dtime, l_cband,
                              l_dbz, l_gpm, l_atten, l_mw, gr_offset, l_write, rid, orbit, outdir))

    if len(args_list) == 0:
        print_red("Nothing to do. Is the configuration file correct?")
        return None

    # Start multiprocessing.
    with Pool(ncpu) as pool:
        pool.starmap(multiprocessing_driver, args_list)

    #remove tar temp path
    if l_tar:
       shutil.rmtree(tempdir_path)
    
    print('finished')
    
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

    parser.add_argument('-c', '--config', type=str, dest='config_file', help='Configuration file.', default=None, required=True)
    parser.add_argument('-j', '--cpu', dest='ncpu', default=16, type=int, help='Number of process')
    
    args = parser.parse_args()
    
    CONFIG_FILE = args.config_file
    ncpu        = args.ncpu

    warnings.simplefilter("ignore")
    main()
