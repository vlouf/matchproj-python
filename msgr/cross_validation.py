# Serious business is done here.
import datetime
import warnings

import pyproj
import numpy as np
from numpy import sqrt, cos, sin, pi, exp

# Custom modules
from . import volume_matching
from .core import Radar, Satellite
from .io.read_radar import read_radar
from .utils import reflectivity_conversion
from .utils.misc import *  # bunch of useful functions
from .instruments.satellite import correct_parallax


def match_volumes(configuration_file, radar_file_list, sat_file_1, sat_file_2A25_trmm=None, dtime=None,
                  l_cband=True, l_dbz=True, l_gpm=True, l_atten=True):
    '''
    MATCHPROJ_FUN

    Parameters
    ==========
    configuration_file: str
        Configuration file.
    radar_file_list: List
        List of radar files for the current date.
    sat_file_1: str.
        GPM file or TRMM 2A23 file.
    sat_file_2A25_trmm: str
        TRMM 2A25 files (None for GPM).
    dtime: str
        Date of current processing.
    l_cband, l_dbz, l_gpm, l_atten: bool
        Switches for C-band, use of natural reflectivity, is this GPM, and
        attenuation correction

    Returns
    =======
    match_vol: dict
        A dictionnary structure containing the comparable reflectivities.
    '''
    # Spawning Radar and Satellite
    cpol = Radar(configuration_file)
    satellite = Satellite(configuration_file, sat_file_1, sat_file_2A25_trmm)

    # Projecting on a WGS84 grid.
    pyproj_config = "+proj=tmerc +lon_0=%f +lat_0=%f +ellps=WGS84" % (cpol.longitude, cpol.latitude)
    smap = pyproj.Proj(pyproj_config)

    day_of_treatment = dtime

    # Convert to Cartesian coordinates
    satellite_proj_cart = smap(satellite.lon, satellite.lat)
    xproj_sat = satellite_proj_cart[0]
    yproj_sat = satellite_proj_cart[1]

    # Identify profiles withing the domnain
    ioverx, iovery = np.where((xproj_sat >= cpol.xmin) & (xproj_sat <= cpol.xmax) & (yproj_sat >= cpol.ymin) & (yproj_sat <= cpol.ymax))
    if len(ioverx) == 0:
        print_red("Insufficient satellite rays in domain for " + day_of_treatment.strftime("%d %b %Y"))
        return None

    # Note the first and last scan indices
    i1x, i1y = np.min(ioverx), np.min(iovery)
    i2x, i2y = np.max(ioverx), np.max(iovery)

    # Determine the datetime of the closest approach of TRMM to the GR
    xclose_sat = xproj_sat[:, 24]  # Grid center
    yclose_sat = yproj_sat[:, 24]
    iclose = np.argmin(sqrt(xclose_sat**2 + yclose_sat**2))

    dtime_sat = satellite.datetime_for_bin(iclose)
    date = dtime_sat.strftime("%Y%m%d")

    # Compute the distance of every ray to the radar
    dist_to_gr_rays = sqrt(xproj_sat**2 + yproj_sat**2)

    # Identify precipitating profiles within the radaar range limits
    iscan, iray = np.where((dist_to_gr_rays >= cpol.rmin) &
                           (dist_to_gr_rays <= cpol.rmax) &
                           (satellite.pflag == 2))
    nprof = len(iscan)
    if nprof < satellite.min_prof_nb:
        print_red('Insufficient precipitating satellite rays in domain %i.' % (nprof))
        return None

    # Extract data for these rays
    xproj_sat = xproj_sat[iscan, iray]
    yproj_sat = yproj_sat[iscan, iray]
    xclose_sat = xclose_sat[iscan]
    yclose_sat = yclose_sat[iscan]
    ptype = satellite.ptype[iscan, iray]
    zbb = satellite.zbb[iscan, iray]
    bbwidth = satellite.bbwidth[iscan, iray]
    sfc = satellite.sfc[iscan, iray]
    quality = satellite.quality[iscan, iray]
    alpha = satellite.alpha[iray]

    tmp = np.zeros((nprof, satellite.nbin), dtype=float)
    for k in range(0, satellite.nbin):
        tmp[:, k] = (satellite.refl[:, :, k])[iscan, iray]
    dbz_sat = tmp

    # we want a shape of (nprof, satellite.nbin)
    range_sat_2d = np.zeros((nprof, satellite.nbin))
    for idx in range(0, nprof):
        range_sat_2d[idx, :] = satellite.dr * np.arange(satellite.nbin)

    # Correct coordinates from parallax errors.
    rslt = correct_parallax(xclose_sat, yclose_sat, xproj_sat, yproj_sat, alpha, range_sat_2d)
    xproj_sat_pxcorr, yproj_sat_pxcorr, z_sat_pxcorr, ds_pxcorr, alpha_pxcorr = rslt

    if len(ds_pxcorr) == 0:
        return None
    if np.min(ds_pxcorr) < 0:
        return None

    # Compute the (approximate) volume of each PR bin
    rt = satellite.altitude / cos(pi / 180 * alpha_pxcorr) - range_sat_2d

    # Compute the ground-radar coordinates of the PR pixels
    gamma = sqrt(xproj_sat_pxcorr**2 + yproj_sat_pxcorr**2) / cpol.gaussian_radius
    elev_pr_grref = 180 / pi * np.arctan((cos(gamma) - (cpol.gaussian_radius + cpol.altitude) / (cpol.gaussian_radius + z_sat_pxcorr)) / sin(gamma))
    range_pr_grref = (cpol.gaussian_radius + z_sat_pxcorr) * sin(gamma) / cos(pi / 180 * elev_pr_grref)  # Not used
    azi_pr_grref = 90 - 180 / pi * np.arctan2(yproj_sat_pxcorr, xproj_sat_pxcorr)  # Not used

    # Determine the median brightband height
    ibb = np.where((zbb > 0) & (bbwidth > 0) & (quality == 1))[0]
    nbb = len(ibb)
    if nbb >= satellite.min_prof_nb:
        zbb = np.median(zbb[ibb])
        bbwidth = np.median(bbwidth[ibb])
    else:
        print_red('Insufficient bright band rays %i for ' % (nbb) + day_of_treatment.strftime("%d %b %Y"))
        return None

    if radar_file_list is None:
        print_green("Satellite side OK for this date {}. You can fetch the data for it.".format(dtime_sat.strftime("%Y%m%d_%H%M")))
        return None

    print_green("Satellite side OK. Looking at the ground radar data now.")

    # Set all values less than satellite.min_refl_thrld as missing
    dbz_sat = np.ma.masked_where(dbz_sat < satellite.min_refl_thrld, dbz_sat)

    # Convert to S-band using method of Cao et al. (2013)
    if l_cband:
        refp_ss, refp_sh = reflectivity_conversion.convert_to_Cband(dbz_sat, z_sat_pxcorr, zbb, bbwidth)
    else:
        refp_ss, refp_sh = reflectivity_conversion.convert_to_Sband(dbz_sat, z_sat_pxcorr, zbb, bbwidth)

    # Get the datetime for each radar files
    dtime_radar = [get_time_from_filename(radfile, date) for radfile in radar_file_list]
    dtime_radar = list(filter(None, dtime_radar))  # Removing None values

    if len(dtime_radar) == 0:
        print_red("No corresponding ground radar files for {}.".format(date))
        return None

    # Find the nearest scan time
    # Looking at the time difference between satellite and radar
    closest_dtime_rad = get_closest_date(dtime_radar, dtime_sat)
    time_difference = np.abs(dtime_sat - closest_dtime_rad)
    if time_difference.seconds > satellite.max_time_delta:
        print_red('Time difference is %is while the maximum accpetable value is %is.' %
                  (time_difference.seconds, satellite.max_time_delta), bold=True)
        return None

    # Radar file corresponding to the nearest scan time
    radfile = get_filename_from_date(radar_file_list, closest_dtime_rad)
    time = closest_dtime_rad  # Keeping the IDL program notation

    print_yellow("Reading {}.".format(radfile))
    radar = read_radar(radfile, l_atten, cpol.offset)
    cpol.set_fields(radar)
    print_yellow("Ground radar data loaded.")

    # The Call.
    print_magenta("Starting volume matching.")
    rslt = volume_matching.process(satellite, cpol, nprof, dbz_sat,
                                   refp_ss, refp_sh, xproj_sat_pxcorr, yproj_sat_pxcorr, z_sat_pxcorr,
                                   rt, elev_pr_grref, alpha_pxcorr, zbb, l_dbz)
    x, y, z, dz, ds, r, ref1, ref2, ref3, ref4, ref5, iref1, iref2, stdv1, stdv2, ntot1, nrej1, ntot2, nrej2, vol1, vol2 = rslt
    print_magenta("Volume matching done.")

    # Correct std
    stdv1[np.isnan(stdv1)] = 0
    stdv2[np.isnan(stdv2)] = 0

    ref2[ref2 < cpol.min_refl_thrld] = np.NaN

    # Extract comparison pairs
    ipairx, ipairy = np.where((~np.isnan(ref1)) & (~np.isnan(ref2)))
    if len(ipairx) < satellite.min_pair_nb:
        print_red(
            'Insufficient comparison pairs for ' +
            day_of_treatment.strftime("%d %b %Y"))
        return None

    # Save structure
    match_vol = dict()

    match_vol['zbb'] = zbb
    match_vol['date'] = day_of_treatment
    match_vol['bbwidth'] = bbwidth
    match_vol['dt'] = time_difference.seconds

    match_vol['x'] = x[ipairx, ipairy]
    match_vol['y'] = y[ipairx, ipairy]
    match_vol['z'] = z[ipairx, ipairy]
    match_vol['dz'] = dz[ipairx, ipairy]
    match_vol['ds'] = ds[ipairx, ipairy]
    match_vol['r'] = r[ipairx, ipairy]
    match_vol['el'] = cpol.fields['elang'][ipairy]

    match_vol['ref1'] = ref1[ipairx, ipairy]
    match_vol['ref2'] = ref2[ipairx, ipairy]
    match_vol['ref3'] = ref3[ipairx, ipairy]
    match_vol['ref4'] = ref4[ipairx, ipairy]
    match_vol['ref5'] = ref5[ipairx, ipairy]
    match_vol['iref1'] = iref1[ipairx, ipairy]
    match_vol['iref2'] = iref2[ipairx, ipairy]
    match_vol['ntot1'] = ntot1[ipairx, ipairy]
    match_vol['nrej1'] = nrej1[ipairx, ipairy]
    match_vol['ntot2'] = ntot2[ipairx, ipairy]
    match_vol['nrej2'] = nrej2[ipairx, ipairy]

    match_vol['sfc'] = sfc[ipairx]
    match_vol['ptype'] = ptype[ipairx]
    match_vol['iray'] = iray[ipairx]
    match_vol['iscan'] = iscan[ipairx]

    match_vol['stdv1'] = stdv1[ipairx, ipairy]
    match_vol['stdv2'] = stdv2[ipairx, ipairy]
    match_vol['vol1'] = vol1[ipairx, ipairy]
    match_vol['vol2'] = vol2[ipairx, ipairy]

    return match_vol
