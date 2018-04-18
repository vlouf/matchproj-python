# Serious business is done here.
import logging
import datetime
import warnings

import pyproj
import numpy as np
from numpy import sqrt, cos, sin, tan, pi, exp

# Custom modules
from . import volume_matching
from .core import Radar, Satellite
from .io.read_radar import read_radar
from .utils import reflectivity_conversion
from .utils.misc import *  # bunch of useful functions


def correct_parallax(xc, yc, xp, yp, alpha, the_range):
    """
    Correct for parallax to get x, y, z coordinates.

    alpha dim is nprof x 1 and now we want nprof x nbin
    xc, yc, xp, yp dimensions are nprof x 1
    """

    nprof, nbin = the_range.shape
    alpha_corr = np.zeros((nprof, nbin))
    xc0 = np.zeros((nprof, nbin))
    yc0 = np.zeros((nprof, nbin))
    xp0 = np.zeros((nprof, nbin))
    yp0 = np.zeros((nprof, nbin))
    for idx in range(0, nbin):
        alpha_corr[:, idx] = alpha[:]
        xc0[:, idx] = xc[:]
        yc0[:, idx] = yc[:]
        xp0[:, idx] = xp[:]
        yp0[:, idx] = yp[:]

    alpha = alpha_corr
    zp = the_range * cos(pi / 180. * alpha_corr)
    ds = the_range * sin(pi / 180. * alpha_corr)
    ang = np.arctan2(yp0 - yc0, xp0 - xc0)
    dx = ds * cos(ang)
    dy = ds * sin(ang)
    xp = xp0 + dx
    yp = yp0 + dy

    return xp, yp, zp, ds, alpha


def match_volumes(configuration_file, radfile, sat_file_1, sat_file_2A25_trmm=None, dtime_sat=None,
                  l_cband=True, l_dbz=True, l_gpm=True, l_atten=True, gr_offset=0):
    '''
    MATCHPROJ_FUN

    Parameters
    ==========
    configuration_file: str
        Configuration file.
    radfile: st
        Ground radar file corresponding to satellite pass.
    sat_file_1: str.
        GPM file or TRMM 2A23 file.
    sat_file_2A25_trmm: str
        TRMM 2A25 files (None for GPM).
    dtime_sat: str
        Date of current processing.
    l_cband, l_dbz, l_gpm, l_atten: bool
        Switches for C-band, use of natural reflectivity, is this GPM, and
        attenuation correction

    Returns
    =======
    match_vol: dict
        A dictionnary structure containing the comparable reflectivities.
    '''
    logging.basicConfig(filename="log_matchvol_{}.log".format(dtime_sat.strftime("%Y%m%d")), level=logging.DEBUG)
    # Spawning Radar and Satellite
    cpol = Radar(configuration_file, gr_offset=gr_offset)
    satellite = Satellite(configuration_file, sat_file_1, sat_file_2A25_trmm)

    date = dtime_sat.strftime("%Y%m%d")

    # Projecting on a WGS84 grid.
    pyproj_config = "+proj=tmerc +lon_0=%f +lat_0=%f +ellps=WGS84" % (cpol.longitude, cpol.latitude)
    smap = pyproj.Proj(pyproj_config)

    day_of_treatment = dtime_sat

    # Convert to Cartesian coordinates
    satellite_proj_cart = smap(satellite.lon, satellite.lat)
    xproj_sat = satellite_proj_cart[0]
    yproj_sat = satellite_proj_cart[1]

    # Identify profiles withing the domnain
    ioverx, iovery = np.where((xproj_sat >= cpol.xmin) & (xproj_sat <= cpol.xmax) & (yproj_sat >= cpol.ymin) & (yproj_sat <= cpol.ymax))
    if len(ioverx) == 0:
        print_red("Insufficient satellite rays in domain for " + dtime_sat.strftime("%d %b %Y"))
        logging.error("Insufficient satellite rays in domain for " + dtime_sat.strftime("%d %b %Y"))
        return None

    # Note the first and last scan indices
    i1x, i1y = np.min(ioverx), np.min(iovery)
    i2x, i2y = np.max(ioverx), np.max(iovery)

    # Determine the datetime of the closest approach of TRMM to the GR
    xclose_sat = xproj_sat[:, 24]  # Grid center
    yclose_sat = yproj_sat[:, 24]
    iclose = np.argmin(sqrt(xclose_sat**2 + yclose_sat**2))

    # Compute the distance of every ray to the radar
    dist_to_gr_rays = sqrt(xproj_sat**2 + yproj_sat**2)

    # Identify precipitating profiles within the radaar range limits
    iscan, iray = np.where((dist_to_gr_rays >= cpol.rmin) & (dist_to_gr_rays <= cpol.rmax) & (satellite.pflag == 2))
    nprof = len(iscan)
    if nprof < satellite.min_prof_nb:
        print_red('Insufficient precipitating satellite rays in domain %i.' % (nprof))
        logging.error('Insufficient precipitating satellite rays in domain %i.' % (nprof))
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
        logging.error('Insufficient bright band rays %i for ' % (nbb) + day_of_treatment.strftime("%d %b %Y"))
        return None

    print_green("Satellite side OK.")

    # Set all values less than satellite.min_refl_thrld as missing
    dbz_sat = np.ma.masked_where(dbz_sat < satellite.min_refl_thrld, dbz_sat)

    # Convert to S-band using method of Cao et al. (2013)
    if l_cband:
        refp_ss, refp_sh = reflectivity_conversion.convert_to_Cband(dbz_sat, z_sat_pxcorr, zbb, bbwidth)
    else:
        refp_ss, refp_sh = reflectivity_conversion.convert_to_Sband(dbz_sat, z_sat_pxcorr, zbb, bbwidth)

    print_yellow("Reading {}.".format(radfile))
    radar = read_radar(radfile, l_atten, cpol.offset)
    if radar is None:
        print_red("Could not read the ground radar file. Doing nothing.")
        return None
    dtime_radar = radar['time']
    time_difference = np.abs(dtime_sat - dtime_radar)

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
        print_red('Insufficient comparison pairs for ' + day_of_treatment.strftime("%d %b %Y"))
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
