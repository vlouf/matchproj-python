# Python standard library
import logging
import datetime
import warnings
import itertools

# Other libraries
import pyproj
import numpy as np
from numpy import sqrt, cos, sin, tan, pi, exp

# Custom modules
from .core import Radar, Satellite
from .io.read_radar import read_radar
from .utils import reflectivity_conversion
from .utils.misc import *


def _matching(satellite, cpol, nprof, reflectivity_satellite,
              refp_ss, refp_sh, xp, yp, zp, rt, ep, alpha, zbb, l_dbz=True):
    """
    The volume matching is done here

    Parameters:
    ===========
    satellite: Object
        Satellite object defined in core.py
    cpol: Object
        Ground radar object defined in core.py
    nprof: int
        Number of precipitating satellite rays in domain
    reflectivity_satellite: ndarray
        Satellite reflectivity.
    refp_ss: ndarray
        Satellite stratiform reflectivity converted to S/C band.
    refp_sh: ndarray
        Satellite convective reflectivity converted to S/C band.
    xp: ndarray
        x-cartesian coordinates of satellite data corrected from parallax with
        respect to the ground radar.
    yp: ndarray
        y-cartesian coordinates of satellite data corrected from parallax with
        respect to the ground radar.
    zp: ndarray
        z-cartesian coordinates of satellite data corrected from parallax with
        respect to the ground radar.
    rt: ndarray
        Approximate volume of each satellite bin
    ep: ndarray
        elev_pr_grref
    alpha: ndarray
        Angle
    zbb: float
        Bright band altitude.
    l_dbz: bool
        Are the statistics over reflectivity done in natural units or in dBZ.

    Returns:
    ========
    """
    zt = satellite.altitude
    bwt = satellite.beamwidth
    drt = satellite.dr

    xg, yg, zg = cpol.get_cartesian_coordinates()
    rg = cpol.fields['range']
    reflectivity_ground_radar = cpol.fields['reflec']
    elang = cpol.fields['elang']
    ntilt = cpol.fields['ntilt']
    dr = cpol.fields['dr']
    refg_ku = cpol.convert_refl_ku(zbb)

    bwr = cpol.beamwidth
    earth_gaussian_radius = cpol.gaussian_radius
    rmax = cpol.rmax

    try:
        reflectivity_satellite = reflectivity_satellite.filled(np.NaN)
    except AttributeError:
        pass
    try:
        reflectivity_ground_radar = reflectivity_ground_radar.filled(np.NaN)
    except AttributeError:
        pass

    # Create arrays to store comparison variables
    '''Coordinates'''
    x = np.zeros((nprof, ntilt))  # x coordinate of sample
    y = np.zeros((nprof, ntilt))  # y coordinate of sample
    z = np.zeros((nprof, ntilt))  # z coordinate of sample
    dz = np.zeros((nprof, ntilt))  # depth of sample
    ds = np.zeros((nprof, ntilt))  # width of sample
    r = np.zeros((nprof, ntilt))  # range of sample from ground radar

    '''Reflectivities'''
    ref1 = np.zeros((nprof, ntilt)) + np.NaN  # PR reflectivity
    ref2 = np.zeros((nprof, ntilt)) + np.NaN  # GR reflectivity
    ref3 = np.zeros((nprof, ntilt)) + np.NaN  # PR reflec S-band, snow
    ref4 = np.zeros((nprof, ntilt)) + np.NaN  # PR reflec S-band, hail
    ref5 = np.zeros((nprof, ntilt)) + np.NaN  # GR reflectivity Ku-band
    iref1 = np.zeros((nprof, ntilt)) + np.NaN  # path-integrated PR reflec
    iref2 = np.zeros((nprof, ntilt)) + np.NaN  # path-integrated GR reflec
    stdv1 = np.zeros((nprof, ntilt)) + np.NaN  # STD of PR reflectivity
    stdv2 = np.zeros((nprof, ntilt)) + np.NaN  # STD of GR reflectivity

    '''Number of bins in sample'''
    ntot1 = np.zeros((nprof, ntilt), dtype=int)  # Total nb of PR bin in sample
    # Nb of rejected PR bin in sample
    nrej1 = np.zeros((nprof, ntilt), dtype=int)
    ntot2 = np.zeros((nprof, ntilt), dtype=int)  # Total nb of GR bin in sample
    # Nb of rejected GR bin in sample
    nrej2 = np.zeros((nprof, ntilt), dtype=int)
    # Total volume of PR bins in sample
    vol1 = np.zeros((nprof, ntilt)) + np.NaN
    # Total volume of GR bins in sample
    vol2 = np.zeros((nprof, ntilt)) + np.NaN

    # Compute the volume of each radar (sat/ground) bin
    volp = 1.e-9 * np.pi * drt * (rt * np.pi / 180 * bwt / 2.)**2
    volg = 1e-9 * np.pi * dr * (rg * np.pi / 180 * bwr / 2)**2

    # Compute the path-integrated reflectivities at every points
    nat_refp = 10**(reflectivity_satellite / 10.0)  # In natural units
    nat_refg = 10**(reflectivity_ground_radar / 10.0)
    irefp = np.fliplr(np.nancumsum(np.fliplr(nat_refp), axis=1))
    irefg = np.nancumsum(nat_refg, axis=0)
    irefp = drt * (irefp - nat_refp / 2)
    irefg = dr * (irefg - nat_refg / 2)
    irefp = 10 * np.log10(irefp)
    irefg = 10 * np.log10(irefg)

    # Convert to linear units
    if not l_dbz:
        reflectivity_satellite = 10**(reflectivity_satellite / 10.0)
        reflectivity_ground_radar = 10**(reflectivity_ground_radar / 10.0)
        refp_ss = 10**(refp_ss / 10.0)
        refp_sh = 10**(refp_sh / 10.0)
        refg_ku = 10**(refg_ku / 10.0)

    irefp = 10**(irefp / 10.0)
    irefg = 10**(irefg / 10.0)

    # Loop over the TRMM/GPM profiles and Loop over the GR elevation scan
    for ii, jj in itertools.product(range(nprof), range(ntilt)):
        # Identify those PR bins which fall within the GR sweep
        ip = np.where((ep[ii, :] >= elang[jj] - bwr / 2) & (ep[ii, :] <= elang[jj] + bwr / 2))

        # Store the number of bins
        ntot1[ii, jj] = len(ip)
        if len(ip) == 0:
            continue

        x[ii, jj] = np.mean(xp[ii, ip])
        y[ii, jj] = np.mean(yp[ii, ip])
        z[ii, jj] = np.mean(zp[ii, ip])

        # Compute the thickness of the layer
        nip = len(ip)
        dz[ii, jj] = nip * drt * cos(pi / 180 * alpha[ii, 0])

        # Compute the PR averaging volume
        vol1[ii, jj] = np.sum(volp[ii, ip])

        # Note the mean TRMM beam diameter
        ds[ii, jj] = pi / 180 * bwt * \
            np.mean((zt - zp[ii, ip]) / cos(pi / 180 * alpha[ii, ip]))

        # Note the radar range
        s = sqrt(x[ii, jj]**2 + y[ii, jj]**2)
        r[ii, jj] = (earth_gaussian_radius + z[ii, jj]) * \
            sin(s / earth_gaussian_radius) / cos(pi / 180 * elang[jj])

        # Check that sample is within radar range
        if r[ii, jj] + ds[ii, jj] / 2 > rmax:
            continue

        # Extract the relevant PR data
        refp1 = reflectivity_satellite[ii, ip].flatten()
        refp2 = refp_ss[ii, ip].flatten()
        refp3 = refp_sh[ii, ip].flatten()
        irefp1 = irefp[ii, ip].flatten()

        # Average over those bins that exceed the reflectivity
        # threshold (linear average)

        try:
            ref1[ii, jj] = np.nanmean(refp1)
        except ValueError:
            pass
        try:
            ref3[ii, jj] = np.nanmean(refp2)
        except ValueError:
            pass
        try:
            ref4[ii, jj] = np.nanmean(refp3)
        except ValueError:
            pass
        try:
            iref1[ii, jj] = np.nanmean(irefp1)
        except ValueError:
            pass

        try:
            if not l_dbz:
                stdv1[ii, jj] = np.nanstd(10 * np.log10(refp1))
            else:
                stdv1[ii, jj] = np.nanstd(refp1)
        except ValueError:
            pass

        # Note the number of rejected bins
        nrej1[ii, jj] = int(np.sum(np.isnan(refp1)))
        if ~np.isnan(stdv1[ii, jj]) and nip - nrej1[ii, jj] > 1:
            continue

        # Compute the horizontal distance to all the GR bins
        d = sqrt((xg[:, :, jj] - x[ii, jj])**2 + (yg[:, :, jj] - y[ii, jj])**2)

        # Find all GR bins within the SR beam
        igx, igy = np.where(d <= ds[ii, jj] / 2)

        # Store the number of bins
        ntot2[ii, jj] = len(igx)
        if len(igx) == 0:
            continue

        # Extract the relevant GR data
        refg1 = reflectivity_ground_radar[:, :, jj][igx, igy].flatten()
        refg2 = refg_ku[:, :, jj][igx, igy].flatten()
        volg1 = volg[:, :, jj][igx, igy].flatten()
        irefg1 = irefg[:, :, jj][igx, igy].flatten()

        #  Comupte the GR averaging volume
        vol2[ii, jj] = np.sum(volg1)

        # Average over those bins that exceed the reflectivity
        # threshold (exponential distance and volume weighting)
        w = volg1 * exp(-1 * (d[igx, igy] / (ds[ii, jj] / 2.))**2)
        w = w * refg1 / refg2

        ref2[ii, jj] = np.nansum(w * refg1) / np.nansum(w)

        # if ref2[ii, jj] < minrefp:
        #     ref2[ii, jj] = np.NaN

        ref5[ii, jj] = np.nansum(w * refg2) / np.nansum(w)
        iref2[ii, jj] = np.nansum(w * irefg1) / np.nansum(w)

        if not l_dbz:
            stdv2[ii, jj] = np.nanstd(10 * np.log10(refg1))
        else:
            stdv2[ii, jj] = np.nanstd(refg1)

        # Note the number of rejected bins
        nrej2[ii, jj] = int(np.sum(np.isnan(refg1)))
    # END FOR (satellite profiles, radar elevation)

    # Convert back to dBZ
    iref1 = 10 * np.log10(iref1)
    iref2 = 10 * np.log10(iref2)

    if not l_dbz:
        ref1 = 10 * np.log10(ref1)
        ref2 = 10 * np.log10(ref2)
        ref3 = 10 * np.log10(ref3)
        ref4 = 10 * np.log10(ref4)
        ref5 = 10 * np.log10(ref5)

    # Correct std
    stdv1[np.isnan(stdv1)] = 0
    stdv2[np.isnan(stdv2)] = 0

    ref2[ref2 < cpol.min_refl_thrld] = np.NaN

    # Extract comparison pairs
    ipairx, ipairy = np.where((~np.isnan(ref1)) & (~np.isnan(ref2)))
    if len(ipairx) < satellite.min_pair_nb:
        print_red('Insufficient comparison pairs.')
        return None

    # Save structure
    match_vol = dict()
    match_vol['zbb'] = zbb  # Average bright band height
    match_vol['x'] = x[ipairx, ipairy]  # x coordinate of sample
    match_vol['y'] = y[ipairx, ipairy]  # y coordinate of sample
    match_vol['z'] = z[ipairx, ipairy]  # z coordinate of sample
    match_vol['dz'] = dz[ipairx, ipairy]  # depth of sample
    match_vol['ds'] = ds[ipairx, ipairy]  # width of sample
    match_vol['r'] = r[ipairx, ipairy]  # range of sample from GR
    match_vol['el'] = cpol.fields['elang'][ipairy]  # Elevation angle

    match_vol['ref1'] = ref1[ipairx, ipairy]  # PR reflectivity (Ku-band)
    match_vol['ref2'] = ref2[ipairx, ipairy]  # GR reflectivity (S/C-band)
    match_vol['ref3'] = ref3[ipairx, ipairy]  # PR reflectivity (S-band stratiform)
    match_vol['ref4'] = ref4[ipairx, ipairy]  # PR reflectivity (S-band convective)
    match_vol['ref5'] = ref5[ipairx, ipairy]  # GR reflectivity (Ku-band)
    match_vol['iref1'] = iref1[ipairx, ipairy]  # path-integrated PR reflectivity
    match_vol['iref2'] = iref2[ipairx, ipairy]  # path-integrated GR reflectivity
    match_vol['ntot1'] = ntot1[ipairx, ipairy]  # total number of PR bins in sample
    match_vol['nrej1'] = nrej1[ipairx, ipairy]  # number of rejected PR bins in sample
    match_vol['ntot2'] = ntot2[ipairx, ipairy]  # total number of GR bins in sample
    match_vol['nrej2'] = nrej2[ipairx, ipairy]  # number of rejected GR bins in sample

    match_vol['stdv1'] = stdv1[ipairx, ipairy]  # std. dev. of PR reflectivity in sample
    match_vol['stdv2'] = stdv2[ipairx, ipairy]  # std. dev. of GR reflectivity in sample
    match_vol['vol1'] = vol1[ipairx, ipairy]  # total volume of PR bins in sample
    match_vol['vol2'] = vol2[ipairx, ipairy]  # total volume of GR bins in sample

    return match_vol, ipairx


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
                  l_cband=True, l_dbz=True, l_atten=True, gr_offset=0):
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
    l_cband, l_dbz, l_atten: bool
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

    # Projecting on a WGS84 grid.
    pyproj_config = "+proj=tmerc +lon_0=%f +lat_0=%f +ellps=WGS84" % (cpol.longitude, cpol.latitude)
    smap = pyproj.Proj(pyproj_config)

    # Convert to Cartesian coordinates
    satellite_proj_cart = smap(satellite.lon, satellite.lat)
    xproj_sat = satellite_proj_cart[0]
    yproj_sat = satellite_proj_cart[1]

    # Identify profiles withing the domnain
    ioverx, iovery = np.where((xproj_sat >= cpol.xmin) & (xproj_sat <= cpol.xmax) &
                              (yproj_sat >= cpol.ymin) & (yproj_sat <= cpol.ymax))
    if len(ioverx) == 0:
        print_red("Insufficient satellite rays in domain for " + dtime_sat.strftime("%d %b %Y"))
        logging.error("Insufficient satellite rays in domain for " + dtime_sat.strftime("%d %b %Y"))
        return None

    # Note the first and last scan indices
    # i1x, i1y = np.min(ioverx), np.min(iovery)
    # i2x, i2y = np.max(ioverx), np.max(iovery)

    # Determine the datetime of the closest approach of TRMM to the GR
    xclose_sat = xproj_sat[:, 24]  # Grid center
    yclose_sat = yproj_sat[:, 24]
    # iclose = np.argmin(sqrt(xclose_sat**2 + yclose_sat**2))

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
    # range_pr_grref = (cpol.gaussian_radius + z_sat_pxcorr) * sin(gamma) / cos(pi / 180 * elev_pr_grref)  # Not used
    # azi_pr_grref = 90 - 180 / pi * np.arctan2(yproj_sat_pxcorr, xproj_sat_pxcorr)  # Not used

    # Determine the median brightband height
    ibb = np.where((zbb > 0) & (bbwidth > 0) & (quality == 1))[0]
    nbb = len(ibb)
    if nbb >= satellite.min_prof_nb:
        zbb = np.median(zbb[ibb])
        bbwidth = np.median(bbwidth[ibb])
    else:
        print_red('Insufficient bright band rays %i for ' % (nbb) + dtime_sat.strftime("%d %b %Y"))
        logging.error('Insufficient bright band rays %i for ' % (nbb) + dtime_sat.strftime("%d %b %Y"))
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
    radar = read_radar(radfile, offset=cpol.offset)
    if radar is None:
        print_red("Could not read the ground radar file. Doing nothing.")
        return None
    dtime_radar = radar['time']
    time_difference = np.abs(dtime_sat - dtime_radar)

    cpol.set_fields(radar)
    print_yellow("Ground radar data loaded.")

    # The Call.
    print_magenta("Starting volume matching.")
    match_vol, ipairx = _matching(satellite, cpol, nprof, dbz_sat, refp_ss, refp_sh, xproj_sat_pxcorr,
                                  yproj_sat_pxcorr, z_sat_pxcorr, rt, elev_pr_grref, alpha_pxcorr, zbb, l_dbz)
    if match_vol is None:
        return None

    print_magenta("Volume matching done.")

    match_vol['date'] = dtime_sat
    match_vol['bbwidth'] = bbwidth
    match_vol['dt'] = time_difference.seconds

    match_vol['sfc'] = sfc[ipairx]
    match_vol['ptype'] = ptype[ipairx]
    match_vol['iray'] = iray[ipairx]
    match_vol['iscan'] = iscan[ipairx]

    return match_vol
