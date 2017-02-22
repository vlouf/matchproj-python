# Serious business is done here.

import datetime
import warnings
import numpy as np
import itertools
from numpy import sqrt, cos, sin, pi, exp

# Custom modules
from . import reflectivity_conversion

from .util_fun import * # bunch of useful functions
from .io.read_gpm import read_gpm
from .io.read_trmm import read_trmm
from .io.read_radar import read_radar
from .instruments.satellite import correct_parallax    # functions related to the satellite data


def matchproj_fun(PATH_params, PROJ_params, RADAR_params, SAT_params,
                  SWITCH_params, THRESHOLDS_params, sat_file_1,
                  sat_file_2A25_trmm=None, dtime=None):
    '''
    MATCHPROJ_FUN

    Parameters
    ==========
        PATH_params: dict
            Dictionnary containing paths information.
        PROJ_params: dict
            Dictionnary containing map projection information.
        RADAR_params: dict
            Dictionnary containing the ground radar parameters.
        SAT_params: dict
            Dictionnary containing the satellite radar parameters.
        SWITCH_params: dict
            Dictionnary containing the information about different switches,
            i.e. DBZ statistics or natural units, attenuation correction,
            ground radar is C-Band, or satellite is GPM.
        THRESHOLDS_params: dict
            Dictionnary containing

        sat_file_1:
            Satellite data filename (just one file for GPM, or 2A23 for TRMM).
        sat_file_2A25_trmm:
            Second satellite file (TRMM only).
        dtime: dateime
            Date of the day of comparison.

    Returns
    =======
        match_vol: dict
            A dictionnary structure containing the comparable reflectivities.
    '''

    raddir = PATH_params['raddir']

    earth_gaussian_radius = PROJ_params['earth_gaussian_radius']
    smap = PROJ_params['smap']

    xmin = RADAR_params['xmin']
    xmax = RADAR_params['xmax']
    ymin = RADAR_params['ymin']
    ymax = RADAR_params['ymax']
    rmin = RADAR_params['rmin']
    rmax = RADAR_params['rmax']
    z0 = RADAR_params['z0']
    bwr = RADAR_params['bwr']
    gr_reflectivity_offset = RADAR_params['gr_reflectivity_offset']

    zt = SAT_params['zt']
    drt = SAT_params['drt']
    bwt = SAT_params['bwt']

    l_cband = SWITCH_params['l_cband']
    l_dbz = SWITCH_params['l_dbz']
    l_gpm = SWITCH_params['l_gpm']
    l_atten = SWITCH_params['l_atten']

    minprof = THRESHOLDS_params['minprof']
    maxdt = THRESHOLDS_params['maxdt']
    minrefg = THRESHOLDS_params['minrefg']
    minrefp = THRESHOLDS_params['minrefp']
    minpair = THRESHOLDS_params['minpair']

    julday = dtime
    if l_gpm:
        sat = read_gpm(sat_file_1)
        txt = 'READING ' + sat_file_1
        print_with_time(txt)
    else:
        sat = read_trmm(sat_file_1, sat_file_2A25_trmm)
        print_with_time("READING " + sat_file_1)
        print_with_time("READING " + sat_file_2A25_trmm)

    if sat is None:
        print_red('Bad satellite data')
        return None

    # Extracting lat/lon just to make a check
    lonp = sat['lon']
    latp = sat['lat']

    # Convert to Cartesian coordinates
    res = smap(lonp, latp)
    xp = res[0]
    yp = res[1]

    # Identify profiles withing the domnain
    ioverx, iovery = np.where((xp >= xmin) & (xp <= xmax) &
                              (yp >= ymin) & (yp <= ymax))

    if len(ioverx) == 0:
        print_red("Insufficient satellite rays in domain for " + julday.strftime("%d %b %Y"))
        return None

    # Extracting the rest of the satellite data
    nscan = sat['nscan']
    nray = sat['nray']
    nbin = sat['nbin']
    yearp = sat['year']
    monthp = sat['month']
    dayp = sat['day']
    hourp = sat['hour']
    minutep = sat['minute']
    secondp = sat['second']
    pflag = sat['pflag']
    ptype = sat['ptype']
    zbb = sat['zbb']
    bbwidth = sat['bbwidth']
    sfc = sat['sfc']
    quality = sat['quality']
    reflectivity_satellite = sat['refl']

    # Note the first and last scan indices
    i1x, i1y = np.min(ioverx), np.min(iovery)
    i2x, i2y = np.max(ioverx), np.max(iovery)

    # Identify the coordinates of these points
    xf = xp[i1x:i2x]
    yf = yp[i1y:i2y]

    # Determine the date and time (in seconds since the start of the day)
    # of the closest approach of TRMM to the GR
    xc = xp[:, 24]  # Grid center
    yc = yp[:, 24]
    dc = sqrt(xc**2 + yc**2)
    iclose = np.argmin(dc)

    year = yearp[iclose]
    month = monthp[iclose]
    day = dayp[iclose]
    hour = hourp[iclose]
    minute = minutep[iclose]
    second = secondp[iclose]

    date = "%i%02i%02i" % (year, month, day)
    timep = "%02i%02i%02i" % (hour, minute, second)
    dtime_sat = datetime.datetime(year, month, day, hour, minute, second)
    # dtime_sat corresponds to the julp/tp stuff in the IDL code

    # Compute the distance of every ray to the radar
    d = sqrt(xp**2 + yp**2)

    # Identify precipitating profiles within the radaar range limits
    iscan, iray = np.where((d >= rmin) & (d <= rmax) & (pflag == 2))
    nprof = len(iscan)
    if nprof < minprof:
        print_red('Insufficient precipitating satellite rays in domain %i.' % (nprof))
        return None

    # Extract data for these rays
    xp = xp[iscan, iray]
    yp = yp[iscan, iray]
    xc = xc[iscan]
    yc = yc[iscan]
    ptype = ptype[iscan, iray]
    zbb = zbb[iscan, iray]
    bbwidth = bbwidth[iscan, iray]
    sfc = sfc[iscan, iray]
    quality = quality[iscan, iray]

    tmp = np.zeros((nprof, nbin), dtype=float)
    for k in range(0, nbin):
        tmp[:, k] = (reflectivity_satellite[:, :, k])[iscan, iray]

    reflectivity_satellite = tmp

    # Note the scan angle for each ray
    alpha = np.abs(-17.04 + np.arange(nray)*0.71)
    alpha = alpha[iray]

    # the_range shape is (nbin, ), and we now wnat to copy it for (nprof, nbin)
    the_range_1d = np.arange(nbin)*drt
    the_range = np.zeros((nprof, nbin))
    for idx in range(0, nprof):
        the_range[idx, :] = the_range_1d[:]

    xp, yp, zp, ds, the_alpha = correct_parallax(xc, yc, xp, yp, alpha, the_range)
    alpha = the_alpha

    if len(ds) == 0:
        return None
    if np.min(ds) < 0:
        return None

    # Compute the (approximate) volume of each PR bin
    rt = zt/cos(pi/180*alpha) - the_range
    volp = drt*(1.e-9)*pi*(rt*pi/180*bwt/2.)**2

    # Compute the ground-radar coordinates of the PR pixels
    sp = sqrt(xp**2 + yp**2)
    gamma = sp/earth_gaussian_radius
    ep = 180/pi*np.arctan((cos(gamma) - \
         (earth_gaussian_radius + z0)/(earth_gaussian_radius + zp))/sin(gamma))
    # rp = (earth_gaussian_radius + zp)*sin(gamma)/cos(pi/180*ep)  # Not used
    # ap = 90-180/pi*np.arctan2(yp, xp)  # Shape (nprof x nbin)  # Not used

    # Determine the median brightband height
    # 1D arrays
    ibb = np.where((zbb > 0) & (bbwidth > 0) & (quality == 1))[0]
    nbb = len(ibb)
    if nbb >= minprof:
        zbb = np.median(zbb[ibb])
        bbwidth = np.median(bbwidth[ibb])
    else:
        print_red('Insufficient bright band rays %i for ' % (nbb) +
                  julday.strftime("%d %b %Y"))
        return None

    # Set all values less than minrefp as missing
    # ibadx, ibady = np.where(reflectivity_satellite < minrefp)
    # if len(ibadx) > 0:
    #     reflectivity_satellite[ibadx, ibady] = np.NaN
    reflectivity_satellite[reflectivity_satellite < minrefp] = np.NaN

    # Convert to S-band using method of Cao et al. (2013)
    if l_cband:
        refp_ss, refp_sh = reflectivity_conversion.convert_to_Cband(reflectivity_satellite, zp, zbb, bbwidth)
    else:
        refp_ss, refp_sh = reflectivity_conversion.convert_to_Sband(reflectivity_satellite, zp, zbb, bbwidth)

    # Get the ground radar file lists (next 20 lines can be a function)
    radar_file_list = get_files(raddir + '/', julday)
    if len(radar_file_list) == 0:
        print_red('No radar file found for this date '+ julday.strftime("%d %b %Y"))
        return None

    print_yellow("%i radar files found." % (len(radar_file_list)))

    # Get the datetime for each radar files
    dtime_radar = [None]*len(radar_file_list)  # Allocate empty list
    for cnt, radfile in enumerate(radar_file_list):
        dtime_radar[cnt] = get_time_from_filename(radfile, date)

    dtime_radar = list(filter(None, dtime_radar))  # Removing None values

    if len(dtime_radar) == 0:
        print_red("No corresponding ground radar files for this date " + julday.strftime("%d %b %Y"))
        return None

    # Find the nearest scan time    )
    closest_dtime_rad = get_closest_date(dtime_radar, dtime_sat)

    if dtime_sat >= closest_dtime_rad:
        time_difference = dtime_sat - closest_dtime_rad
    else:
        time_difference = closest_dtime_rad - dtime_sat

    # Looking at the time difference between satellite and radar
    if time_difference.seconds > maxdt:
        print_red('Time difference is of %i s.' % (time_difference.seconds), bold=True)
        print_red('This time difference is bigger' +
              ' than the acceptable value of %i s.' % (maxdt))
        return None  # To the next satellite file

    # Radar file corresponding to the nearest scan time
    radfile = get_filename_from_date(radar_file_list, closest_dtime_rad)
    time = closest_dtime_rad  # Keeping the IDL program notation

    print_with_time('READING ' + radfile)
    radar = read_radar(radfile, l_atten, gr_reflectivity_offset)

    if radar is None:
        return None  # Bad dimensions message already printed

    ngate = radar['ngate']
    nbeam = radar['nbeam']
    ntilt = radar['ntilt']
    rg = radar['range']
    ag = radar['azang']
    eg = radar['elev_3d']
    elang = radar['elang']
    dr = radar['dr']
    reflectivity_ground_radar = radar['reflec']

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        reflectivity_ground_radar[reflectivity_ground_radar < 10] = np.NaN

    # Determine the Cartesian coordinates of the ground radar's pixels
    # rg, ag, eg = np.meshgrid(r_range, azang, elang, indexing='ij')
    zg = sqrt(rg**2 + (earth_gaussian_radius + z0)**2 + \
         2*rg*(earth_gaussian_radius + z0)*sin(pi/180*eg)) - earth_gaussian_radius
    sg = earth_gaussian_radius*np.arcsin(rg*cos(pi/180*eg)/(earth_gaussian_radius + zg))
    xg = sg*cos(pi/180*(90 - ag))
    yg = sg*sin(pi/180*(90 - ag))

    # Compute the volume of each radar bin
    volg = 1e-9*pi*dr*(pi/180*bwr/2*rg)**2

    # Convert S-band GR reflectivities to Ku-band
    refg_ku = reflectivity_conversion.convert_to_Ku(reflectivity_ground_radar, zg, zbb, l_cband)

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
    ref2 = np.zeros((nprof, ntilt)) + np.NaN  # PR reflec S-band, snow
    ref3 = np.zeros((nprof, ntilt)) + np.NaN  # PR reflec S-band, hail
    ref4 = np.zeros((nprof, ntilt)) + np.NaN  # GR reflectivity
    ref5 = np.zeros((nprof, ntilt)) + np.NaN  # GR reflectivity Ku-band
    iref1 = np.zeros((nprof, ntilt)) + np.NaN  # path-integrated PR reflec
    iref2 = np.zeros((nprof, ntilt)) + np.NaN  # path-integrated GR reflec
    stdv1 = np.zeros((nprof, ntilt)) + np.NaN  # STD of PR reflectivity
    stdv2 = np.zeros((nprof, ntilt)) + np.NaN  # STD of GR reflectivity

    '''Number of bins in sample'''
    ntot1 = np.zeros((nprof, ntilt), dtype=int)  # Total nb of PR bin in sample
    nrej1 = np.zeros((nprof, ntilt), dtype=int)  # Nb of rejected PR bin in sample
    ntot2 = np.zeros((nprof, ntilt), dtype=int)  # Total nb of GR bin in sample
    nrej2 = np.zeros((nprof, ntilt), dtype=int)  # Nb of rejected GR bin in sample
    vol1 = np.zeros((nprof, ntilt)) + np.NaN  # Total volume of PR bins in sample
    vol2 = np.zeros((nprof, ntilt)) + np.NaN  # Total volume of GR bins in sample

    # Compute the path-integrated reflectivities at every points
    nat_refp = 10**(reflectivity_satellite/10.0)  # In natural units
    nat_refg = 10**(reflectivity_ground_radar/10.0)
    irefp = np.fliplr(nancumsum(np.fliplr(nat_refp), 1))
    irefg = nancumsum(nat_refg)
    irefp = drt*(irefp - nat_refp/2)
    irefg = dr*(irefg - nat_refg/2)
    irefp = 10*np.log10(irefp)
    irefg = 10*np.log10(irefg)

    # Convert to linear units
    if not l_dbz:
        reflectivity_satellite = 10**(reflectivity_satellite/10.0)
        reflectivity_ground_radar = 10**(reflectivity_ground_radar/10.0)
        refp_ss = 10**(refp_ss/10.0)
        refp_sh = 10**(refp_sh/10.0)
        refg_ku = 10**(refg_ku/10.0)

    irefp = 10**(irefp/10.0)
    irefg = 10**(irefg/10.0)

    print_green("Starting comparison")

    # Loop over the TRMM/GPM profiles and Loop over the GR elevation scan
    for ii, jj in itertools.product(range(nprof), range(ntilt)):
        # Temporally kill warnings (because of nanmean)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            # Identify those PR bins which fall within the GR sweep
            ip = np.where((ep[ii, :] >= elang[jj] - bwr/2) &
                          (ep[ii, :] <= elang[jj] + bwr/2))

            # Store the number of bins
            ntot1[ii, jj] = len(ip)
            if len(ip) == 0:
                continue

            x[ii, jj] = np.mean(xp[ii, ip])
            y[ii, jj] = np.mean(yp[ii, ip])
            z[ii, jj] = np.mean(zp[ii, ip])

            # Compute the thickness of the layer
            nip = len(ip)
            dz[ii, jj] = nip*drt*cos(pi/180*alpha[ii, 0])

            # Compute the PR averaging volume
            vol1[ii, jj] = np.sum(volp[ii, ip])

            # Note the mean TRMM beam diameter
            ds[ii, jj] = pi/180*bwt*np.mean((zt - zp[ii, ip])/cos(pi/180*alpha[ii, ip]))

            # Note the radar range
            s = sqrt(x[ii, jj]**2 + y[ii, jj]**2)
            r[ii, jj] = (earth_gaussian_radius + z[ii, jj])*sin(s/earth_gaussian_radius)/cos(pi/180*elang[jj])

            # Check that sample is within radar range
            if r[ii, jj] + ds[ii, jj]/2 > rmax:
                continue

            # Extract the relevant PR data
            refp1 = reflectivity_satellite[ii, ip].flatten()
            refp2 = refp_ss[ii, ip].flatten()
            refp3 = refp_sh[ii, ip].flatten()
            irefp1 = irefp[ii, ip].flatten()

            # Average over those bins that exceed the reflectivity
            # threshold (linear average)

            ref1[ii, jj] = np.nanmean(refp1)
            ref3[ii, jj] = np.nanmean(refp2)
            ref4[ii, jj] = np.nanmean(refp3)
            iref1[ii, jj] = np.nanmean(irefp1)

            if not l_dbz:
                stdv1[ii, jj] = np.nanstd(10*np.log10(refp1))
            else:
                stdv1[ii, jj] = np.nanstd(refp1)

            # Note the number of rejected bins
            nrej1[ii, jj] = int(np.sum(np.isnan(refp1)))
            if ~np.isnan(stdv1[ii, jj]) and nip - nrej1[ii, jj] > 1:
                continue

            # Compute the horizontal distance to all the GR bins
            d = sqrt((xg[:, :, jj] - x[ii, jj])**2 + (yg[:, :, jj] - y[ii, jj])**2)

            # Find all GR bins within the SR beam
            igx, igy = np.where(d <= ds[ii, jj]/2)

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
            w = volg1*exp(-1*(d[igx, igy]/(ds[ii, jj]/2.))**2)
            w = w*refg1/refg2

            ref2[ii, jj] = np.nansum(w*refg1)/np.nansum(w)

            # if ref2[ii, jj] < minrefp:
            #     ref2[ii, jj] = np.NaN

            ref5[ii, jj] = np.nansum(w*refg2)/np.nansum(w)
            iref2[ii, jj] = np.nansum(w*irefg1)/np.nansum(w)

            if not l_dbz:
                stdv2[ii, jj] = np.nanstd(10*np.log10(refg1))
            else:
                stdv2[ii, jj] = np.nanstd(refg1)

            # Note the number of rejected bins
            nrej2[ii, jj] = int(np.sum(np.isnan(refg1)))

        # END WITH (RuntimeWarning ignore)
    # END FOR (satellite profiles, radar elevation)

    # Correct std
    stdv1[np.isnan(stdv1)] = 0
    stdv2[np.isnan(stdv2)] = 0

    # Convert back to dBZ
    iref1 = 10*np.log10(iref1)
    iref2 = 10*np.log10(iref2)

    if not l_dbz:
        reflectivity_satellite = 10*np.log10(reflectivity_satellite)
        reflectivity_ground_radar = 10*np.log10(reflectivity_ground_radar)
        refp_ss = 10*np.log10(refp_ss)
        refp_sh = 10*np.log10(refp_sh)
        refg_ku = 10*np.log10(refg_ku)
        ref1 = 10*np.log10(ref1)
        ref2 = 10*np.log10(ref2)
        ref3 = 10*np.log10(ref3)
        ref4 = 10*np.log10(ref4)
        ref5 = 10*np.log10(ref5)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ref2[ref2 < minrefg] = np.NaN

    # Extract comparison pairs
    ipairx, ipairy = np.where((~np.isnan(ref1)) & (~np.isnan(ref2)))
    if len(ipairx) < minpair:
        print_red('Insufficient comparison pairs for ' + julday.strftime("%d %b %Y"))
        return None    

    iprof = ipairx
    itilt = ipairy

    # Save structure
    match_vol = dict()

    match_vol['zbb'] = zbb
    match_vol['date'] = julday
    match_vol['bbwidth'] = bbwidth
    match_vol['dt'] = time_difference.seconds  # TODO CHECK!

    match_vol['x'] = x[ipairx, ipairy]
    match_vol['y'] = y[ipairx, ipairy]
    match_vol['z'] = z[ipairx, ipairy]
    match_vol['dz'] = dz[ipairx, ipairy]
    match_vol['ds'] = ds[ipairx, ipairy]
    match_vol['r'] = r[ipairx, ipairy]
    match_vol['el'] = elang[itilt]

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

    match_vol['sfc'] = sfc[iprof]
    match_vol['ptype'] = ptype[iprof]
    match_vol['iray'] = iray[iprof]
    match_vol['iscan'] = iscan[iprof]

    match_vol['stdv1'] = stdv1[ipairx, ipairy]
    match_vol['stdv2'] = stdv2[ipairx, ipairy]
    match_vol['vol1'] = vol1[ipairx, ipairy]
    match_vol['vol2'] = vol2[ipairx, ipairy]

    return match_vol
