import numpy as np
import pyproj
import pandas as pd
import glob
import re
import datetime
from numpy import sqrt, cos, sin, tan, pi
from read_gpm import read_gpm


def convert_reflectivity_ku_to_sband(refp, zp, zbb, bbwidth, l_cband=1):
    """Convert to S-band using method of Cao et al. (2013)"""
    # Set coefficients for conversion from Ku-band to S-band
    #        Rain      90%      80%      70%      60%      50%      40%      30%      20%      10%     Snow
    as0=[ 4.78e-2, 4.12e-2, 8.12e-2, 1.59e-1, 2.87e-1, 4.93e-1, 8.16e-1, 1.31e+0, 2.01e+0, 2.82e+0, 1.74e-1]
    as1=[ 1.23e-2, 3.66e-3, 2.00e-3, 9.42e-4, 5.29e-4, 5.96e-4, 1.22e-3, 2.11e-3, 3.34e-3, 5.33e-3, 1.35e-2]
    as2=[-3.50e-4, 1.17e-3, 1.04e-3, 8.16e-4, 6.59e-4, 5.85e-4, 6.13e-4, 7.01e-4, 8.24e-4, 1.01e-3,-1.38e-3]
    as3=[-3.30e-5,-8.08e-5,-6.44e-5,-4.97e-5,-4.15e-5,-3.89e-5,-4.15e-5,-4.58e-5,-5.06e-5,-5.78e-5, 4.74e-5]
    as4=[ 4.27e-7, 9.25e-7, 7.41e-7, 6.13e-7, 5.80e-7, 6.16e-7, 7.12e-7, 8.22e-7, 9.39e-7, 1.10e-6, 0.00e+0]
    #        Rain      90%      80%      70%      60%      50%      40%      30%      20%      10%     Hail
    ah0=[ 4.78e-2, 1.80e-1, 1.95e-1, 1.88e-1, 2.36e-1, 2.70e-1, 2.98e-1, 2.85e-1, 1.75e-1, 4.30e-2, 8.80e-2]
    ah1=[ 1.23e-2,-3.73e-2,-3.83e-2,-3.29e-2,-3.46e-2,-2.94e-2,-2.10e-2,-9.96e-3,-8.05e-3,-8.27e-3, 5.39e-2]
    ah2=[-3.50e-4, 4.08e-3, 4.14e-3, 3.75e-3, 3.71e-3, 3.22e-3, 2.44e-3, 1.45e-3, 1.21e-3, 1.66e-3,-2.99e-4]
    ah3=[-3.30e-5,-1.59e-4,-1.54e-4,-1.39e-4,-1.30e-4,-1.12e-4,-8.56e-5,-5.33e-5,-4.66e-5,-7.19e-5, 1.90e-5]
    ah4=[ 4.27e-7, 1.59e-6, 1.51e-6, 1.37e-6, 1.29e-6, 1.15e-6, 9.40e-7, 6.71e-7, 6.33e-7, 9.52e-7, 0.00e+0]

    refp_ss = np.zeros(refp.shape) + np.NaN  # snow
    refp_sh = np.zeros(refp.shape) + np.NaN  # hail
    zmlt = zbb + bbwidth/2.    # APPROXIMATION!
    zmlb = zbb - bbwidth/2.    # APPROXIMATION!
    ratio = (zp-zmlb)/(zmlt-zmlb)

    iax, iay = np.where(ratio >= 1)
    if len(iax) > 0: # above melting layer
        dfrs = as0[10] + as1[10]*refp[iax, iay] + as2[10]*refp[iax, iay]**2 \
               + as3[10]*refp[iax, iay]**3 + as4[10]*refp[iax, iay]**4
        dfrh = ah0[10] + ah1[10]*refp[iax, iay] + ah2[10]*refp[iax, iay]**2 \
               + ah3[10]*refp[iax, iay]**3 + ah4[10]*refp[iax, iay]**4
        refp_ss[iax, iay] = refp[iax, iay] + dfrs
        refp_sh[iax, iay] = refp[iax, iay] + dfrh

    ibx, iby = np.where(ratio <= 0)
    if len(ibx) > 0: # below the melting layer
        dfrs = as0[0] + as1[0]*refp[ibx, iby] + as2[0]*refp[ibx, iby]**2 + \
               as3[0]*refp[ibx, iby]**3 + as4[0]*refp[ibx, iby]**4
        dfrh = ah0[0] + ah1[0]*refp[ibx, iby] + ah2[0]*refp[ibx, iby]**2 + \
               ah3[0]*refp[ibx, iby]**3 + ah4[0]*refp[ibx, iby]**4
        refp_ss[ibx, iby] = refp[ibx, iby] + dfrs
        refp_sh[ibx, iby] = refp[ibx, iby] + dfrh

    imx, imy = np.where((ratio > 0) & (ratio < 1))
    if len(imx) > 0: # within the melting layer
        ind = np.round(ratio[imx, imy]).astype(int)
        dfrs = as0[ind] + as1[ind]*refp[imx, imy] + as2[ind]*refp[imx, imy]**2 + \
               as3[ind]*refp[imx, imy]**3 + as4[ind]*refp[imx, imy]**4
        dfrh = ah0[ind] + ah1[ind]*refp[imx, imy] + ah2[ind]*refp[imx, imy]**2+ \
               ah3[ind]*refp[imx, imy]**3 + ah4[ind]*refp[imx, imy]**4
        refp_ss[imx, imy] = refp[imx, imy] + dfrs
        refp_sh[imx, imy] = refp[imx, imy] + dfrh

    # Jackson Tan's fix for C-band
    if l_cband == 1:
        deltas = (refp_ss - refp)*5.3/10.0
        refp_ss = refp + deltas
        deltah = (refp_sh - refp)*5.3/10.0
        refp_sh = refp + deltah

    return refp_ss, refp_sh


def radar_gaussian_curve(lat0):
    '''Determine the Earth's Gaussian radius of curvature at the radar'''
    '''https://en.wikipedia.org/wiki/Earth_radius#Radii_of_curvature'''

    # Major and minor radii of the Ellipsoid
    a = 6378137.0 # in meters
    e2 = 0.0066943800
    b = a*sqrt(1-e2)

    tmp = (a*cos(pi/180*lat0))**2+(b*sin(pi/180*lat0))**2
    an = (a**2)/sqrt(tmp)
    am = (a*b)**2/tmp**(3/2.)
    ag = sqrt(an*am)
    ae = (4/3.)*ag
    return ae


def get_orbit_number(the_file):
    '''Get the orbit number from filename (last serie of 6 consecutives numbers in filename)'''

    to_return = re.findall("[0-9]{6}", the_file)[-1] #Get orbit number
    return to_return


l_write=1    # Switch for writing out volume-matched data
l_cband=1    # Switch for C-band GR
l_netcdf=1   # Switch for NetCDF GR data
l_dbz=0      # Switch for averaging in dBZ
l_dp=1       # Switch for dual-pol data
l_gpm=1      # Switch for GPM PR data

if l_gpm == 0:
    satstr='trmm'
else:
    satstr='gpm'

# Set the data directories
raddir='/data/vlouf/cpol/20150219'
satdir='/data/vlouf/GPM_DATA'

#CPOL parameters
rname = 'CPOL'
rid = 'IDR59'
lon0 = 131.0440
lat0 = -12.2490
z0 = 50.
bwr = 1.0

# Start and end dates
date1='20150218'
date2='20150218'

if l_dbz == 0:
    outdir = raddir+'/'+rid+'/'+satstr+'_comp'
else:
    outdir = raddir+'/'+rid+'/'+satstr+'_comp_dbz'

outdir = outdir + '_new'

# Orbit parameters
if l_gpm == 0:
    zt = 402500.   # orbital height of TRMM (post boost)
    drt = 250.     # gate spacing of TRMM
else:
    zt = 407000.   # orbital height of GPM
    drt = 125.     # gate spacing of GPM
bwt=0.71

# Algorithm parameters and thresholds
rmin = 15000. # minimum GR range (m)
rmax = 150000. # maximum GR range (m)
minprof = 10   # minimum number of PR profiles with precip
maxdt = 300.   # maximum PR-GR time difference (s)
tscan = 90.    # approx. time to do first few tilts (s)
minrefg = 0.   # minimum GR reflectivity
minrefp = 18.  # minimum PR reflectivity
minpair = 10   # minimum number of paired samples


# Initialise error counters
ntot = 0
nerr = np.zeros((8,), dtype=int)

# Map Projection
# Options: projection transverse mercator, lon and lat of radar, and ellipsoid WGS84
smap = pyproj.Proj('+proj=tmerc +lon_0=131.0440 +lat_0=-12.2490 +ellps=WGS84')

# Note the lon,lat limits of the domain
xmin = -1*rmax
xmax = rmax
ymin = -1*rmax
ymax = rmax
lonmin, latmin = smap(xmin, ymin, inverse=True)
lonmax, latmax = smap(xmax, ymax, inverse=True)

# Gaussian radius of curvatur for the radar's position
ae = radar_gaussian_curve(lat0)

# Determine the Julian days to loop over
jul1 = datetime.datetime.strptime(date1, '%Y%m%d')
jul2 = datetime.datetime.strptime(date2, '%Y%m%d')
nday = jul2-jul1
#juls = jul1+LINDGEN(nday)

# Date loop
for the_date in pd.date_range(jul1, jul2):
    year = the_date.year
    month = the_date.month
    day = the_date.day
    date = "%i%02i%02i" % (year, month, day)

    # Note the Julian day corresponding to 00 UTC
    jul0 = datetime.datetime(year, month, day, 0, 0, 0)

    # Note the number of satellite overpasses on this day
    satfiles = glob.glob(satdir + '/*' + date + '*.HDF5')
    nswath = len(satfiles)

    if nswath == 0:
        print('No satellite swaths')
        nerr[0] += 1
        continue

    # File loop
    for the_file in satfiles:
        ntot += 1
        orbit = get_orbit_number(the_file)

        print("Orbit", orbit)

        sat = read_gpm(the_file)
        if sat is None:
            print('Bad satellite data')
            continue

        nscan = sat['nscan']
        nray = sat['nray']
        nbin = sat['nbin']
        yearp = sat['year']
        monthp = sat['month']
        dayp = sat['day']
        hourp = sat['hour']
        minutep = sat['minute']
        secondp = sat['second']
        lonp = sat['lon']
        latp = sat['lat']
        pflag = sat['pflag']
        ptype = sat['ptype']
        zbb = sat['zbb']
        bbwidth = sat['bbwidth']
        sfc = sat['sfc']
        quality = sat['quality']
        refp = sat['refl']

        # Convert to Cartesian coordinates
        res = smap(lonp, latp)
        xp = res[0]
        yp = res[1]

        # Identify profiles withing the domnain
        ioverx, iovery = np.where((xp >= xmin) & (xp <= xmax) &
                                  (yp >= ymin) & (yp <= ymax))

        if len(ioverx) == 0:
            nerr[1] += 1
            print("Insufficient satellite rays in domain.")
            continue

        # Note the first and last scan indices
        i1x, i1y = np.min(ioverx), np.min(iovery)
        i2x, i2y = np.max(ioverx), np.max(iovery)

        # Identify the coordinates of these points
        xf = xp[i1x:i2x]
        yf = yp[i1y:i2y]

        # Determine the date and time (in seconds since the start of the day) of the closest approach of TRMM to the GR
        xc = xp[:, 24] #Grid center
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
        timep= "%02i%02i%02i" % (hour, minute, second)
        julp = datetime.datetime(year, month, day, hour, minute, second)

        # Some julp bullshit

        # Compute the distance of every ray to the radar
        d = sqrt(xp**2 + yp**2)

        # Identify precipitating profiles within the radaar range limits
        iscan, iray = np.where((d >= rmin) & (d <= rmax) & (pflag == 2))
        nprof = len(iscan)
        if nprof < minprof:
            nerr[2] += 1
            print('Insufficient precipitating satellite rays in domain', nprof)

        # Note the scan and ray indices for these rays
        # iscan, iray = np.unravel_index(iprof, d.shape)

        #Extract data for these rays
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
            tmp[:, k] = (refp[:, :, k])[iscan, iray]

        refp = tmp

        # Note the scan angle for each ray
        alpha = np.abs(-17.04 + np.arange(nray)*0.71)
        alpha = alpha[iray]

        # Correct for parallax to get x, y, z coordinates
        # The next 30 lines took 2 days to write....

        # Remember Python's ways: unlike IDL, rebin cannot change the number of dimension.
        # the_range dimension is equal to nbin, and we nw wnat to copy it for nprof x nbin
        the_range_1d = np.arange(nbin)*drt
        the_range = np.zeros((nprof, nbin))
        for idx in range(0, nprof):
            the_range[idx, :] = the_range_1d[:]

        # alpha dim is nprof x 1 and now we want nprof x nbin
        # xc, yc, xp, yp dimensions are nprof x 1
        the_alpha = np.zeros((nprof, nbin))
        xc0 = np.zeros((nprof, nbin))
        yc0 = np.zeros((nprof, nbin))
        xp0 = np.zeros((nprof, nbin))
        yp0 = np.zeros((nprof, nbin))
        for idx in range(0, nbin):
            the_alpha[:, idx] = alpha[:]
            xc0[:, idx] = xc[:]
            yc0[:, idx] = yc[:]
            xp0[:, idx] = xp[:]
            yp0[:, idx] = yp[:]

        alpha = the_alpha
        zp = the_range*cos(pi/180.*the_alpha)
        ds = the_range*sin(pi/180.*the_alpha)
        ang = np.arctan2(yp0-yc0, xp0-xc0)
        dx = ds*cos(ang)
        dy = ds*sin(ang)
        xp = xp0+dx
        yp = yp0+dy

        if np.min(ds) < 0:
            continue

        # Compute the (approximate) volume of each PR bin
        rt = zt/cos(pi/180*alpha) - the_range
        volp = drt*(1.e-9)*pi*(rt*pi/180*bwt/2.)**2

        # Compute the ground-radar coordinates of the PR pixels
        sp = sqrt(xp**2 + yp**2)
        gamma = sp/ae
        ep = 180/pi*np.arctan((cos(gamma) - (ae + z0)/(ae + zp))/sin(gamma))
        rp = (ae + zp)*sin(gamma)/cos(pi/180*ep)
        ap = 90-180/pi*np.arctan2(yp, xp)

        # Determine the median brightband height
        ibb = np.where((zbb > 0) & (bbwidth > 0) & (quality == 1))[0] #1D arrays
        nbb = len(ibbx)
        if nbb >= minprof:
            zbb = np.median(zbb[ibb])
            bbwidth = np.median(bbwidth[ibb])
        else:
            nerr[3] += 1
            print('Insufficient bright band rays', nbb)
            continue

        # Set all values less than minrefp as missing
        ibadx, ibady = np.where(refp < minrefp) #WHERE(refp lt minrefp,nbad)
        if len(ibadx) > 0:
            refp[ibadx, ibady] = np.NaN

        # Convert to S-band using method of Cao et al. (2013)
        refp_ss, refp_sh = convert_reflectivity_ku_to_sband(refp, zp, zbb, bbwidth, l_cband)
