import numpy as np
import pyproj
import pandas as pd
import glob
import re
import os
import datetime
import pyart
from numpy import sqrt, cos, sin, tan, pi

# Custom modules
import reflectivity_conversion
from read_gpm import read_gpm
from ground_radar import *
from satellite import *


def get_files(inpath):
    '''GET_FILES'''
    '''Returns a list of with the supported extension (netcdf) in the given
    path. Will recursively search in subdirectories too.'''

    supported_extension = ['.nc', '.NC', '.cdf']
    flist = []

    for dirpath, dirnames, filenames in os.walk(inpath):
        for filenames_slice in filenames:
            file_extension = os.path.splitext(str(filenames_slice))[1]
            # Get extension

            if np.any(np.in1d(supported_extension, file_extension)):
                # Check if file extension is in the list of supported ones
                the_path = os.path.join(dirpath, filenames_slice)
            else:  # If not test next file.
                continue

            # File does have the supported extension, we keep it for returning
            # list
            flist.append(the_path)

    to_return = flist

    return sorted(to_return)  # Type: List[str, ...]


def get_time_from_filename(filename, date):
    '''GET_TIME_FROM_FILENAME'''
    '''Capture the time string inside the filename and returns it'''

    # Looking for date followed by underscore and 6 consecutives number (i.e.
    # the time)
    date_time_str = re.findall(date + '_[0-9]{6}', filename)[0]
    # Turn it into a datetime object
    to_return = datetime.datetime.strptime(date_time_str, '%Y%m%d_%H%M%S')

    return to_return  # Type: str


def get_closest_date(list_date, base_time):
    '''GET_CLOSEST_DATE'''
    # from:  http://stackoverflow.com/a/17249470/846892

    b_d = base_time
    def func(x):
        dd =  x
        delta =  dd - b_d if dd > b_d else datetime.timedelta.max
        return delta

    # There is some black magic going on here...
    return min(list_date, key=func)  # Type: datetime


def get_filename_from_date(file_list, the_date):
    '''GET_FILENAME_FROM_DATE'''
    '''Looks for a file in a list of file with the exact corresponding date and
       returns it'''

    rt_str = the_date.strftime("%Y%m%d_%H%M%S")
    for the_file in file_list:
        try:
            re.findall(rt_str, the_file)[0]
            to_return = the_file
            break
        except IndexError:
            continue

    return to_return  # Type: str


""" SECTION of user-defined parameters """
l_write = 1    # Switch for writing out volume-matched data
l_cband = 1    # Switch for C-band GR
l_netcdf = 1   # Switch for NetCDF GR data
l_dbz = 0      # Switch for averaging in dBZ
l_dp = 1       # Switch for dual-pol data
l_gpm = 1      # Switch for GPM PR data

# Start and end dates
date1 = '20150211'
date2 = '20150211'

# Set the data directories
raddir = '/g/ns/cw/arm/data-1/vlouf/cpol_season_1415'
satdir = '/data/vlouf/GPM_DATA'

# Algorithm parameters and thresholds
rmin = 15000.  # minimum GR range (m)
rmax = 150000  # maximum GR range (m)
minprof = 10   # minimum number of PR profiles with precip
maxdt = 300.   # maximum PR-GR time difference (s)
tscan = 90.    # approx. time to do first few tilts (s)
minrefg = 0.   # minimum GR reflectivity
minrefp = 18.  # minimum PR reflectivity
minpair = 10   # minimum number of paired samples
""" End of the section for user-defined parameters """

if l_gpm == 0:
    satstr = 'trmm'
    raise ValueError("TRMM not yet implemented")
else:
    satstr = 'gpm'

# Ground radar parameters
GR_param = ground_radar_params('CPOL')
rid = GR_param['rid']
lon0 = GR_param['lon0']
lat0 = GR_param['lat0']
z0 = GR_param['z0']
bwr = GR_param['bwr']

sat_params = satellite_params(satstr)
zt = sat_params['zt']
drt = sat_params['drt']
bwt = sat_params['bwt']

# Output directory
if l_dbz == 0:
    outdir = raddir+'/'+rid+'/'+satstr+'_comp'
else:
    outdir = raddir+'/'+rid+'/'+satstr+'_comp_dbz'

outdir = outdir + '_new'

# Initialise error counters
ntot = 0
nerr = np.zeros((8,), dtype=int)

# Map Projection
# Options: projection transverse mercator, lon and lat of radar, and ellipsoid
# WGS84
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

        print("Orbit " + orbit + " -- " + jul0.strftime("%d %B %Y"))

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
            nerr[2] += 1
            print('Insufficient precipitating satellite rays in domain', nprof)

        # Note the scan and ray indices for these rays
        # iscan, iray = np.unravel_index(iprof, d.shape)

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
            tmp[:, k] = (refp[:, :, k])[iscan, iray]

        refp = tmp

        # Note the scan angle for each ray
        alpha = np.abs(-17.04 + np.arange(nray)*0.71)
        alpha = alpha[iray]

        # Correct for parallax to get x, y, z coordinates
        # The next 30 lines took 2 days to write....

        # Remember Python's ways: unlike IDL, rebin cannot change the number
        # of dimension. the_range dimension is equal to nbin, and we nw wnat
        # to copy it for nprof x nbin
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

        if len(ds) == 0:
            continue
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
        # 1D arrays
        ibb = np.where((zbb > 0) & (bbwidth > 0) & (quality == 1))[0]
        nbb = len(ibb)
        if nbb >= minprof:
            zbb = np.median(zbb[ibb])
            bbwidth = np.median(bbwidth[ibb])
        else:
            nerr[3] += 1
            print('Insufficient bright band rays', nbb)
            continue

        # Set all values less than minrefp as missing
        ibadx, ibady = np.where(refp < minrefp)  # WHERE(refp lt minrefp,nbad)
        if len(ibadx) > 0:
            refp[ibadx, ibady] = np.NaN

        # Convert to S-band using method of Cao et al. (2013)
        if l_cband:
            refp_ss, refp_sh = reflectivity_conversion.convert_to_Cband(refp, zp, zbb, bbwidth)
        else:
            refp_ss, refp_sh = reflectivity_conversion.convert_to_Sband(refp, zp, zbb, bbwidth)

        # Get the ground radar file lists (next 20 lines can be a function)
        radar_file_list = get_files(raddir + '/' + date + '/')

        # Get the datetime for each radar files
        dtime_radar = [None]*len(radar_file_list)  # Allocate empty list
        for cnt, radfile in enumerate(radar_file_list):
            dtime_radar[cnt] = get_time_from_filename(radfile, date)

        # Find the nearest scan time
        closest_dtime_rad = get_closest_date(dtime_radar, dtime_sat)

        if dtime_sat >= closest_dtime_rad:
            time_difference = dtime_sat - closest_dtime_rad
        else:
            time_difference = closest_dtime_rad - dtime_sat

        # Looking at the time difference between satellite and radar
        if time_difference.seconds > maxdt:
            print('Time difference is of %i.' % (time_difference.seconds))
            print('This time difference is bigger than the acceptable value of ', maxdt)            
            nerr[5] += 1
            continue  # To the next satellite file

        # Radar file corresponding to the nearest scan time
        radfile = get_filename_from_date(radar_file_list, closest_dtime_rad)
        time = closest_dtime_rad  # Keeping the IDL program notation

        radar = pyart.io.read(radfile)
