import datetime
import configparser

import numpy as np

from numpy import sqrt, cos, sin, tan, pi

from .io.read_gpm import read_gpm
from .io.read_trmm import read_trmm


class Radar:
    def __init__(self, config_file):
        # Read radar config file.
        config = self._read_configfile(config_file)
        GR_param = config['radar']

        self._check_keys(GR_param)  # raise error if wrong key.

        # Radar name
        self.name = GR_param.get('radar_name')
        self.id = GR_param.get('radar_id')

        # Radar lat/lon
        self.longitude = GR_param.getfloat('longitude')
        self.latitude = GR_param.getfloat('latitude')

        # Range
        self.rmin = GR_param.getfloat('rmin')
        self.rmax = GR_param.getfloat('rmax')

        # x/y min/max
        self.xmin = -1 * self.rmax
        self.xmax = self.rmax
        self.ymin = -1 * self.rmax
        self.ymax = self.rmax

        # others infos.
        self.altitude =  GR_param.getfloat('altitude')
        self.beamwidth = GR_param.getfloat('beamwidth')
        self.min_refl_thrld = config['thresholds'].getfloat('min_gr_reflec')
        
        try:
            self.offset = GR_param.getfloat('offset')
        except KeyError:
            self.offset = 0

        # Compute Earth gaussian radius.
        self.gaussian_radius = self._radar_gaussian_curvature()

    def _read_configfile(self, config_file):
        config = configparser.ConfigParser()
        config.read(config_file)
        return config

    def _check_keys(self, GR_param):
        keytab = ['radar_name', 'rmin', 'rmax', 'radar_id', 'longitude', 'latitude', 'altitude', 'beamwidth', 'offset']
        for mykey in keytab:
            try:
                GR_param[mykey]
            except KeyError:
                raise KeyError("Problem with configuration file, key: '{}' is missing and/or invalid.".format(mykey))
        return None

    def _radar_gaussian_curvature(self):
        '''
        Determine the Earth's Gaussian radius of curvature at the radar
        https://en.wikipedia.org/wiki/Earth_radius#Radii_of_curvature
        '''
        lat0 = self.latitude

        # Major and minor radii of the Ellipsoid
        a = 6378137.0  # Earth radius in meters
        e2 = 0.0066943800
        b = a * sqrt(1 - e2)

        tmp = (a * cos(pi / 180 * lat0))**2 + (b * sin(pi / 180 * lat0))**2   # Denominator
        an = (a**2) / sqrt(tmp)  # Radius of curvature in the prime vertical (east–west direction)
        am = (a * b)**2 / tmp**1.5  # Radius of curvature in the north–south meridian
        ag = sqrt(an * am)  # Earth's Gaussian radius of curvature
        ae = (4 / 3.) * ag

        return ae


class Satellite:
    def __init__(self, config_file, sat_file_1, sat_file_2A25_trmm=None):
        config = self._read_configfile(config_file)
        thresholds = config['thresholds']

        self.l_gpm = config['switch'].getboolean('gpm')
        if not self.l_gpm and sat_file_2A25_trmm is None:
            raise ValueError("Configuration file says that the satellite is TRMM but no TRMM 2A25 files given.")

        self.min_prof_nb = thresholds.getint('min_profiles')  # minprof
        self.max_time_delta = thresholds.getfloat('max_time_delta')  # maxdt
        self.min_refl_thrld = thresholds.getfloat('min_sat_reflec')  # minrefp
        self.min_pair_nb = thresholds.getint('min_pair')  # minpair

        sname = sname.upper()
        # Orbit parameters
        if self.l_gpm:
            self.altitude = 407000.   # orbital height of GPM (zt)
            self.dr = 125.            # gate spacing of GPM (drt)
            satdata = read_gpm(sat_file_1)
        else:
            self.altitude = 402500.   # orbital height of TRMM (post boost)
            self.dr = 250.            # gate spacing of TRMM
            satdata = read_trmm(sat_file_1, sat_file_2A25_trmm)
        self.beamwidth = 0.71

        if satdata is None:
            raise ValueError("Incorrect satellite data file.")

        self.__dict__.update(**satdata)

        self.alpha = np.abs(-17.04 + np.arange(self.nray) * 0.71)

    def _read_configfile(self, config_file):
        config = configparser.ConfigParser()
        config.read(config_file)
        return config

    def datetime_for_bin(self, position):
        dtime_sat = datetime.datetime(self.year[position], self.month[position],
                                      self.day[position], self.hour[position],
                                      self.minute[position], self.second[position])
        return dtime_sat
