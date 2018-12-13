import datetime
import configparser

import numpy as np

from numpy import sqrt, cos, sin, tan, pi

from .io.read_gpm import read_gpm
from .io.read_trmm import read_trmm
from .utils import reflectivity_conversion


class Radar:
    def __init__(self, config_file, gr_offset=0):
        self.offset = gr_offset

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
        self.altitude = GR_param.getfloat('altitude')
        self.beamwidth = GR_param.getfloat('beamwidth')
        self.min_refl_thrld = config['thresholds'].getfloat('min_gr_reflec')
        self.l_cband = config['switch'].getboolean('cband')

        # Compute Earth gaussian radius.
        self.gaussian_radius = self._radar_gaussian_curvature()

        self.fields = dict()

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

    def set_fields(self, mydict):
        """
        Populate field dictionnary
        """
        for k, v in mydict.items():
            self.fields[k] = v

    def get_cartesian_coordinates(self):
        rg = self.fields['range']
        ag = self.fields['azang']
        eg = self.fields['elev_3d']
        # Determine the Cartesian coordinates of the ground radar's pixels
        zg = sqrt(rg**2 + (self.gaussian_radius + self.altitude)**2 +
                  2 * rg * (self.gaussian_radius + self.altitude) * sin(pi / 180 * eg)) - self.gaussian_radius
        sg = self.gaussian_radius * np.arcsin(rg * cos(pi / 180 * eg) / (self.gaussian_radius + zg))
        xg = sg * cos(pi / 180 * (90 - ag))
        yg = sg * sin(pi / 180 * (90 - ag))

        return xg, yg, zg

    def convert_refl_ku(self, zbb):
        xg, yg, zg = self.get_cartesian_coordinates()
        refg_ku = reflectivity_conversion.convert_to_Ku(self.fields['reflec'], zg, zbb, self.l_cband)
        return refg_ku


class Satellite:
    def __init__(self, config_file, sat_file_1, sat_file_2A25_trmm=None):
        config = self._read_configfile(config_file)
        thresholds = config['thresholds']
        try:
            sat_offset = config['satellite'].getfloat("sat_offset")
        except KeyError:
            sat_offset = None
            pass

        self.l_gpm = config['switch'].getboolean('gpm')
        if not self.l_gpm and sat_file_2A25_trmm is None:
            self.TRMM_NEW_VERSION = True
            # raise ValueError("Configuration file says that the satellite is TRMM but no TRMM 2A25 files given.")

        self.min_prof_nb = thresholds.getint('min_profiles')  # minprof
        self.max_time_delta = thresholds.getfloat('max_time_delta')  # maxdt
        self.min_refl_thrld = thresholds.getfloat('min_sat_reflec')  # minrefp
        self.min_pair_nb = thresholds.getint('min_pair')  # minpair

        # Orbit parameters
        if self.l_gpm:
            self.altitude = 407000.   # orbital height of GPM (zt)
            self.dr = 125.            # gate spacing of GPM (drt)
            satdata = read_gpm(sat_file_1, sat_offset)
        elif self.TRMM_NEW_VERSION:
            self.altitude = 402500.   # orbital height of TRMM (post boost)
            self.dr = 250.            # gate spacing of TRMM
            satdata = read_gpm(sat_file_1, sat_offset)
        else:
            self.altitude = 402500.   # orbital height of TRMM (post boost)
            self.dr = 250.            # gate spacing of TRMM
            satdata = read_trmm(sat_file_1, sat_file_2A25_trmm, sat_offset)

        self.beamwidth = 0.71

        if satdata is None:
            raise ValueError("Incorrect satellite data file.")

        self.__dict__.update(**satdata)

        # Determine the direction of the scans
        self.alpha = np.abs(-17.04 + np.arange(self.nray) * 0.71)

    def _read_configfile(self, config_file):
        config = configparser.ConfigParser()
        config.read(config_file)
        return config
