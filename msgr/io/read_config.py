# PSL
import time
import datetime
import configparser

# Other libraries.
import pyproj
import crayons

from numpy import sqrt, cos, sin, tan, pi


def _radar_gaussian_curve(lat0):
    '''
    RADAR_GAUSSIAN_CURVE
    Determine the Earth's Gaussian radius of curvature at the radar
    https://en.wikipedia.org/wiki/Earth_radius#Radii_of_curvature
    '''

    # Major and minor radii of the Ellipsoid
    a = 6378137.0  # Earth radius in meters
    e2 = 0.0066943800
    b = a * sqrt(1 - e2)

    tmp = (a * cos(pi / 180 * lat0))**2 + (b * sin(pi / 180 * lat0))**2   # Denominateur
    an = (a**2) / sqrt(tmp)  # Radius of curvature in the prime vertical (east–west direction)
    am = (a * b)**2 / tmp**1.5  # Radius of curvature in the north–south meridian
    ag = sqrt(an * am)  # Earth's Gaussian radius of curvature
    ae = (4 / 3.) * ag

    return ae


def _satellite_params(sname='GPM'):
    """
    SATELLITE_PARAMS
    return the gate spacing and the orbital height of the satellite
    """

    sname = sname.upper()
    # Orbit parameters
    if sname == 'GPM':
        zt = 407000.   # orbital height of GPM
        drt = 125.     # gate spacing of GPM
    elif sname == 'TRMM':
        zt = 402500.   # orbital height of TRMM (post boost)
        drt = 250.     # gate spacing of TRMM
    else:
        raise ValueError("The available satellites are GPM or TRMM.")
    bwt = 0.71

    return {'zt': zt, 'drt': drt, 'bwt': bwt}


def read_configuration_file(config_file):
    """
    Reading the configuration file.

    Parameters:
    ===========
        config_file: str
            Path to the configuration file.

    Returns:
    ========
        start_date: datetime
            Start date for treatment.
        end_date: datetime
            End date for treatment.
        ncpu: int
            Number of process for treatment.
        PARAMETERS_dict: dict
            Dictionnary containing all the configuration parameters.
    """
    #  Reading configuration file
    config = configparser.ConfigParser()
    config.read(config_file)

    general = config['general']
    ncpu = general.getint('ncpu')
    date1 = general.get('start_date')
    date2 = general.get('end_date')

    switch = config['switch']
    # Switch for writing out volume-matched data
    l_write = switch.getboolean('write')
    l_cband = switch.getboolean('cband')   # Switch for C-band GR
    l_dbz = switch.getboolean('dbz')       # Switch for averaging in dBZ
    l_gpm = switch.getboolean('gpm')       # Switch for GPM PR data
    # Switch for GPM PR data
    l_atten = switch.getboolean('correct_gr_attenuation')

    try:  # Optionnal parameters
        l_intermediary = switch.getboolean('intermediary')
        if l_intermediary is None:
            l_intermediary = False
    except KeyError:
        l_intermediary = False

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
    # minimum number of PR profiles with precip
    minprof = thresholds.getint('min_profiles')
    # maximum PR-GR time difference (s)
    maxdt = thresholds.getfloat('max_time_delta')
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
    pyproj_config = "+proj=tmerc +lon_0=%f +lat_0=%f +ellps=WGS84" % (
        lon0, lat0)
    smap = pyproj.Proj(pyproj_config)

    # Gaussian radius of curvatur for the radar's position
    earth_gaussian_radius = _radar_gaussian_curve(lat0)

    # Stocking parameters in dictionnaries
    if l_gpm:
        SAT_params = _satellite_params('gpm')
    else:
        SAT_params = _satellite_params('trmm')

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
    SWITCH_params['l_intermediary'] = l_intermediary

    PATH_params['raddir'] = raddir
    PATH_params['satdir'] = satdir
    PATH_params['outdir'] = outdir

    RADAR_params['xmin'] = -1 * rmax
    RADAR_params['xmax'] = rmax
    RADAR_params['ymin'] = -1 * rmax
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


def welcome_message(l_gpm, l_atten, l_dbz, l_write, outdir, satdir, raddir,
                    ncpu, start_date, end_date):
    '''
    WELCOME_MESSAGE
    Print a welcome message with a recap on the main global variables status
    '''

    msg = " " * 38 + "MSGR\n" + " " * 22 + "Matching Satellite and Ground Radar"

    print("#" * 80)
    print(crayons.blue("\n" + msg + "\n", bold=True))
    print("Volume matching program between GPM/TRMM spaceborne radar and ground radars.")
    if l_gpm:
        print("The spaceborne instrument used is GPM.")
    else:
        print("The spaceborne instrument used is TRMM.")
    print("The volume matching will be executed between " +
          start_date.strftime('%d %b %Y') + ' and ' + end_date.strftime('%d %b %Y'))
    if l_atten:
        print("Ground radar attenuation will be corrected.")
    else:
        print("Ground radar attenuation will NOT be corrected. I suppose it has already been done.")
    if l_dbz:
        print("The statistics will be done in dBZ.")
    else:
        print("The statistics will be done in natural units.")
    if l_write:
        print("The results will be saved in " + outdir)
    else:
        print("The results won't be saved.")
    print("This program will look for satellite data in " + satdir)
    print("This program will look for ground radar data in " + raddir)
    print("This program will run on %i cpu(s)." % (ncpu))
    print("#" * 80)
    print("\n\n")

    return None
