# Python Standard Library
import datetime

# Other libraries.
import netCDF4
import numpy as np
# from pyhdf.SD import SD, SDC

def read_date_from_TRMM(hdf_file1, radar_lat, radar_lon):
    """
    Extract datetime from TRMM HDF files.

    Parameters:
    ===========
    hdf_file1: str
        File name for TRMM satellite file.
    radar_lat: float
        Latitude of ground radar
    radar_lon: float
        Longitude of ground radar 
    
    Returns:
    ========
    trmm_date: datetime
        Datetime of satellite data at ground radar position.
    min_dist: float
        Minimal distance between satellite swath and ground radar, i.e. 
        is satellite swath are in ground radar domain?
    """
    with netCDF4.Dataset(hdf_file1, 'r') as ncid:
        year = ncid['Year'][:]
        month = ncid['Month'][:]
        day = ncid['DayOfMonth'][:]
        hour = ncid['Hour'][:]
        minute = ncid['Minute'][:]
        second = ncid['Second'][:]
        latitude = ncid['Latitude'][:]
        longitude = ncid['Longitude'][:]

    # Using distance, find min to radar
    dist = np.sqrt((latitude - radar_lat)**2 + (longitude - radar_lon)**2)
    dist_atrack = np.amin(dist, axis=1)  # Min distance along track axis
    radar_center = np.argmin(dist_atrack)
    min_dist = np.amin(dist_atrack)
    trmm_date = datetime.datetime(year[radar_center], month[radar_center], day[radar_center],
                                  hour[radar_center], minute[radar_center], second[radar_center])
    return trmm_date, min_dist


def read_trmm(hdf_file1, hdf_file2, sat_offset=None):
    '''
    Reads TRMM 2A23 and 2A25 data files.

    Parameters:
    ===========
    hdf_file1: str
        File name for TRMM 2A23 data.
    hdf_file2: str
        File name for TRMM 2A25 data.

    Returns:
    ========
    data_dict: dict
        Dictionnary containing all the needed data from the 2A23 and 2A25 files.
    '''
    with netCDF4.Dataset(hdf_file1, 'r') as ncid:
        year = ncid['Year'][:]
        month = ncid['Month'][:]
        day = ncid['DayOfMonth'][:]
        hour = ncid['Hour'][:]
        minute = ncid['Minute'][:]
        second = ncid['Second'][:]
        Latitude = ncid['Latitude'][:]
        Longitude = ncid['Longitude'][:]
        bbwidth = ncid['BBwidth'][:]
        HBB = ncid['HBB'][:]
        dataQuality = ncid['dataQuality'][:]
        rainFlag = ncid['rainFlag'][:]
        rainType = ncid['rainType'][:]
        status = ncid['status'][:]

    if dataQuality.max() != 0:
        return None

    with netCDF4.Dataset(hdf_file2, 'r') as ncid:
        # Latitude25 = ncid['Latitude'][:]
        # Longitude25 = ncid['Longitude'][:]
        correctZFactor = ncid['correctZFactor'][:]
        nscan, nray, nbin = correctZFactor.shape

    reflectivity = correctZFactor / 100.0

    #  Reverse direction along the beam
    if sat_offset is not None:
        print(f"{sat_offset}dB offset applied to TRMM reflectivity.")
        reflectivity = reflectivity[:, :, ::-1] + sat_offset
    else:
        reflectivity = reflectivity[:, :, ::-1]

    rainFlag[(rainFlag >= 10) & (rainFlag < 20)] = 1
    rainFlag[rainFlag >= 20] = 2

    ptype = rainType.copy()
    ptype[rainType >= 300] = 3
    ptype[(rainType >= 200) & (rainType < 300)] = 2
    ptype[(rainType >= 100) & (rainType < 200)] = 1
    ptype[rainType == -88] = 0
    ptype[rainType == -99] = 1

    # Extract the surface Type
    sfc = np.zeros(status.shape, dtype=int) - 1
    sfc[status == 168] = 0
    sfc[status % 10 == 0] = 1
    sfc[(status - 1) % 10 == 0] = 2
    sfc[(status - 2) % 10 == 0] = 3
    sfc[(status - 3) % 10 == 0] = 4
    sfc[(status - 4) % 10 == 0] = 5
    sfc[(status - 9) % 10 == 0] = 9

    # Extract 2A23 quality
    quality = np.zeros(status.shape, dtype=int) - 1
    quality[status == 168] = 0
    quality[status < 50] = 1
    quality[status >= 50] = 2

    data_dict = dict()
    data_dict = {'nscan': nscan, 'nray': nray, 'nbin': nbin, 'year': year,
                 'month': month, 'day': day, 'hour': hour, 'minute': minute,
                 'second': second, 'lon': Longitude, 'lat': Latitude,
                 'pflag': rainFlag, 'ptype': ptype, 'zbb': HBB, 'bbwidth': bbwidth,
                 'sfc': sfc, 'quality': quality, 'refl': reflectivity}

    return data_dict
