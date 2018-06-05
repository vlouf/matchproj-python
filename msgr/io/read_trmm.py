# Python Standard Library
import datetime

# Other libraries.
import numpy as np
import netCDF4
from pyhdf.SD import SD, SDC

def read_date_from_TRMM(hdf_file1, radar_lat, radar_lon):
    """
    Extract datetime from TRMM HDF files.
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
    is_dataQuality = True
    
    hdf = SD(hdf_file1, SDC.READ)
    year = hdf.select('Year').get()
    month = hdf.select('Month').get()
    day = hdf.select('DayOfMonth').get()
    hour = hdf.select('Hour').get()
    minute = hdf.select('Minute').get()
    second = hdf.select('Second').get()
    Latitude = hdf.select('Latitude').get()
    Longitude = hdf.select('Longitude').get()
    bbwidth = hdf.select('BBwidth').get()
    HBB = hdf.select('HBB').get()
    rainFlag = hdf.select('rainFlag').get()
    rainType = hdf.select('rainType').get()
    status = hdf.select('status').get()
    try:
        dataQuality = hdf.select('dataQuality').get()
    except:
        is_dataQuality = False
    hdf.end()    
    
    hdf_25 = SD(hdf_file2, SDC.READ)
    Latitude25 = hdf_25.select('Latitude').get()
    Longitude25 = hdf_25.select('Longitude').get()
    correctZFactor = hdf_25.select('correctZFactor').get()
    if not is_dataQuality:
        dataQuality = hdf_25.select('dataQuality').get()
    nscan = hdf_25.select('correctZFactor').dimensions()['nscan']
    nray = hdf_25.select('correctZFactor').dimensions()['nray']
    nbin = hdf_25.select('correctZFactor').dimensions()['ncell1']
    hdf_25.end()

    if dataQuality.max() != 0:
        raise ValueError('TRMM data quality are bad.')

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
