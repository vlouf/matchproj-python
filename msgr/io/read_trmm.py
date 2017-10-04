# Python Standard Library
import copy

# Other libraries.
import numpy as np
from pyhdf.SD import SD, SDC


def read_trmm(hdf_file1, hdf_file2):
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
    dataQuality = hdf.select('dataQuality').get()
    rainFlag = hdf.select('rainFlag').get()
    rainType = hdf.select('rainType').get()
    status = hdf.select('status').get()
    hdf.end()

    if dataQuality.max() != 0:
        return None

    hdf_25 = SD(hdf_file2, SDC.READ)
    Latitude25 = hdf_25.select('Latitude').get()
    Longitude25 = hdf_25.select('Longitude').get()
    correctZFactor = hdf_25.select('correctZFactor').get()
    # dataQuality = hdf_25.select('dataQuality').get()

    nscan = hdf_25.select('correctZFactor').dimensions()['nscan']
    nray = hdf_25.select('correctZFactor').dimensions()['nray']
    nbin = hdf_25.select('correctZFactor').dimensions()['ncell1']
    hdf_25.end()

    reflectivity = correctZFactor / 100.0

    #  Reverse direction along the beam
    reflectivity = reflectivity[:, :, ::-1]

    rainFlag[(rainFlag >= 10) & (rainFlag < 20)] = 1
    rainFlag[rainFlag >= 20] = 2

    ptype = copy.deepcopy(rainType)
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

# Pfiou!
