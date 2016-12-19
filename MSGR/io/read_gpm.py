import h5py
import numpy as np


def read_gpm(infile):
    """
    READ_GPM
    Read HDF5 GPM file with these parameters:
    - dataQuality
    - landSurfaceType
    - flagPrecip
    - flagBB
    - heightBB
    - widthBB
    - qualityBB
    - typePrecip
    - qualityTypePrecip
    - zFactorCorrected

    It will reverse direction along the beam for reflectivity so that the first
    value along the z-axis of the reflectivity array corresponds to the ground
    value and not the top of the atm.

    Returns a dictionnary containing the data.
    """

    with h5py.File(infile   , 'r') as file_id:
        obj_id = file_id['NS']

        lat = obj_id['Latitude'].value
        lon = obj_id['Longitude'].value

        # Read time data
        mem_id = obj_id['ScanTime']
        year = mem_id['Year'].value
        month = mem_id['Month'].value
        day = mem_id['DayOfMonth'].value
        hour = mem_id['Hour'].value
        minute = mem_id['Minute'].value
        second = mem_id['Second'].value

        # Read in the surface type and precipitation flag
        mem_id = obj_id['PRE']
        sfc = mem_id['landSurfaceType'].value
        pflag = mem_id['flagPrecip'].value

        # Read in the brightband and precipitation type data
        mem_id = obj_id['CSF']
        zbb= mem_id['heightBB'].value
        qbb = mem_id['qualityBB'].value
        qtype = mem_id['qualityTypePrecip'].value
        ptype = mem_id['typePrecip'].value
        bbwidth = mem_id['widthBB'].value

        # Read in the data quality
        mem_id = obj_id['scanStatus']
        quality = mem_id['dataQuality'].value

        # Read in the reflectivity data
        mem_id = obj_id['SLV']
        refl = mem_id['zFactorCorrected'].value

    # Determine the dimensions
    if refl.ndim != 3:
        print("Invalid number of dimensions")
        return None

    nscan, nray, nbin = refl.shape

    #  Reverse direction along the beam
    refl = refl[:, :, ::-1]

    # Change pflag=1 to pflag=2 to be consistent with 'Rain certain' in TRMM
    pflag[pflag == 1] = 2

    # Simplify the precipitation types
    ptype = ptype/10000000.0

    # Simplify the surface types
    imissx, imissy = np.where(sfc == -9999)
    nmiss = len(imissx)
    sfc = sfc/100+1
    if nmiss > 0:
        sfc[imissx, imissy] = 0

    # Set a quality indicator for the BB and precip type data
    quality = np.zeros((nscan, nray))
    quality[((qbb == 0) | (qbb == 1)) & (qtype == 1)] = 1
    quality[(qbb > 1) | (qtype > 1)] = 2

    # Store data in a dict
    to_return = dict()
    to_return = {'nscan':nscan, 'nray':nray, 'nbin':nbin, 'year':year, 'month':month,
                 'day':day,  'hour':hour, 'minute':minute, 'second':second, 'lon':lon,
                 'lat':lat, 'pflag':pflag, 'ptype':ptype, 'zbb':zbb, 'bbwidth':bbwidth,
                 'sfc':sfc, 'quality':quality, 'refl':refl}

    return to_return
