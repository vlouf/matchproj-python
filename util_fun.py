import re
import os
import datetime
import numpy as np
import pickle
import copy


def save_data(out_file, data):
    '''SAVE_DATA'''
    '''Save data in file'''

    metadat = dict()
    metadat['iscan'] = {'long_name': 'PR scan index', 'units': None}
    metadat['iray'] = {'long_name': 'PR ray index', 'units': None}
    metadat['itilt'] = {'long_name': 'GR tilt index', 'units': None}
    metadat['x'] = {'long_name': 'W-E location w.r.t. radar', 'units': 'm'}
    metadat['y'] = {'long_name': 'S-N location w.r.t. radar', 'units': 'm'}
    metadat['z'] = {'long_name': 'Height above MSL', 'units': 'm'}
    metadat['r'] = {'long_name': 'Range from GR', 'units': 'm'}
    metadat['el'] = {'long_name': 'Elevation angle', 'units': 'deg'}
    metadat['dz'] = {'long_name': 'Depth of averaging volume', 'units': 'm'}
    metadat['ds'] = {'long_name': 'Diameter of averaging volume', 'units': 'm'}
    metadat['dt'] = {'long_name': 'Time difference between GR and PR samples', 'units': 's'}
    metadat['ntot1'] = {'long_name': 'Number of PR points in averaging volume', 'units': None}
    metadat['ntot2'] = {'long_name': 'Number of GR points in averaging volume', 'units': None}
    metadat['nrej1'] = {'long_name': 'Number of rejected PR points in averaging volume', 'units': None}
    metadat['nrej2'] = {'long_name': 'Number of rejected GR points in averaging volume', 'units': None}
    metadat['sfc'] = {'long_name': 'Surface type (1=Ocean, 2=Land, 3=Coast, 4=Lake, 5=Unknown)', 'units': None}
    metadat['ptype'] = {'long_name': 'Precipitation type (1=Strat, 2=Conv, 3=Other)', 'units': None}
    metadat['ref1'] = {'long_name': 'PR reflectivity', 'units': 'dBZ' }
    metadat['ref2'] = {'long_name': 'GR reflectivity', 'units': 'dBZ'}
    metadat['ref3'] = {'long_name': 'PR reflectivity (S-band, Snow)', 'units': 'dBZ'}
    metadat['ref4'] = {'long_name': 'PR reflectivity (S-band, Hail)', 'units': 'dBZ'}
    metadat['ref5'] = {'long_name': 'GR reflectivity (Ku-band)', 'units': 'dBZ'}
    metadat['stdv1'] = {'long_name': 'Standard deviation of PR reflectivity', 'units': 'dB'}
    metadat['stdv2'] = {'long_name': 'Standard deviation of GR reflectivity', 'units': 'dB'}
    metadat['iref1'] = {'long_name': 'Path-integrated PR reflectivity', 'units': 'dB'}
    metadat['iref2'] = {'long_name': 'Path-integrated GR reflectivity', 'units': 'dB'}
    metadat['zbb'] = {'long_name': 'Average bright band height', 'units': 'm'}
    metadat['bbwidth'] = {'long_name': 'Average bright band width', 'units': 'm'}
    metadat['nbb'] = {'long_name': 'Number of profiles with a bright band', 'units': None}
    metadat['vol1'] = {'long_name': 'PR averaging volume', 'units': 'km^3'}
    metadat['vol2'] = {'long_name': 'GR averaging volume', 'units': 'km^3'}

    to_save = dict()
    for k in data.keys():
        try:
            to_save[k] = {'data': data[k], 'long_name': metadat[k]['long_name'], 'units':metadat[k]['units']}
        except KeyError:
            to_save[k] = {'data': data[k], 'long_name': None, 'units':None}

    with open(out_file + ".pkl", 'wb') as fid:
        pickle.dump(to_save, fid)

    return None


def nancumsum(a, ax=0):
    '''NANCUMSUM'''
    '''Cumsum in numpy does not ignore the NaN values, this one does'''
    '''Note that nancumsum will be implemented in numpy v1.12'''

    tmp = copy.deepcopy(a)
    tmp[np.isnan(tmp)] = 0
    rslt = np.cumsum(tmp, axis=ax)
    rslt[np.isnan(a)] = np.NaN

    return rslt


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
        dd = x
        delta = dd - b_d if dd > b_d else datetime.timedelta.max
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
            re.findall(rt_str, the_file)[0]  # If does not exist it raises an error
            to_return = the_file
            break
        except IndexError:
            continue

    return to_return  # Type: str
