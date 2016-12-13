"""
This module regroups a bunch of usefull functions:
    - nancumsum (numpy's cumsum with support of NaN values)
    - get_files (walks into directories and returns all the files with supported extension)
    - get_time_from_filename (uses regular expression to catch the date in a filename)
    - get_closest_date (get the closest date in a list)
    - get_filename_from_date (takes a list of files and a date and returns the file closest to this date)
"""

import re
import os
import datetime
import numpy as np
import copy


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

    supported_extension = ['.nc', '.NC', '.cdf', '.hdf5', '.h5', '.HDF5', '.H5']
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
    try:
        date_time_str = re.findall(date + '_[0-9]{6}', filename)[0]
    except IndexError:
        return None
    # Turn it into a datetime object
    to_return = datetime.datetime.strptime(date_time_str, '%Y%m%d_%H%M%S')

    return to_return  # Type: str


def get_closest_date(list_date, base_time):
    '''GET_CLOSEST_DATE'''
    """from:  http://stackoverflow.com/a/17249470/846892"""

    b_d = base_time

    def func(x):
        dd = x
        delta = dd - b_d if dd > b_d else datetime.timedelta.max
        return delta

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
