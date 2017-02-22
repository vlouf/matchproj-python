"""
This module regroups a bunch of usefull functions:
    - find_file_with_string (return the element of a list containing a specific string)
    - nancumsum (numpy's cumsum with support of NaN values)
    - get_files (walks into directories and returns all the files with supported extension)
    - get_time_from_filename (uses regular expression to catch the date in a filename)
    - get_closest_date (get the closest date in a list)
    - get_filename_from_date (takes a list of files and a date and returns the file closest to this date)
"""

import re
import os
import copy
import crayons  # Color terminal
import datetime
import numpy as np
from dateutil import parser


def find_file_with_string(flist, orb):
    """
    FIND_FILE_WITH_STRING
    Return the element of a list flist containing the value of orb

    Parameters
    ==========
        flist: list[str]
            List of file list
        orb: str
            String we are looking for in the list

    Returns
    =======
        Element of a list flist containing the value of orb
    """
    return [fd for fd in flist if orb in fd][0]


def nancumsum(a, ax=0):
    '''
    NANCUMSUM
    Cumsum in numpy does not ignore the NaN values, this one does
    Note that nancumsum will be implemented in numpy v1.12
    '''

    tmp = copy.deepcopy(a)
    tmp[np.isnan(tmp)] = 0
    rslt = np.cumsum(tmp, axis=ax)
    rslt[np.isnan(a)] = np.NaN

    return rslt

def get_files(inpath, date=None):
    '''
    GET_FILES
    Returns a list of with the supported extension (netcdf) in the given
    path. Will recursively search in subdirectories too. If provided a date
    (string or datetime object) it will only returns the files whose
    filename matches.
    '''

    supported_extension = ['.nc', '.NC', '.cdf', '.hdf5', '.h5', '.HDF5',
                           '.H5', '.lassen', '.PPI', '.UF']
    flist = []

    # Check date type
    if type(date) == datetime.datetime:
        date = date.strftime("%Y%m%d")

    for dirpath, dirnames, filenames in os.walk(inpath):
        for filenames_slice in filenames:

            # If no date provided, nothing new under the sun
            if date is None:
                pass  # pretends there was no if statement
            elif date in filenames_slice:
                pass  # pretends there was no if statement
            else:
                continue

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
    '''
    GET_TIME_FROM_FILENAME
    Capture the time string inside the filename and returns it.
    '''

    # Looking for date followed by underscore and 6 consecutives number (i.e.
    # the time)
    try:
        # There is maybe an optionnal character (like _) between date and time
        date_time_str = re.findall(date + ".?[0-9]{6}", filename)[0]
        to_return = parser.parse(date_time_str, fuzzy=True)
    except IndexError:
        to_return = None

    return to_return  # Type: str


def get_closest_date(list_date, base_time):
    '''
    GET_CLOSEST_DATE
    from:  http://stackoverflow.com/a/17249470/846892
    '''

    b_d = base_time

    def func(x):
        dd = x
        delta = dd - b_d if dd > b_d else datetime.timedelta.max
        return delta

    return min(list_date, key=func)  # Type: datetime


def get_filename_from_date(file_list, the_date):
    '''
    GET_FILENAME_FROM_DATE
    Looks for a file in a list of file with the exact corresponding date and
    returns it.
    '''

    # There is maybe an optionnal character(underscore) between date and time
    rt_str = the_date.strftime("%Y%m%d.?%H%M%S")
    for the_file in file_list:
        try:
            re.findall(rt_str, the_file)[0]  # If does not exist it raises an error
            to_return = the_file
            break  # We found what we are looking for, exiting the loop
        except IndexError:
            continue

    return to_return  # Type: str


def chunks(l, n):
    """
    Yield successive n-sized chunks from l.
    Use it to cut a big list into smaller chunks. => memory efficient
    """
    for i in range(0, len(l), n):
        yield l[i:i + n] # type: Generator[list of strings]


def print_with_time(txt):
    '''
    PRINT_WITH_TIME
    '''
    pfix = "[" + str(datetime.datetime.now().isoformat()) + "]\t"
    print(crayons.blue(pfix, bold=True) + txt)
    return None


# To print in color in the terminal. Pretty much self-explanatory.
def print_red(txt, bold=False):
    print_with_time(crayons.red(txt, bold))
    return None


def print_green(txt, bold=False):
    print_with_time(crayons.green(txt, bold))
    return None


def print_yellow(txt, bold=False):
    print_with_time(crayons.yellow(txt, bold))
    return None


def print_blue(txt, bold=False):
    print(crayons.blue(txt, bold))
    return None


def print_magenta(txt, bold=False):
    print(crayons.magenta(txt, bold))
    return None


def welcome_message(l_gpm, l_atten, l_dbz, l_write, outdir, satdir, raddir,
                    ncpu, start_date, end_date):
    '''
    WELCOME_MESSAGE
    Print a welcome message with a recap on the main global variables status
    '''

    msg = " "*38 + "MSGR\n" + " "*22 + "Matching Satellite and Ground Radar"

    print("#"*80)
    print_magenta("\n" + msg + "\n", bold=True)
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
    print("#"*80)
    print("\n\n")

    return None
