import os
import datetime

import netCDF4
import numpy as np


def rem_outliers(x):
    """
    Remove outliers of a distribution.
    """
    ll = np.percentile(x, 5)
    ul = np.percentile(x, 95)
    return x[(x >= ll) & (x <= ul)]


def compute_offset(infile):
    with netCDF4.Dataset(infile) as ncid:
        ref1 = ncid['ref1'][:]
        ref2 = ncid['ref2'][:]
        ref5 = ncid['ref5'][:]
        std_sat = ncid['stdv1'][:]
        std_rad = ncid['stdv2'][:]
        z = ncid['z'][:]

    pos = (ref1 >= 40) | (std_sat > 4) | (std_rad > 4) | (ref2 >= 36) | (ref1 == 0) | (ref2 < 21)
    ref1[pos] = np.NaN

    dref_ku = (ref5 - ref1)
    dref_ku = dref_ku[~np.isnan(dref_ku)]
    # dref_ku = rem_outliers(dref_ku)

    if len(dref_ku) <= 1:
        return np.Nan
    else:
        offset = - np.median(dref_ku)  # !!! Note the MINUS sign !!!
        return offset
