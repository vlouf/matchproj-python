import numpy as np
from numba import jit  # It improves speed, as somebody told me, not yet convinced
from numpy import sqrt, cos, sin, tan, pi

@jit
def correct_parallax(xc, yc, xp, yp, alpha, the_range):
    """
    CORRECT_PARALLAX
    Correct for parallax to get x, y, z coordinates
    This 'simple' function took 2 days to write....

    Remember Python's ways: unlike IDL, rebin cannot change the number
    of dimension. the_range dimension is equal to nbin, and we now wnat
    to copy it for nprof x nbin

    alpha dim is nprof x 1 and now we want nprof x nbin
    xc, yc, xp, yp dimensions are nprof x 1
    """

    nprof, nbin = the_range.shape
    the_alpha = np.zeros((nprof, nbin))
    xc0 = np.zeros((nprof, nbin))
    yc0 = np.zeros((nprof, nbin))
    xp0 = np.zeros((nprof, nbin))
    yp0 = np.zeros((nprof, nbin))
    for idx in range(0, nbin):
        the_alpha[:, idx] = alpha[:]
        xc0[:, idx] = xc[:]
        yc0[:, idx] = yc[:]
        xp0[:, idx] = xp[:]
        yp0[:, idx] = yp[:]

    alpha = the_alpha
    zp = the_range*cos(pi/180.*the_alpha)
    ds = the_range*sin(pi/180.*the_alpha)
    ang = np.arctan2(yp0-yc0, xp0-xc0)
    dx = ds*cos(ang)
    dy = ds*sin(ang)
    xp = xp0+dx
    yp = yp0+dy

    return xp, yp, zp, ds, alpha
