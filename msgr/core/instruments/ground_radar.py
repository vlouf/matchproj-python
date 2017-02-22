from numpy import sqrt, cos, sin, tan, pi
from numba import jit

@jit
def radar_gaussian_curve(lat0):
    '''
    RADAR_GAUSSIAN_CURVE
    Determine the Earth's Gaussian radius of curvature at the radar
    https://en.wikipedia.org/wiki/Earth_radius#Radii_of_curvature
    '''

    # Major and minor radii of the Ellipsoid
    a = 6378137.0 # in meters
    e2 = 0.0066943800
    b = a*sqrt(1-e2)

    tmp = (a*cos(pi/180*lat0))**2+(b*sin(pi/180*lat0))**2   # Le denominateur
    an = (a**2)/sqrt(tmp)  # Radius of curvature in the prime vertical (east–west direction)
    am = (a*b)**2/tmp**(3/2.)  # Radius of curvature in the north–south meridian
    ag = sqrt(an*am)  # Earth's Gaussian radius of curvature
    ae = (4/3.)*ag
    return ae
