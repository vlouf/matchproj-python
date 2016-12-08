def radar_gaussian_curve(lat0):
    '''RADAR_GAUSSIAN_CURVE'''
    '''Determine the Earth's Gaussian radius of curvature at the radar'''
    '''https://en.wikipedia.org/wiki/Earth_radius#Radii_of_curvature'''

    # Major and minor radii of the Ellipsoid
    a = 6378137.0 # in meters
    e2 = 0.0066943800
    b = a*sqrt(1-e2)

    tmp = (a*cos(pi/180*lat0))**2+(b*sin(pi/180*lat0))**2
    an = (a**2)/sqrt(tmp)
    am = (a*b)**2/tmp**(3/2.)
    ag = sqrt(an*am)
    ae = (4/3.)*ag
    return ae


def ground_radar_params(rname='CPOL'):
    """GROUND_RADAR_PARAMS"""
    """Input parameters:
            - rname: radar location name, type string
       Output:
            - Dictionnary containing the radar id (rid), latitude and longitude
              of the radar (lat0 and lon0, respectivelly), the radar height (z0),
              and the bwr
    """

    if rname == 'CPOL':
        rid='IDR59'
        lon0=131.0440
        lat0=-12.2490
        z0=50.
        bwr=1.0
    elif rname =='Berrimah':
        rid='IDR63'
        lon0=130.9250
        lat0=-12.4570
        z0=51.
        bwr=1.0
    elif rname == 'Brisbane':
        rid='IDR66'
        lon0=153.2400
        lat0=-27.7178
        z0=174.
        bwr=1.0
    elif rname == 'Marburg':
        rid='IDR50'
        lon0=152.5390
        lat0=-27.6080
        z0=372.0
        bwr=1.9
    elif rname == 'Grafton':
        rid='IDR28'
        lon0=152.951
        lat0=-29.622
        z0=40.0
        bwr=1.9
    elif rname == 'Gympie':
        rid='IDR08'
        lon0=152.577
        lat0=-25.957
        z0=375.0
        bwr=2.0
    else:
        raise ValueError("The available radars are: CPOL, Berrimah, Brisbane," +
                         " Marburg, Grafton, and Gympie")

    to_return = {'rid': rid, 'lon0': lon0, 'lat0': lat0, 'z0': z0, 'bwr': bwr}
    return to_return
