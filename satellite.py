def satellite_params(sname='GPM'):
    """SATELLITE_PARAMS"""
    
    sname = sname.upper()
    # Orbit parameters
    if sname == 'GPM':
        zt = 402500.   # orbital height of TRMM (post boost)
        drt = 250.     # gate spacing of TRMM
    elif sname == 'TRMM':
        zt = 407000.   # orbital height of GPM
        drt = 125.     # gate spacing of GPM
    else:
        raise ValueError("The available satellites are GPM or TRMM.")
    bwt=0.71

    return {'zt':zt, 'drt':drt, 'bwt':bwt}


def get_orbit_number(the_file):
    """GET_ORBIT_NUMBER"""
    '''Get the orbit number from filename (last serie of 6 consecutives numbers
       in filename)'''

    to_return = re.findall("[0-9]{6}", the_file)[-1] #Get orbit number
    return to_return
