import pickle

def save_data(out_file, data):
    '''
    SAVE_DATA
    Dumps data in a python's pickle file
    Will try to populate the metadata based on the key name.
    '''

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

    # Opening file and dumping data in it.
    with open(out_file + ".pkl", 'wb') as fid:
        pickle.dump(to_save, fid)

    return None
