import netCDF4


def _get_metadata():
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
    metadat['ref1'] = {'long_name': 'PR reflectivity', 'units': 'dBZ'}
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

    return metadat


def save_data(outfilename, data, date, offset1=None, offset2=None, nb_pass=0):
    """
    SAVE_DATA
    Dumps data in a python's pickle file
    Will try to populate the metadata based on the key name.

    Parameters:
    ===========
        outfilename: str
            Output file name.
        data: dict
            Dictionnary of data to save.
        do_hdf: bool
            Save as a HDF file.
    """
    metadat = _get_metadata()

    xdim = len(data['ref1'])
    tiltdim = len(data['el'])
    profdim = len(data['sfc'])

    with netCDF4.Dataset(outfilename, "w", format="NETCDF4") as rootgrp:
        # Create dimension
        rootgrp.createDimension("x", xdim)
        rootgrp.createDimension('tilt', tiltdim)
        rootgrp.createDimension('profile', profdim)
        rootgrp.createDimension('time', 1)
        time = rootgrp.createVariable("time", "f8", ('time'))
        time.units = "seconds since 1970-01-01T00:00:00Z"
        time[:] = netCDF4.date2num(date, "seconds since 1970-01-01T00:00:00Z")

        if nb_pass == 0:
            ncoff = rootgrp.createVariable("offset1", "f8", ("time"))
            ncoff[:] = offset1
            ncoff.setncattr_string("description", "Difference reflectivity Satellite - Ground Radar. PASS 1")
        else:
            ncoff = rootgrp.createVariable("offset1", "f8", ("time"))
            ncoff[:] = offset1
            ncoff.setncattr_string("description", "Difference reflectivity Satellite - Ground Radar. PASS 1")

            ncoff = rootgrp.createVariable("offset2", "f8", ("time"))
            ncoff[:] = offset2
            ncoff.setncattr_string("description", "Difference reflectivity Satellite - Ground Radar. PASS 2")

            ncoff = rootgrp.createVariable("offset_total", "f8", ("time"))
            ncoff[:] = offset2 + offset1
            ncoff.setncattr_string("description", "Difference reflectivity Satellite - Ground Radar. TOTAL")

        for k, v in data.items():
            if k in ['zbb', 'date', 'bbwidth', 'dt']:
                # rootgrp.setncattr(k, v)
                continue

            if k == "el":
                ncelev = rootgrp.createVariable('elevation', 'f8', ("tilt"), zlib=True)
                ncelev[:] = v
                ncelev.setncattr_string("long_name", metadat[k]['long_name'])
                if metadat[k]['units'] is not None:
                    ncelev.units = metadat[k]['units']
                continue

            if k in ['sfc', 'ptype', 'iray', 'iscan']:
                ncprof = rootgrp.createVariable(k, 'f8', ("profile"), zlib=True)
                ncprof[:] = v
                ncprof.setncattr_string("long_name", metadat[k]['long_name'])
                if metadat[k]['units'] is not None:
                    ncprof.units = metadat[k]['units']
                continue

            ncmoment = rootgrp.createVariable(k, 'f8', ("x",), zlib=True)
            ncmoment[:] = v
            try:
                if metadat[k]['units'] is not None:
                    ncmoment.units = metadat[k]['units']
                ncmoment.setncattr_string("long_name", metadat[k]['long_name'])
            except KeyError:
                pass

    return None
