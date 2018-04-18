# Python Standard Library
import copy
import warnings
import traceback

# Other libraries.
import pyart
import netCDF4
import numpy as np


def get_reflectivity_field_name(radar):
    '''
    GET_REFLECTIVITY_FIELD_NAME
    Because of different conventions for naming fields, it will try a variety
    of different reflectvity name and return the first one that works
    '''

    potential_name = ['DBZ_F', 'DBZ', 'reflectivity', 'total_power', 'Refl', 'corrected_reflectivity', 'DBZH']
    for key in potential_name:
        try:
            rd = radar.fields[key]
            return key
        except KeyError:
            continue

    return None


def get_azimuth_resolution(azimuth):
    ra = np.diff(azimuth)
    ra[(ra < 0) | (ra > 10)] = np.NaN
    rslt = np.round(np.nanmean(ra), 1)
    return rslt


def read_radar(infile, attenuation_correction=True, reflec_offset=0):
    '''
    READ_RADAR
    Read a radar data file, will correct the attenuation if kindly asked.
    Returns a dictionnary containing the necessary parameters.
    '''

    if ".h5" in infile or ".H5" in infile:
        try:
            radar = pyart.aux_io.read_odim_h5(infile)
        except TypeError:
            print("This file is corrupted. It has not been properly created and cannot be read {}.".format(infile))
            traceback.print_exc()
            return None
    else:
        radar = pyart.io.read(infile)

    dtime_radar = netCDF4.num2date(radar.time['data'][0], radar.time['units'])

    refl_field_name = get_reflectivity_field_name(radar)

    if refl_field_name is None:
        print('Reflectivity field not found.')
        return None

    rg = radar.range['data']  # Extract range
    ngate = radar.ngates
    ntilt = radar.nsweeps
    nrays = radar.nrays

    sweep_number = radar.sweep_number['data']  # Extract number of tilt

    if (nrays / ntilt).is_integer():
        # Checking that the number of rays does not change with azimuth.
        # If it does we will have to create (or remove) the missing
        # (or the extra) dimension.
        nbeam = int(nrays / ntilt)
    else:
        azi = radar.azimuth['data'][radar.get_slice(0)]
        res_azi = get_azimuth_resolution(azi)
        nbeam = int(360 / res_azi)

    el = np.zeros((nbeam, ntilt))
    az = np.zeros((nbeam, ntilt))

    reflec = np.zeros((ngate, nbeam, ntilt))
    elevation = np.zeros((ntilt, ))

    if sweep_number[0] == 1:  # New pyart version causing problems ?
        sweep_number = sweep_number - 1

    # Radar fields have [range, time] dimensions, we want [range, azimuth, elevation]
    for cnt, sw in enumerate(sweep_number):
        sweep_slice = radar.get_slice(sw)  # Get indices of given slice

        raw_azi = radar.azimuth['data'][sweep_slice]
        raw_azi[raw_azi == 360] = 0
        raw_elev = radar.elevation['data'][sweep_slice]

        elevation[cnt] = np.mean(radar.elevation['data'][sweep_slice])

        # Extracting reflectivity
        refl_slice = radar.fields[refl_field_name]['data'][sweep_slice]

        _, uniq_index = np.unique(raw_azi, return_index=True)

        azi = raw_azi[uniq_index]
        elev = raw_elev[uniq_index]

        relf_sort_uniq_slice = refl_slice[uniq_index, :]  # Shape  (azi, r)

        if len(azi) > nbeam:
            # In case that the number of rays change with the elevation
            val_pos = range(0, nbeam)
            elev = elev[val_pos]
            azi = azi[val_pos]
            relf_sort_uniq_slice = relf_sort_uniq_slice[val_pos, :]

        while len(azi) < nbeam:
            # Idem
            azi = np.append(azi, np.NaN)
            elev = np.append(elev, elev[0])
            tmp_refl = copy.deepcopy(relf_sort_uniq_slice)
            empty_arr = np.zeros((ngate, ))
            tmp_refl = np.vstack((tmp_refl, empty_arr))

            relf_sort_uniq_slice = tmp_refl

        el[:, cnt] = elev
        az[:, cnt] = azi

        reflec[:, :, sw] = relf_sort_uniq_slice.T  # Shape  (r, azi, elev)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        reflec[reflec >= 100] = np.NaN  # NaNing the weird values
        reflec[reflec <= -20] = np.NaN

    # Applying the offset
    reflec = reflec + reflec_offset

    # Make 3D matrices for coordinates shape (r, azi, elev)
    rg2d = np.repeat(rg[:, np.newaxis], nbeam, axis=1)
    rg3d = np.repeat(rg2d[:, :, np.newaxis], ntilt, axis=2)
    az3d = np.repeat(az[np.newaxis, :, :], ngate, axis=0)
    el3d = np.repeat(el[np.newaxis, :, :], ngate, axis=0)

    data_dict = dict()
    data_dict = {'ngate': ngate, 'nbeam': nbeam, 'ntilt': ntilt, 'azang': az3d,
                 'elev_3d': el3d, 'range': rg3d, 'elang': elevation,
                 'dr': rg[2] - rg[1], 'reflec': reflec, 'time': dtime_radar}

    return data_dict  # Type: Dict
