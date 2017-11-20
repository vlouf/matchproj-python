# Python Standard Library
import copy
import warnings
import traceback

# Other libraries.
import pyart
import numpy as np

from scipy.integrate import cumtrapz


def correct_attenuation(radar, refl_field_name='DBZ_F', rhv_field_name='RHOHV_F', phidp_field_name='PHIDP_F'):
    """
    Use of pyart correction capabilities to estimate and correct the attenuation.

    Parameters:
    ===========
        radar: struct
            Py-ART radar structure
        method: str
            Can be 'pyart' or 'bringi'
    """

    print("Correcting ground radar attenuation.")

    try:
        # Tries if normalised coherent poiwer field exist
        ncp = radar.fields['NCP']
    except KeyError:
        # Creates a dummy NCP field.
        reflec = radar.fields[refl_field_name]['data']
        ncp = np.zeros_like(reflec) + 1
        radar.add_field_like(refl_field_name, 'NCP', ncp)

        spec_at, cor_z = pyart.correct.calculate_attenuation(radar, 0,
                                                             refl_field=refl_field_name,
                                                             ncp_field='NCP',
                                                             rhv_field=rhv_field_name,
                                                             phidp_field=phidp_field_name)

    return cor_z


def do_gatefilter(radar, dbz_name, zdr_name, rhohv_name):
    gf = pyart.filters.GateFilter(radar)

    gf.exclude_outside(dbz_name, 0, 80)
    gf.exclude_outside(zdr_name, -3.0, 8.0)
    gf.exclude_below(rhohv_name, 0.5)
    gf.include_above(rhohv_name, 0.8)

    gf_despeckeld = pyart.correct.despeckle_field(radar, dbz_name, gatefilter=gf)

    return gf_despeckeld


def filter_hardcoding(my_array, nuke_filter, bad=-9999):
    """
    Harcoding GateFilter into an array.

    Parameters:
    ===========
        my_array: array
            Array we want to clean out.
        nuke_filter: gatefilter
            Filter we want to apply to the data.
        bad: float
            Fill value.

    Returns:
    ========
        to_return: masked array
            Same as my_array but with all data corresponding to a gate filter
            excluded.
    """
    filt_array = np.ma.masked_where(nuke_filter.gate_excluded, my_array)
    filt_array.set_fill_value(bad)
    filt_array = filt_array.filled(fill_value=bad)
    to_return = np.ma.masked_where(filt_array == bad, filt_array)
    return to_return


def get_reflectivity_field_name(radar):
    '''
    GET_REFLECTIVITY_FIELD_NAME
    Because of different conventions for naming fields, it will try a variety
    of different reflectvity name and return the first one that works
    '''

    potential_name = ['DBZ_F', 'DBZ', 'reflectivity', 'total_power', 'Refl', 'corrected_reflectivity']
    for pn in potential_name:
        try:
            rd = radar.fields[pn]
            return pn
        except KeyError:
            continue

    return None


def get_phidb_field_name(radar):
    '''
    GET_PHIDB_FIELD_NAME
    Because of different conventions for naming fields, it will try a variety
    of different reflectvity name and return the first one that works
    '''

    potential_name = ['PHIDP_F', 'PHIDP', 'differential_phase']
    for pn in potential_name:
        try:
            rd = radar.fields[pn]
            return pn
        except KeyError:
            continue

    return None


def get_rhohv_field_name(radar):
    '''
    GET_RHOHV_FIELD_NAME
    Because of different conventions for naming fields, it will try a variety
    of different reflectvity name and return the first one that works
    '''

    potential_name = ['RHOHV_F', 'RHOHV', 'cross_correlation_ratio']
    for pn in potential_name:
        try:
            rd = radar.fields[pn]
            return pn
        except KeyError:
            continue

    return None


def get_zdr_field_name(radar):
    '''
    GET_RHOHV_FIELD_NAME
    Because of different conventions for naming fields, it will try a variety
    of different reflectvity name and return the first one that works
    '''

    potential_name = ['ZDR', 'ZDR_F', 'differential_reflectivity']
    for pn in potential_name:
        try:
            rd = radar.fields[pn]
            return pn
        except KeyError:
            continue

    return None


def get_azimuth_resolution(azimuth):
    ra = np.diff(azimuth)
    ra[(ra < 0) | (ra > 10)] = np.NaN
    rslt = np.round(np.nanmean(ra), 1)
    return rslt


def populate_missing_azimuth(azi, refl_slice, ngate):
    '''
    POPULATE_MISSING_AZIMUTH
    If the number of azimuth of one sweep is lower than 360, create empty
    columns corresponding to the missing azimuth.

    Returns the new azimuth and the reflectvity field
    '''

    a = azi.tolist()
    tmp_refl = refl_slice
    for e in range(0, 360):
        if a.count(e) == 0:
            azi = np.append(azi, e)
            empty_arr = np.zeros((ngate, ))
            tmp_refl = np.vstack((tmp_refl, empty_arr))

    return azi, tmp_refl


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

    refl_field_name = get_reflectivity_field_name(radar)
    # zdr_name = get_zdr_field_name(radar)
    # phidp_name = get_phidb_field_name(radar)
    # rhohv_name = get_rhohv_field_name(radar)
    if refl_field_name is None:
        print('Reflectivity field not found.')
        return None

    # if zdr_name is not None:
    #     try:
    #         gf = do_gatefilter(radar, refl_field_name, zdr_name, rhohv_name)
    #         dbz = filter_hardcoding(radar.fields[refl_field_name]['data'], gf)
    #         radar.add_field_like(refl_field_name, "NEW_DBZ", dbz, replace_existing=True)
    #         refl_field_name = "NEW_DBZ"
    #     except Exception:
    #         pass

    # if attenuation_correction:
    #     if phidp_name is None or rhohv_name is None or kdp_name is None:
    #         attenuation_correction = False
    #         print("Attenuation correction impossible, missing dualpol field.")
    #     else:
    #         radar = correct_attenuation(radar, refl_field_name=refl_field_name,
    #                                     rhv_field_name=rhohv_name, phidp_field_name=phidp_name)
    #
    #         radar.add_field('corrected_reflectivity_horizontal', cor_z)
    #         refl_field_name = "corrected_reflectivity_horizontal"

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
    data_dict = {'ngate': ngate,        # Number of gates
                 'nbeam': nbeam,        # Number of azimuth by elevation
                 'ntilt': ntilt,        # Number of sweeps
                 'azang': az3d,         # Azimuth
                 'elev_3d': el3d,       # Elevation angle
                 'range': rg3d,         # Radar's range
                 'elang': elevation,    # Radar's elevation
                 'dr': rg[2] - rg[1],     # Gate spacing
                 'reflec': reflec}      # Reflectivity as shape of (r, azi, ele)

    return data_dict  # Type: Dict
