import pyart
import numpy as np
import copy
from ..util_fun import print_yellow


def correct_attenuation(radar, refl_field_name='DBZ_F',
                        rhv_field_name='RHOHV_F', phidp_field_name='PHIDP_F'):
    """
    CORRECT_ATTENUATION
    Use of pyart correction capabilities to estimate the attenuation and
    correct it
    Returns a pyart radar structure.
    """

    the_radar = copy.deepcopy(radar)
    print_yellow("Correcting ground radar attenuation.")

    try:
        # Tries if normalised coherent poiwer field exist
        ncp = radar.fields['NCP']
    except KeyError:
        # Creates a dummy NCP field.
        reflec = the_radar.fields[refl_field_name]['data']
        ncp = np.zeros_like(reflec) + 1
        the_radar.add_field_like(refl_field_name, 'NCP', ncp)

    spec_at, cor_z = pyart.correct.calculate_attenuation(the_radar, 0,
                     refl_field=refl_field_name, ncp_field='NCP',
                     rhv_field=rhv_field_name, phidp_field=phidp_field_name)

    the_radar.add_field('specific_attenuation', spec_at)
    the_radar.add_field('corrected_reflectivity_horizontal', cor_z)

    return the_radar


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


def what_is_the_reflectivity_field_name(radar):
    '''
    WHAT_IS_THE_REFLECTIVITY_FIELD_NAME
    Because of different conventions for naming fields, it will try a variety
    of different reflectvity name and return the first one that works
    '''

    potential_name = ['DBZ_F', 'DBZ', 'reflectivity']
    for pn in potential_name:
        try:
            rd = radar.fields[pn]
            return pn
        except KeyError:
            continue

    return None


def read_radar(infile, attenuation_correction=True):
    '''
    READ_RADAR
    Read a radar data file, will correct the attenuation if kindly asked.
    Returns a dictionnary containing the necessary parameters.
    '''

    try:
        radar = pyart.io.read(infile)
    except KeyError:
        radar = pyart.aux_io.read_odim_h5(infile)

    refl_field_name = what_is_the_reflectivity_field_name(radar)

    if attenuation_correction:
        radar = correct_attenuation(radar, refl_field_name=refl_field_name)

    rg = radar.range['data']  # Extract range
    sweep_number = radar.sweep_number['data']  # Extract number of tilt

    ngate = radar.ngates
    ntilt = radar.nsweeps
    nbeam = 360

    elevation = np.zeros((len(sweep_number), ))

    # Allocate reflectivity array with dimensions (range, azimuth, elevation)
    reflec = np.zeros((ngate, nbeam, ntilt))

    # Rearranging arrays
    for cnt, sw in enumerate(sweep_number):
        sweep_slice = radar.get_slice(sw)  # Get indices of given slice
        azi = radar.azimuth['data'][sweep_slice].astype(int)  # Extract azimuth

        # Extracting reflectivity
        if attenuation_correction:
            refl_slice = radar.fields['corrected_reflectivity_horizontal']['data'][sweep_slice]
        else:
            refl_slice = radar.fields[refl_field_name]['data'][sweep_slice]

        elevation[cnt] = np.mean(radar.elevation['data'][sweep_slice])

        if len(azi) < 360:
            azi, refl_slice = populate_missing_azimuth(azi, refl_slice, ngate)

        _, uniq_index = np.unique(azi, return_index=True)  #  Sorted unique azimuth

        if len(uniq_index) != 360:
            raise ValueError('Wrong size.')

        relf_sort_uniq_slice = refl_slice[uniq_index, :]
        reflec[:, :, sw] = relf_sort_uniq_slice.T

        # print(azi.shape, azi[uniq_index].shape)
        # print(azi[uniq_index][0], azi[uniq_index][-1])

    azimuth = np.arange(0, 360, dtype=int)

    to_return = dict()
    to_return = {'ngate': ngate,        # Number of gates
                 'nbeam': nbeam,        # Number of azimuth by elevation
                 'ntilt': ntilt,        # Number of sweeps
                 'azang': azimuth,      # Azimuth
                 'elang': elevation,    # Elevation angle
                 'range': rg,           # Radar's range
                 'dr': rg[2]-rg[1],     # Gate spacing
                 'reflec': reflec}      # Reflectivity as shape of (r, azi, ele)

    return to_return  # Type: Dict
