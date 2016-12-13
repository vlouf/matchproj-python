import pyart
import numpy as np


def populate_missing_azimuth(azi, refl_slice, ngate):
    '''POPULATE_MISSING_AZIMUTH'''

    a = azi.tolist()
    tmp_refl = refl_slice
    for e in range(0, 360):
        if a.count(e) == 0:
            azi = np.append(azi, e)
            empty_arr = np.zeros((ngate, ))
            tmp_refl = np.vstack((tmp_refl, empty_arr))

    return azi, tmp_refl


def read_radar(infile):
    '''READ_RADAR'''

    try:
        radar = pyart.io.read(infile)
    except KeyError:
        radar = pyart.aux_io.read_odim_h5(infile)

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
        try:
            refl_slice = radar.fields['DBZ']['data'][sweep_slice]  # Reflectivity
        except KeyError:
            refl_slice = radar.fields['DBZ_F']['data'][sweep_slice]  # Reflectivity
        except KeyError:
            refl_slice = radar.fields['reflectivity']['data'][sweep_slice]  # Reflectivity

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
