import traceback

# Other libraries.
import pyart
import netCDF4
import numpy as np

from scipy.interpolate import griddata


def _transform_grid(r, azi, dbz):
    # Sort reflectivity field such as dbz[0, :] corresponds to azimuth 0 degree.
    pos = np.argsort(azi)
    azi = azi[pos]
    dbz = dbz[pos, :]
    if len(azi) == 360:
        return dbz

    # The azimuth length is not equal to 360, we'll have to interpolate it.
    yazi = np.arange(0, 360)
    XR, YA = np.meshgrid(r, yazi)
    Xi, Yi = np.meshgrid(r, azi)
    vq = griddata((Xi.flatten(), Yi.flatten()), dbz.flatten(),
                  (XR.flatten(), YA.flatten())).reshape(YA.shape)

    return vq


def transform_reflectivity(radar, refl_name):
    """
    Transform reflectivity from a 2D array of dimensions (time, range) to (range, azimuth, elevation)

    Parameters:
    ===========
    radar: Py-ART obj
    refl_name: str
        Name of the reflectivity field.

    Returns:
    ========
    reflectivity_3D: ndarray <range, azimuth, elevation.
        Return a 3D reflectivity array. Azimuth length is 360.
    """
    ngate = radar.ngates
    ntilt = radar.nsweeps

    reflectivity_in = radar.fields[refl_name]['data'].filled(np.NaN).copy()
    reflectivity_3D = np.zeros(ngate, 360, ntilt)

    r = radar.range['data']
    azimuth = radar.azimuth['data']
    for idx in range(ntilt):
        sl = radar.get_slice(idx)
        azi_slice = azimuth[sl]
        refl_slice = reflectivity_in[sl]
        tmp = _transform_grid(r, azi_slice, refl_slice)
        reflectivity_3D[:, :, idx] = tmp.T  # From (azi, r) to (r, azi,)

    return reflectivity_3D


def _read_radar_pyart(filename):
    try:
        radar = pyart.io.read(filename)
    except Exception:
        radar = pyart.aux_io.read_odim_h5(filename)

    return radar


def read_radar(filename, refl_name='reflectivity', offset=None):
    """
    Read input ground radar files and format range, azimuth, elevation and
    the reflectivity field the way the volume matching code wants it.

    Parameters:
    ===========
    filename: str
        Radar file name.
    refl_name: str
        Reflectivity field name.
    offset: float
        Not used anymore. Kept for compatibility issue

    Returns:
    ========
    data_dict: dictionnary
        Data structure.
    """
    radar = _read_radar_pyart(filename)
    dtime_radar = netCDF4.num2date(radar.time['data'][0], radar.time['units'])

    try:
        radar.fields[refl_name]
    except KeyError:
        raise KeyError(f"Wrong reflectivity field name provided. You said it is '{refl_name}', but I couldn't find it.")

    ngate = radar.ngates
    nbeam = 360
    ntilt = radar.nsweeps

    reflectivity = transform_reflectivity(radar, refl_name)
    reflectivity[(reflectivity < 0) | (reflectivity > 60)] = np.NaN  # Remove extreme values.

    r = radar.range['data']
    azimuth = np.arange(0, 360)
    elevation = np.zeros(ntilt)
    for idx in range(ntilt):
        sl = radar.get_slice(idx)
        elevation[idx] = radar.elevation['data'][sl].mean()

    # Make 3D matrices for coordinates shape (r, azi, elev)
    rg2d = np.repeat(r[:, np.newaxis], nbeam, axis=1)
    rg3d = np.repeat(rg2d[:, :, np.newaxis], ntilt, axis=2)
    # Azimuth
    az2d = np.repeat(azimuth[:, np.newaxis], ntilt, axis=1)
    az3d = np.repeat(az2d[np.newaxis, :, :], ngate, axis=0)
    # Elevation
    el2d = np.repeat(elevation[np.newaxis, :], nbeam, axis=0)
    el3d = np.repeat(el2d[np.newaxis, :, :], ngate, axis=0)

    data_dict = dict()
    data_dict['ngate'] = ngate
    data_dict['nbeam'] = nbeam
    data_dict['ntilt'] = ntilt
    data_dict['range'] = rg3d
    data_dict['azang'] = az3d
    data_dict['elev_3d'] = el3d
    data_dict['elang'] = elevation
    data_dict['dr'] = r[1] - r[0]
    data_dict['reflec'] = reflectivity
    data_dict['time'] = dtime_radar

    return data_dict
