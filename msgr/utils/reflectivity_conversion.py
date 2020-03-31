import copy
import numpy as np


def convert_sat_refl_to_gr_band(refp, zp, zbb, bbwidth, radar_band='S'):
    """
    Convert the satellite reflectivity to S, C, or X-band using the Cao et al.
    (2013) method.

    Parameters
    ==========
    refp:
        Satellite reflectivity field.
    zp:
        Altitude.
    zbb:
        Bright band height.
    bbwidth:
        Bright band width.
    radar_band: str
        Possible values are 'S', 'C', or 'X'

    Return
    ======
    refp_ss:
        Stratiform reflectivity conversion from Ku-band to S-band
    refp_sh:
        Convective reflectivity conversion from Ku-band to S-band
    """
    # Set coefficients for conversion from Ku-band to S-band
    #        Rain      90%      80%      70%      60%      50%      40%      30%      20%      10%     Snow
    as0=[ 4.78e-2, 4.12e-2, 8.12e-2, 1.59e-1, 2.87e-1, 4.93e-1, 8.16e-1, 1.31e+0, 2.01e+0, 2.82e+0, 1.74e-1]
    as1=[ 1.23e-2, 3.66e-3, 2.00e-3, 9.42e-4, 5.29e-4, 5.96e-4, 1.22e-3, 2.11e-3, 3.34e-3, 5.33e-3, 1.35e-2]
    as2=[-3.50e-4, 1.17e-3, 1.04e-3, 8.16e-4, 6.59e-4, 5.85e-4, 6.13e-4, 7.01e-4, 8.24e-4, 1.01e-3,-1.38e-3]
    as3=[-3.30e-5,-8.08e-5,-6.44e-5,-4.97e-5,-4.15e-5,-3.89e-5,-4.15e-5,-4.58e-5,-5.06e-5,-5.78e-5, 4.74e-5]
    as4=[ 4.27e-7, 9.25e-7, 7.41e-7, 6.13e-7, 5.80e-7, 6.16e-7, 7.12e-7, 8.22e-7, 9.39e-7, 1.10e-6, 0]
    #        Rain      90%      80%      70%      60%      50%      40%      30%      20%      10%     Hail
    ah0=[ 4.78e-2, 1.80e-1, 1.95e-1, 1.88e-1, 2.36e-1, 2.70e-1, 2.98e-1, 2.85e-1, 1.75e-1, 4.30e-2, 8.80e-2]
    ah1=[ 1.23e-2,-3.73e-2,-3.83e-2,-3.29e-2,-3.46e-2,-2.94e-2,-2.10e-2,-9.96e-3,-8.05e-3,-8.27e-3, 5.39e-2]
    ah2=[-3.50e-4, 4.08e-3, 4.14e-3, 3.75e-3, 3.71e-3, 3.22e-3, 2.44e-3, 1.45e-3, 1.21e-3, 1.66e-3,-2.99e-4]
    ah3=[-3.30e-5,-1.59e-4,-1.54e-4,-1.39e-4,-1.30e-4,-1.12e-4,-8.56e-5,-5.33e-5,-4.66e-5,-7.19e-5, 1.90e-5]
    ah4=[ 4.27e-7, 1.59e-6, 1.51e-6, 1.37e-6, 1.29e-6, 1.15e-6, 9.40e-7, 6.71e-7, 6.33e-7, 9.52e-7, 0]

    refp_ss = np.zeros(refp.shape) + np.NaN  # snow
    refp_sh = np.zeros(refp.shape) + np.NaN  # hail
    zmlt = zbb + bbwidth / 2.    # APPROXIMATION!
    zmlb = zbb - bbwidth / 2.    # APPROXIMATION!
    ratio = (zp - zmlb) / (zmlt - zmlb)

    iax, iay = np.where(ratio >= 1)
    # above melting layer
    if len(iax) > 0:
        dfrs = as0[10] + as1[10] * refp[iax, iay] + as2[10] * refp[iax, iay]**2 + as3[10] * refp[iax, iay]**3 + as4[10] * refp[iax, iay]**4
        dfrh = ah0[10] + ah1[10] * refp[iax, iay] + ah2[10] * refp[iax, iay]**2 + ah3[10] * refp[iax, iay]**3 + ah4[10] * refp[iax, iay]**4
        refp_ss[iax, iay] = refp[iax, iay] + dfrs
        refp_sh[iax, iay] = refp[iax, iay] + dfrh

    ibx, iby = np.where(ratio <= 0)
    if len(ibx) > 0:  # below the melting layer
        dfrs = as0[0] + as1[0] * refp[ibx, iby] + as2[0] * refp[ibx, iby]**2 + as3[0] * refp[ibx, iby]**3 + as4[0] * refp[ibx, iby]**4
        dfrh = ah0[0] + ah1[0] * refp[ibx, iby] + ah2[0] * refp[ibx, iby]**2 + ah3[0] * refp[ibx, iby]**3 + ah4[0] * refp[ibx, iby]**4
        refp_ss[ibx, iby] = refp[ibx, iby] + dfrs
        refp_sh[ibx, iby] = refp[ibx, iby] + dfrh

    imx, imy = np.where((ratio > 0) & (ratio < 1))
    if len(imx) > 0:  # within the melting layer
        ind = np.round(ratio[imx, imy]).astype(int)[0]
        dfrs = as0[ind] + as1[ind] * refp[imx, imy] + as2[ind] * refp[imx, imy]**2 + as3[ind] * refp[imx, imy]**3 + as4[ind] * refp[imx, imy]**4
        dfrh = ah0[ind] + ah1[ind] * refp[imx, imy] + ah2[ind] * refp[imx, imy]**2 + ah3[ind] * refp[imx, imy]**3 + ah4[ind] * refp[imx, imy]**4
        refp_ss[imx, imy] = refp[imx, imy] + dfrs
        refp_sh[imx, imy] = refp[imx, imy] + dfrh

    # Jackson Tan's fix for C-band
    if radar_band == 'C':
        deltas = 5.3 / 10.0 * (refp_ss - refp)
        refp_ss = refp + deltas
        deltah = 5.3 / 10.0 * (refp_sh - refp)
        refp_sh = refp + deltah
    elif radar_band == 'X':
        deltas = 3.2 / 10.0 * (refp_ss - refp)
        refp_ss = refp + deltas
        deltah = 3.2 / 10.0 * (refp_sh - refp)
        refp_sh = refp + deltah

    return refp_ss, refp_sh


def convert_gpmrefl_grband_dfr(refp, radar_band=None):
    '''
    Convert GPM reflectivity to ground radar band using the DFR relationship
    found in Louf et al. (2019) paper.
    NOTE: It's a continuous relationship, so here there is no difference here
    between conv/strat.
    Parameters
    ==========
    refp:
        Satellite reflectivity field.
    radar_band: str
        Possible values are 'S', 'C', or 'X'
    Return
    ======
    refl_strat:
        Stratiform reflectivity conversion from Ku-band to S-band
    refl_conv:
        Convective reflectivity conversion from Ku-band to S-band
    '''
    if radar_band == 'S':
        cof = np.array([ 2.01236803e-07, -6.50694273e-06,  1.10885533e-03, -6.47985914e-02, -7.46518423e-02])
        dfr = np.poly1d(cof)
    elif radar_band == 'C':
        cof = np.array([ 1.21547932e-06, -1.23266138e-04,  6.38562875e-03, -1.52248868e-01, 5.33556919e-01])
        dfr = np.poly1d(cof)
    elif radar_band == 'X':
        # Use of C band DFR relationship multiply by ratio
        cof = np.array([ 1.21547932e-06, -1.23266138e-04,  6.38562875e-03, -1.52248868e-01, 5.33556919e-01])
        dfr = 3.2 / 5.5 * np.poly1d(cof)
    else:
         raise ValueError(f'Radar reflectivity band ({radar_band}) not supported.')

    # It's a continuous relationship, so here there is no difference between conv/strat.
    refl_strat = refp + dfr(refp)
    refl_conv = refl_strat

    return refl_strat, refl_conv


def convert_to_Ku(refg, zg, zbb, radar_band='S'):
    '''
    From Liao and Meneghini (2009)

    Parameters
    ==========
    refg:
        Ground radar reflectivity field.
    zg:
        Altitude.
    zbb:
        Bright band height.
    bbwidth:
        Bright band width.
    radar_band: str
        Possible values are 'S', 'C', or 'X'

    Returns
    =======
    refg_ku:
        Ground radar reflectivity field converted to Ku-band.
    '''

    refg_ku = np.zeros(refg.shape) + np.NaN
    iax, iay, iaz = np.where(zg >= zbb)

    #  Above bright band
    if len(iax) > 0:
        refg_ku[iax, iay, iaz] = 0.185074 + 1.01378 * refg[iax, iay, iaz] - 0.00189212 * refg[iax, iay, iaz]**2

    # Below bright band
    ibx, iby, ibz = np.where(zg < zbb)
    if len(ibx) > 0:
        refg_ku[ibx, iby, ibz] = -1.50393 + 1.07274 * refg[ibx, iby, ibz] + 0.000165393 * refg[ibx, iby, ibz]**2

    #  Jackson Tan's fix for C-band
    if radar_band == 'C':
        delta = (refg_ku - refg) * 5.3 / 10.0
        refg_ku = refg + delta
    elif radar_band == 'X':
        delta = (refg_ku - refg) * 3.2 / 10.0
        refg_ku = refg + delta

    return refg_ku
