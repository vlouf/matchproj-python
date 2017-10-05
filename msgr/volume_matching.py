# PSL
import itertools

# Other libraries
import numpy as np

from numpy import sin, cos, tan, sqrt, exp, pi


def process(satellite, cpol, nprof, ntilt, reflectivity_satellite, reflectivity_ground_radar,
            refp_ss, refp_sh, refg_ku, xp, yp, zp, rt, ep, rg, xg, yg, elang, dr, alpha, l_dbz=True):

    zt = satellite.altitude
    bwt = satellite.beamwidth
    drt = satellite.dr

    bwr = cpol.beamwidth
    earth_gaussian_radius = cpol.gaussian_radius
    rmax = cpol.rmax

    reflectivity_satellite = reflectivity_satellite.filled(np.NaN)
    reflectivity_ground_radar = reflectivity_ground_radar.filled(np.NaN)

    # Create arrays to store comparison variables
    '''Coordinates'''
    x = np.zeros((nprof, ntilt))  # x coordinate of sample
    y = np.zeros((nprof, ntilt))  # y coordinate of sample
    z = np.zeros((nprof, ntilt))  # z coordinate of sample
    dz = np.zeros((nprof, ntilt))  # depth of sample
    ds = np.zeros((nprof, ntilt))  # width of sample
    r = np.zeros((nprof, ntilt))  # range of sample from ground radar

    '''Reflectivities'''
    ref1 = np.zeros((nprof, ntilt)) + np.NaN  # PR reflectivity
    ref2 = np.zeros((nprof, ntilt)) + np.NaN  # PR reflec S-band, snow
    ref3 = np.zeros((nprof, ntilt)) + np.NaN  # PR reflec S-band, hail
    ref4 = np.zeros((nprof, ntilt)) + np.NaN  # GR reflectivity
    ref5 = np.zeros((nprof, ntilt)) + np.NaN  # GR reflectivity Ku-band
    iref1 = np.zeros((nprof, ntilt)) + np.NaN  # path-integrated PR reflec
    iref2 = np.zeros((nprof, ntilt)) + np.NaN  # path-integrated GR reflec
    stdv1 = np.zeros((nprof, ntilt)) + np.NaN  # STD of PR reflectivity
    stdv2 = np.zeros((nprof, ntilt)) + np.NaN  # STD of GR reflectivity

    '''Number of bins in sample'''
    ntot1 = np.zeros((nprof, ntilt), dtype=int)  # Total nb of PR bin in sample
    # Nb of rejected PR bin in sample
    nrej1 = np.zeros((nprof, ntilt), dtype=int)
    ntot2 = np.zeros((nprof, ntilt), dtype=int)  # Total nb of GR bin in sample
    # Nb of rejected GR bin in sample
    nrej2 = np.zeros((nprof, ntilt), dtype=int)
    # Total volume of PR bins in sample
    vol1 = np.zeros((nprof, ntilt)) + np.NaN
    # Total volume of GR bins in sample
    vol2 = np.zeros((nprof, ntilt)) + np.NaN

    # Compute the volume of each radar (sat/ground) bin
    volp = 1.e-9 * np.pi * drt * (rt * np.pi / 180 * bwt / 2.)**2
    volg = 1e-9 * np.pi * dr * (rg * np.pi / 180 * bwr / 2)**2

    # Compute the path-integrated reflectivities at every points
    nat_refp = 10**(reflectivity_satellite / 10.0)  # In natural units
    nat_refg = 10**(reflectivity_ground_radar / 10.0)
    irefp = np.fliplr(np.nancumsum(np.fliplr(nat_refp), axis=1))
    irefg = np.nancumsum(nat_refg, axis=0)
    irefp = drt * (irefp - nat_refp / 2)
    irefg = dr * (irefg - nat_refg / 2)
    irefp = 10 * np.log10(irefp)
    irefg = 10 * np.log10(irefg)

    # Convert to linear units
    if not l_dbz:
        reflectivity_satellite = 10**(reflectivity_satellite / 10.0)
        reflectivity_ground_radar = 10**(reflectivity_ground_radar / 10.0)
        refp_ss = 10**(refp_ss / 10.0)
        refp_sh = 10**(refp_sh / 10.0)
        refg_ku = 10**(refg_ku / 10.0)

    irefp = 10**(irefp / 10.0)
    irefg = 10**(irefg / 10.0)

    # Loop over the TRMM/GPM profiles and Loop over the GR elevation scan
    for ii, jj in itertools.product(range(nprof), range(ntilt)):
        # Identify those PR bins which fall within the GR sweep
        ip = np.where((ep[ii, :] >= elang[jj] - bwr / 2) &
                      (ep[ii, :] <= elang[jj] + bwr / 2))

        # Store the number of bins
        ntot1[ii, jj] = len(ip)
        if len(ip) == 0:
            continue

        x[ii, jj] = np.mean(xp[ii, ip])
        y[ii, jj] = np.mean(yp[ii, ip])
        z[ii, jj] = np.mean(zp[ii, ip])

        # Compute the thickness of the layer
        nip = len(ip)
        dz[ii, jj] = nip * drt * cos(pi / 180 * alpha[ii, 0])

        # Compute the PR averaging volume
        vol1[ii, jj] = np.sum(volp[ii, ip])

        # Note the mean TRMM beam diameter
        ds[ii, jj] = pi / 180 * bwt * \
            np.mean((zt - zp[ii, ip]) / cos(pi / 180 * alpha[ii, ip]))

        # Note the radar range
        s = sqrt(x[ii, jj]**2 + y[ii, jj]**2)
        r[ii, jj] = (earth_gaussian_radius + z[ii, jj]) * \
            sin(s / earth_gaussian_radius) / cos(pi / 180 * elang[jj])

        # Check that sample is within radar range
        if r[ii, jj] + ds[ii, jj] / 2 > rmax:
            continue

        # Extract the relevant PR data
        refp1 = reflectivity_satellite[ii, ip].flatten()
        refp2 = refp_ss[ii, ip].flatten()
        refp3 = refp_sh[ii, ip].flatten()
        irefp1 = irefp[ii, ip].flatten()

        # Average over those bins that exceed the reflectivity
        # threshold (linear average)

        ref1[ii, jj] = np.nanmean(refp1)
        ref3[ii, jj] = np.nanmean(refp2)
        ref4[ii, jj] = np.nanmean(refp3)
        iref1[ii, jj] = np.nanmean(irefp1)

        if not l_dbz:
            stdv1[ii, jj] = np.nanstd(10 * np.log10(refp1))
        else:
            stdv1[ii, jj] = np.nanstd(refp1)

        # Note the number of rejected bins
        nrej1[ii, jj] = int(np.sum(np.isnan(refp1)))
        if ~np.isnan(stdv1[ii, jj]) and nip - nrej1[ii, jj] > 1:
            continue

        # Compute the horizontal distance to all the GR bins
        d = sqrt((xg[:, :, jj] - x[ii, jj])**2 +
                 (yg[:, :, jj] - y[ii, jj])**2)

        # Find all GR bins within the SR beam
        igx, igy = np.where(d <= ds[ii, jj] / 2)

        # Store the number of bins
        ntot2[ii, jj] = len(igx)
        if len(igx) == 0:
            continue

        # Extract the relevant GR data
        refg1 = reflectivity_ground_radar[:, :, jj][igx, igy].flatten()
        refg2 = refg_ku[:, :, jj][igx, igy].flatten()
        volg1 = volg[:, :, jj][igx, igy].flatten()
        irefg1 = irefg[:, :, jj][igx, igy].flatten()

        #  Comupte the GR averaging volume
        vol2[ii, jj] = np.sum(volg1)

        # Average over those bins that exceed the reflectivity
        # threshold (exponential distance and volume weighting)
        w = volg1 * exp(-1 * (d[igx, igy] / (ds[ii, jj] / 2.))**2)
        w = w * refg1 / refg2

        ref2[ii, jj] = np.nansum(w * refg1) / np.nansum(w)

        # if ref2[ii, jj] < minrefp:
        #     ref2[ii, jj] = np.NaN

        ref5[ii, jj] = np.nansum(w * refg2) / np.nansum(w)
        iref2[ii, jj] = np.nansum(w * irefg1) / np.nansum(w)

        if not l_dbz:
            stdv2[ii, jj] = np.nanstd(10 * np.log10(refg1))
        else:
            stdv2[ii, jj] = np.nanstd(refg1)

        # Note the number of rejected bins
        nrej2[ii, jj] = int(np.sum(np.isnan(refg1)))
    # END FOR (satellite profiles, radar elevation)

    # Convert back to dBZ
    iref1 = 10 * np.log10(iref1)
    iref2 = 10 * np.log10(iref2)

    if not l_dbz:
        ref1 = 10 * np.log10(ref1)
        ref2 = 10 * np.log10(ref2)
        ref3 = 10 * np.log10(ref3)
        ref4 = 10 * np.log10(ref4)
        ref5 = 10 * np.log10(ref5)

    return x, y, z, dz, ds, r, ref1, ref2, ref3, ref4, ref5, iref1, iref2, stdv1, stdv2, ntot1, nrej1, ntot2, nrej2, vol1, vol2
