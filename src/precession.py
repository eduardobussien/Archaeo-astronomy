"""
Long-term apparent positions of stars for archaeoastronomy.

Precession follows Vondrák, Capitaine & Wallace (2011, A&A 534, A22), valid to
+/-200,000 years. Earth rotation is referred to the Celestial Intermediate
Origin, so the hour angle comes from the Earth Rotation Angle (linear in UT1)
and the CIO locator s, which is integrated numerically from the long-term pole
because the IAU 2006 series for s diverges beyond a few millennia.

Positions are geometric mean places: nutation (<20"), annual aberration (<21"),
polar motion and atmospheric refraction are not applied.
"""
import functools
import warnings

import erfa
import numpy as np

J2000_JD = 2451545.0
EPOCH_RANGE = (-20000.0, 6000.0)

_MAS_TO_RAD = np.pi / 648_000_000.0
_S_GRID_STEP_YR = 10.0


def _pole(epj):
    return np.asarray(erfa.ltpb(epj))[..., 2, :]


@functools.lru_cache(maxsize=1)
def _cio_locator_table():
    epochs = np.arange(EPOCH_RANGE[0], EPOCH_RANGE[1] + _S_GRID_STEP_YR, _S_GRID_STEP_YR)
    x, y, z = _pole(epochs).T
    xm, ym, zm = ((a[1:] + a[:-1]) / 2.0 for a in (x, y, z))
    ds = -(xm * np.diff(y) - ym * np.diff(x)) / (1.0 + zm)
    s = np.concatenate(([0.0], np.cumsum(ds)))
    s -= s[np.searchsorted(epochs, 2000.0)]
    return epochs, s


def julian_epoch(jd):
    epj = erfa.epj(jd, 0.0)
    if np.any((epj < EPOCH_RANGE[0]) | (epj > EPOCH_RANGE[1])):
        raise ValueError(f"Epoch outside supported range {EPOCH_RANGE}")
    return epj


def cio_locator(epj):
    """CIO locator s (radians) for the mean pole at Julian epoch epj."""
    epochs, s = _cio_locator_table()
    return np.interp(epj, epochs, s)


def local_mean_sidereal_time(jd_ut1, lon):
    """Local mean sidereal time in degrees, consistent with the long-term precession model."""
    epj = julian_epoch(jd_ut1)
    gmst = erfa.era00(jd_ut1, 0.0) - erfa.eors(erfa.ltpb(epj), cio_locator(epj))
    return np.degrees(gmst + np.radians(lon)) % 360.0


def gcrs_to_altaz(vectors, jd_ut1, lat, lon):
    """
    Rotate GCRS unit vectors (..., 3) to topocentric (altitude, azimuth) in degrees.

    Azimuth is measured from north through east. The UT1 date also serves as the
    TT argument of precession; the resulting error is below 1" at 10,000 BC.
    """
    epj = julian_epoch(jd_ut1)
    x, y, _ = np.moveaxis(_pole(epj), -1, 0)
    c2i = erfa.c2ixys(x, y, cio_locator(epj))
    ra_cirs, dec_cirs = erfa.c2s(np.einsum('...ij,...j->...i', c2i, vectors))
    ha = erfa.era00(jd_ut1, 0.0) + np.radians(lon) - ra_cirs
    az, el = erfa.hd2ae(ha, dec_cirs, np.radians(lat))
    return np.degrees(el), np.degrees(az)


def star_altaz(ra, dec, pm_ra_cosdec, pm_dec, distance, jd_ut1, lat, lon):
    """
    Altitude and azimuth (degrees) of catalog stars at a UT1 Julian Date.

    ra, dec        ICRS J2000 position, degrees
    pm_ra_cosdec   proper motion in RA including the cos(dec) factor, mas/yr
    pm_dec         proper motion in declination, mas/yr
    distance       parsecs

    Star parameters and jd_ut1 broadcast against each other, so a column of
    dates against a row of stars yields a (dates, stars) grid. Space motion is
    propagated rigorously in 3D with
    zero radial velocity, matching astropy's SkyCoord.apply_space_motion.
    """
    ra, dec = np.radians(ra), np.radians(dec)
    pmr = np.asarray(pm_ra_cosdec) * _MAS_TO_RAD / np.cos(dec)
    pmd = np.asarray(pm_dec) * _MAS_TO_RAD
    parallax = 1.0 / np.asarray(distance)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', erfa.ErfaWarning)
        ra_t, dec_t, *_ = erfa.pmsafe(ra, dec, pmr, pmd, parallax, 0.0,
                                      J2000_JD, 0.0, jd_ut1, 0.0)
    return gcrs_to_altaz(erfa.s2c(ra_t, dec_t), jd_ut1, lat, lon)


def ecliptic_altaz(jd_ut1, lat, lon, longitudes):
    """Altitude and azimuth (degrees) of points on the mean ecliptic of date."""
    lam = np.radians(longitudes)
    ecl = np.stack([np.cos(lam), np.sin(lam), np.zeros_like(lam)], axis=-1)
    gcrs = np.einsum('...ji,...j->...i', erfa.ltecm(julian_epoch(jd_ut1)), ecl)
    return gcrs_to_altaz(gcrs, jd_ut1, lat, lon)
