"""
Sun, Moon and naked-eye planets from ERFA's analytic theories.

Earth: epv00 (VSOP2000). Planets: plan94 (Simon et al. 1994). Moon: moon98
(Meeus 1998, after ELP-2000/82). Checked against JPL DE441, these stay within
0.35 deg of the reference back to 4000 BC; see tests/test_solar_system.py.
Positions are geometric with planetary light time; the topocentric Moon
includes parallax for a spherical Earth.
"""
import warnings

import erfa
import numpy as np

from precession import gcrs_to_altaz
from timescales import ut1_to_tt

BODIES = ('Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn')

_PLANET_ID = {'Mercury': 1, 'Venus': 2, 'Mars': 4, 'Jupiter': 5, 'Saturn': 6}
_LIGHT_SPEED_AU_PER_DAY = 173.1446326846693
_EARTH_RADIUS_AU = 6378.137 / 149_597_870.7


def _earth_heliocentric(jd_tt):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', erfa.ErfaWarning)
        return erfa.epv00(jd_tt, 0.0)[0]['p']


def geocentric_positions(jd_tt):
    """Geocentric GCRS position vectors (au) of each body at a TT Julian Date."""
    earth = _earth_heliocentric(jd_tt)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', erfa.ErfaWarning)
        positions = {'Sun': -earth, 'Moon': erfa.moon98(jd_tt, 0.0)['p']}
        for name, planet_id in _PLANET_ID.items():
            emitted = jd_tt
            for _ in range(2):
                vec = erfa.plan94(emitted, 0.0, planet_id)['p'] - earth
                emitted = jd_tt - np.linalg.norm(vec, axis=-1) / _LIGHT_SPEED_AU_PER_DAY
            positions[name] = vec
    return positions


def body_altaz(jd_ut1, lat, lon):
    """Topocentric {name: (altitude, azimuth)} in degrees at a UT1 Julian Date."""
    positions = geocentric_positions(ut1_to_tt(jd_ut1))
    vectors = np.stack([positions[name] for name in BODIES])
    distance = np.linalg.norm(vectors, axis=-1)
    alt, az = gcrs_to_altaz(vectors / distance[:, None], jd_ut1, lat, lon)
    alt_rad = np.radians(alt)
    alt = np.degrees(np.arctan2(np.sin(alt_rad) - _EARTH_RADIUS_AU / distance, np.cos(alt_rad)))
    return {name: (float(a), float(z)) for name, a, z in zip(BODIES, alt, az)}


def sun_altitude(jd_ut1, lat, lon):
    """Altitude (degrees) of the Sun's centre; jd_ut1 may be an array."""
    sun = -_earth_heliocentric(ut1_to_tt(jd_ut1))
    alt, _ = gcrs_to_altaz(sun / np.linalg.norm(sun, axis=-1, keepdims=True), jd_ut1, lat, lon)
    return alt
