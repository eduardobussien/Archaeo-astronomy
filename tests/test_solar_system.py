import astropy.units as u
import erfa
import numpy as np
import pytest
from astropy.coordinates import AltAz, EarthLocation, get_body, solar_system_ephemeris
from astropy.time import Time

from alignment import _date_to_jd
from solar_system import body_altaz, geocentric_positions
from timescales import _SEGMENTS, delta_t

GIZA_LAT, GIZA_LON = 29.9792, 31.1342

# Geocentric astrometric RA/Dec (degrees) from JPL Horizons, ephemeris DE441, at TT Julian Dates.
JPL_DE441 = [
    (-1565679.0, 1.5, {  # 20 March 9000 BC
        'Sun': (148.32850, 14.14853),
        'Moon': (40.59942, 19.83250),
        'Mercury': (161.23398, 10.98304),
        'Venus': (148.40897, 14.44918),
        'Mars': (130.02200, 19.94565),
        'Jupiter': (314.25105, -18.36220),
        'Saturn': (265.23510, -21.91184),
    }),
    (260534.0, 0.4, {  # 20 March 4000 BC
        'Sun': (79.25957, 23.83762),
        'Moon': (150.62693, 11.32502),
        'Mercury': (70.61136, 22.51952),
        'Venus': (86.77750, 28.24621),
        'Mars': (182.87011, 1.17240),
        'Jupiter': (137.32537, 17.60939),
        'Saturn': (181.85012, 1.35323),
    }),
    (625776.0, 0.25, {  # 20 March 3000 BC
        'Sun': (64.68424, 21.98877),
        'Moon': (228.85273, -19.52420),
        'Mercury': (50.94370, 17.59536),
        'Venus': (62.33022, 20.65478),
        'Mars': (98.15069, 24.66067),
        'Jupiter': (242.94424, -20.06569),
        'Saturn': (164.35565, 8.40389),
    }),
    (1356261.0, 0.25, {  # 20 March 1000 BC
        'Sun': (37.72376, 15.11602),
        'Moon': (34.82652, 13.31045),
        'Mercury': (18.75312, 5.39379),
        'Venus': (35.26428, 13.02779),
        'Mars': (218.95196, -12.84584),
        'Jupiter': (102.39855, 23.32857),
        'Saturn': (122.98717, 20.25363),
    }),
]


def _separation_deg(alt1, az1, alt2, az2):
    alt1, az1, alt2, az2 = map(np.radians, (alt1, az1, alt2, az2))
    cos_sep = (np.sin(alt1) * np.sin(alt2)
               + np.cos(alt1) * np.cos(alt2) * np.cos(az1 - az2))
    return np.degrees(np.arccos(np.clip(cos_sep, -1.0, 1.0)))


@pytest.mark.parametrize('year, expected', [
    (-1000, 25400), (-500, 17190), (0, 10580), (500, 5710), (1000, 1570),
    (1500, 200), (1800, 14), (1900, -3), (1950, 29), (2000, 64),
])
def test_delta_t_matches_nasa_table(year, expected):
    # Published values are rounded to 10 s before 1600 and to 1 s after.
    assert abs(delta_t(year) - expected) < max(2.0, 0.005 * abs(expected))


@pytest.mark.parametrize('boundary', [seg[0] for seg in _SEGMENTS[1:]])
def test_delta_t_is_continuous_at_segment_boundaries(boundary):
    assert abs(delta_t(boundary - 1e-6) - delta_t(boundary + 1e-6)) < 2.0


@pytest.mark.parametrize('jd_tt, tolerance_deg, reference', JPL_DE441)
def test_geocentric_positions_match_jpl_de441(jd_tt, tolerance_deg, reference):
    positions = geocentric_positions(jd_tt)
    for name, (ra, dec) in reference.items():
        ours_ra, ours_dec = erfa.c2s(positions[name])
        sep = np.degrees(erfa.seps(ours_ra, ours_dec, np.radians(ra), np.radians(dec)))
        assert sep < tolerance_deg, f"{name}: {sep:.3f} deg"


def test_topocentric_bodies_match_astropy_near_j2000():
    solar_system_ephemeris.set('builtin')
    jd = _date_to_jd(2000, 3, 20, 20.0)
    t = Time(jd, format='jd', scale='ut1')
    t.delta_ut1_utc = 0.0
    frame = AltAz(obstime=t, location=EarthLocation(lat=GIZA_LAT * u.deg, lon=GIZA_LON * u.deg))

    for name, (alt, az) in body_altaz(jd, GIZA_LAT, GIZA_LON).items():
        ref = get_body(name.lower(), t, location=frame.location).transform_to(frame)
        assert _separation_deg(alt, az, ref.alt.deg, ref.az.deg) < 1.0 / 60.0, name


def test_bodies_available_at_10500_bc():
    bodies = body_altaz(_date_to_jd(-10499, 3, 20, 20.0), GIZA_LAT, GIZA_LON)
    assert set(bodies) == {'Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn'}
    assert all(np.isfinite(v).all() for v in bodies.values())
