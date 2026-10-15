import json
from pathlib import Path

import astropy.units as u
import erfa
import numpy as np
import pytest
from astropy.coordinates import AltAz, EarthLocation, get_body, solar_system_ephemeris
from astropy.time import Time

from alignment import _date_to_jd
from precession import _to_cirs
from solar_system import _sun_geocentric, body_altaz, geocentric_positions
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


# Largest declination of the Sun through the year, from JPL Horizons (DE441) every 6 hours.
@pytest.mark.parametrize('year, de441_max', [(-8999, 24.205), (-2499, 23.977)])
def test_sun_solstice_declination_matches_de441(year, de441_max):
    # Without placing the Sun on the long-term ecliptic, 9000 BC came out 0.19 deg short.
    jd = _date_to_jd(year, 1, 1, 0.0) + np.arange(0, 365, 0.25)
    sun = _sun_geocentric(jd)
    dec = np.degrees(_to_cirs(sun / np.linalg.norm(sun, axis=-1, keepdims=True), jd)[1])
    assert abs(dec.max() - de441_max) < 0.02
    assert abs(-dec.min() - de441_max) < 0.02


# The Sun's geocentric astrometric RA/Dec (degrees) from JPL Horizons (DE441) at 0h TT on the
# first of each month (proleptic Gregorian), from 9998 BC (Horizons' Sun starts in 9999 BC) to
# AD 2000. Horizons was queried with TLIST_TYPE=JD: without it, it read the older dates as
# Modified Julian Dates. Each era's tolerance is about 1.2 times the error measured when this
# test was written, so any loss of accuracy fails it.
SUN_DE441 = {
    -9997: (0.35, [
        (82.75146, 24.67795), (115.11430, 23.09285), (142.96631, 16.18841), (171.73763, 4.55253),
        (199.20737, -8.00818), (229.72698, -19.09976), (262.39561, -24.65642), (296.92586, -22.79733),
        (328.53971, -14.23297), (356.00654, -2.59369), (23.35453, 9.77988), (51.18176, 19.49602),
    ]),
    -9000: (0.28, [
        (68.86565, 23.15852), (101.09314, 24.52035), (129.72575, 19.93837), (159.27112, 9.86039),
        (186.57006, -2.42165), (215.84638, -14.69719), (247.17259, -22.89593), (281.86072, -24.46830),
        (315.05975, -18.49243), (343.78929, -7.94010), (11.34660, 4.60409), (38.38985, 15.58690),
    ]),
    -8000: (0.21, [
        (54.64972, 20.33247), (86.20455, 24.64977), (116.34441, 22.61153), (146.81754, 14.53534),
        (174.39236, 3.06950), (202.86092, -9.69826), (232.73093, -19.85930), (266.73924, -24.66597),
        (301.05363, -21.74932), (331.29280, -12.88634), (359.55224, -0.70277), (26.21007, 11.06694),
    ]),
    -7000: (0.16, [
        (41.88968, 16.73983), (72.49373, 23.49933), (101.60987, 24.23094), (133.03700, 18.75034),
        (161.31500, 8.71049), (189.45179, -3.91836), (217.95752, -15.44901), (250.58902, -23.25594),
        (285.28063, -23.92498), (317.15050, -17.56320), (346.65814, -6.40951), (13.39630, 5.68173),
    ]),
    -6000: (0.11, [
        (28.63744, 12.06688), (58.13629, 21.01449), (87.80311, 24.46554), (119.82368, 21.69147),
        (148.96861, 13.45591), (177.27442, 1.53693), (204.89188, -10.59427), (236.01627, -20.54749),
        (270.19528, -24.49128), (303.22835, -21.00088), (334.22738, -11.45997), (1.53934, 0.40463),
    ]),
    -5000: (0.075, [
        (16.58842, 7.16479), (45.21422, 17.68259), (73.06118, 23.37244), (105.24099, 23.65182),
        (135.34928, 17.79492), (164.28440, 7.19775), (191.45949, -4.93244), (221.10915, -16.43350),
        (254.01330, -23.47774), (287.54002, -23.41140), (320.19629, -16.32075), (348.60698, -5.32159),
    ]),
    -4000: (0.05, [
        (3.75194, 1.54191), (31.77550, 13.21952), (59.65651, 21.16853), (91.48281, 24.23323),
        (122.27669, 20.90657), (152.06854, 12.03145), (179.30613, 0.45830), (207.94579, -11.78996),
        (239.35607, -21.10803), (272.50910, -24.22191), (306.40736, -19.99346), (336.18098, -10.43389),
    ]),
    -3000: (0.03, [
        (351.68635, -3.79030), (19.58804, 8.44833), (45.67509, 17.69272), (76.70237, 23.51317),
        (107.80815, 23.09747), (138.61034, 16.54225), (166.36658, 6.10697), (194.43664, -6.27890),
        (224.31999, -17.30098), (256.33698, -23.48037), (290.86216, -22.71737), (322.20641, -15.39985),
    ]),
    -2000: (0.015, [
        (338.35990, -9.35443), (6.64285, 2.89837), (33.08091, 13.60505), (63.18303, 21.62640),
        (94.10236, 23.92053), (125.69939, 19.87846), (154.24220, 10.97704), (182.27466, -0.96428),
        (211.05429, -12.88039), (241.64996, -21.35198), (275.91557, -23.86291), (308.49223, -19.21513),
    ]),
    -1000: (0.007, [
        (325.42015, -14.09026), (354.53024, -2.42770), (19.86853, 8.52532), (49.02702, 18.43921),
        (79.30518, 23.46645), (111.37741, 22.36949), (140.91269, 15.57868), (169.39015, 4.66685),
        (197.48118, -7.54510), (226.57154, -17.78023), (259.77697, -23.49835), (293.05233, -22.13028),
    ]),
    0: (0.003, [
        (310.82475, -18.37724), (341.22498, -8.04177), (7.73927, 3.38504), (36.26332, 14.55878),
        (65.70417, 21.80973), (97.73747, 23.51112), (128.11376, 19.05593), (157.35108, 9.59506),
        (185.29455, -2.32116), (213.25821, -13.53881), (245.05299, -21.70705), (278.19387, -23.48790),
    ]),
    1000: (0.001, [
        (296.57297, -21.31250), (328.36273, -12.88424), (354.64927, -2.32282), (22.89400, 9.63877),
        (51.41892, 18.83564), (82.92751, 23.41114), (113.89300, 21.74405), (144.14986, 14.32797),
        (172.43501, 3.28062), (199.64611, -8.35082), (229.88252, -18.45424), (262.10581, -23.37206),
    ]),
    2000: (0.001, [
        (280.73585, -23.07203), (313.88462, -17.35272), (342.20134, -7.54913), (10.64110, 4.57708),
        (38.50484, 15.10514), (69.22895, 22.06610), (100.30341, 23.10106), (131.48619, 17.99283),
        (160.46901, 8.24740), (187.45496, -3.21940), (216.48242, -14.45486), (247.39726, -21.81380),
    ]),
}


@pytest.mark.parametrize('year', sorted(SUN_DE441))
def test_sun_matches_de441_through_the_year_in_every_era(year):
    tolerance, reference = SUN_DE441[year]
    jd = np.array([_date_to_jd(year, month, 1, 0.0) for month in range(1, 13)])
    ra, dec = erfa.c2s(_sun_geocentric(jd))
    ref = np.radians(reference)
    sep = np.degrees(erfa.seps(ra, dec, ref[:, 0], ref[:, 1]))
    assert sep.max() < tolerance, f'{sep.max():.4f} deg'
    # Declination sets where the Sun rises on the horizon; it stays far closer than the total.
    assert np.degrees(np.abs(dec - ref[:, 1])).max() < 0.1


# The same comparison for the Moon and planets; JPL's values are in data/de441_moon_planets.json.
# Tolerances (degrees) per era, 9998 BC to AD 2000, are again about 1.2 times the measured error.
# The analytic theories are tuned to recent millennia, so before about 5000 BC the Moon, Jupiter
# and Saturn drift by 1 to 2 degrees.
BODIES_DE441 = json.loads((Path(__file__).parent / 'data' / 'de441_moon_planets.json').read_text())
BODY_TOLERANCES = {
    'Moon': (2.2, 2.1, 1.7, 1.3, 0.7, 0.54, 0.38, 0.2, 0.098, 0.044, 0.026, 0.0092, 0.008),
    'Mercury': (0.57, 0.61, 0.46, 0.32, 0.2, 0.13, 0.078, 0.039, 0.019, 0.0093, 0.0037, 0.002, 0.002),
    'Venus': (0.62, 0.72, 0.19, 0.3, 0.086, 0.12, 0.033, 0.072, 0.011, 0.023, 0.002, 0.0041, 0.002),
    'Mars': (0.22, 0.42, 0.12, 0.14, 0.22, 0.045, 0.034, 0.047, 0.0076, 0.0071, 0.011, 0.0092, 0.002),
    'Jupiter': (1.5, 1.1, 0.28, 0.92, 0.23, 0.12, 0.3, 0.082, 0.089, 0.032, 0.014, 0.002, 0.018),
    'Saturn': (1.8, 0.13, 0.92, 1.1, 0.96, 0.84, 0.47, 0.23, 0.091, 0.13, 0.066, 0.0089, 0.023),
}


@pytest.mark.parametrize('era', range(len(BODIES_DE441['years'])))
@pytest.mark.parametrize('body', sorted(BODY_TOLERANCES))
def test_moon_and_planets_match_de441_in_every_era(body, era):
    year = BODIES_DE441['years'][era]
    jd = np.array([_date_to_jd(year, month, 1, 0.0) for month in range(1, 13)])
    ra, dec = erfa.c2s(geocentric_positions(jd)[body])
    ref = np.radians(BODIES_DE441['bodies'][body][era])
    sep = np.degrees(erfa.seps(ra, dec, ref[:, 0], ref[:, 1]))
    assert sep.max() < BODY_TOLERANCES[body][era], f'{year}: {sep.max():.4f} deg'


def test_bodies_available_at_10500_bc():
    bodies = body_altaz(_date_to_jd(-10499, 3, 20, 20.0), GIZA_LAT, GIZA_LON)
    assert set(bodies) == {'Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn'}
    assert all(np.isfinite(v).all() for v in bodies.values())
