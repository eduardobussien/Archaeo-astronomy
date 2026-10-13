"""The monument data against the publications it comes from, and the engine against them."""
import pytest

from alignment import _date_to_jd
from alignment_finder import find_sun_alignment, horizon_azimuth
from data import MONUMENTS
from precession import mean_obliquity
from solar_system import sun_altaz


def _site(name):
    m = MONUMENTS[name]
    return m['lat'], m['lon'], m['elevation_m'], m['orientation_az'], m['horizon_deg'] or 0.0


def _sun_dates(name, year, event, limb='upper'):
    lat, lon, elevation, azimuth, horizon = _site(name)
    result = find_sun_alignment(lat, lon, azimuth, event, year, horizon, elevation, limb=limb)
    return [(d['month'], d['day']) for d in result['dates']]


@pytest.mark.parametrize('name', sorted(MONUMENTS))
def test_every_axis_is_sourced_or_flagged_approximate(name):
    m = MONUMENTS[name]
    assert -90 <= m['lat'] <= 90 and -180 <= m['lon'] <= 180 and 0 <= m['orientation_az'] < 360
    assert m['horizon_deg'] is None or -2 <= m['horizon_deg'] <= 10
    assert m['source'] or 'approximate' in m['note']


def test_teotihuacan_sunsets_match_sprajc_2000():
    # Šprajc computed the dates for the Sun's centre on the pyramid's axis and horizon.
    assert _sun_dates('Pyramid of the Sun (Teotihuacan)', 150, 'set', 'centre') == [(4, 30), (8, 13)]


def test_el_castillo_sunsets_match_sprajc_and_sanchez_nava():
    assert _sun_dates('El Castillo (Chichen Itza)', 1000, 'set', 'centre') == [(5, 20), (7, 24)]


def test_baalbek_sunrises_match_the_bekaa_valley_study():
    (m1, d1), (m2, d2) = _sun_dates('Baalbek (Temple of Jupiter)', 50, 'rise')
    assert (m1, m2) == (5, 8) or (m1, m2) == (4, 8)
    assert abs(_date_to_jd(50, m1, d1, 0) - _date_to_jd(50, 4, 30, 0)) <= 2
    assert abs(_date_to_jd(50, m2, d2, 0) - _date_to_jd(50, 8, 12, 0)) <= 2


def test_mnajdra_axis_marks_the_equinox_sunrise():
    # Micallef (2000) gives declination 0 for the lower temple axis.
    (m1, d1), (m2, d2) = _sun_dates('Mnajdra Temple (Malta)', -2874, 'rise')
    assert m1 == 3 and 19 <= d1 <= 23 and m2 == 9 and 20 <= d2 <= 25


def test_newgrange_1967_light_falls_inside_patricks_window():
    # The O'Kellys timed the beam on 21 December 1967: in at 08:58 GMT, cut off at 09:15.
    # Patrick (1974) surveyed the passage: sunlight reaches the chamber from 133.7 to 138.4 deg.
    lat, lon, *_ = _site('Newgrange')
    for hour in (8 + 58 / 60, 9 + 15 / 60):
        _, az = sun_altaz(_date_to_jd(1967, 12, 21, hour), lat, lon)
        assert 133.7 <= float(az) <= 138.4


def test_stonehenge_midsummer_sunrise_on_a_level_horizon_in_2500_bc():
    # sarsen.org (2026): 49.6 deg for the Sun's centre on a level horizon without
    # refraction, with the Vondrák obliquity of 23.973 deg.
    lat = MONUMENTS['Stonehenge']['lat']
    eps = float(mean_obliquity(_date_to_jd(-2499, 6, 21, 0.0)))
    assert eps == pytest.approx(23.973, abs=0.001)
    assert float(horizon_azimuth(eps, lat, 0.0)) == pytest.approx(49.6, abs=0.05)
