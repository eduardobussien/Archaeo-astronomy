import numpy as np
import pytest

from alignment import _date_to_jd, _observed
from alignment_finder import event_azimuths, find_axis_crossings
from data import MONUMENTS, STARS
from precession import star_altaz

STONEHENGE = MONUMENTS['Stonehenge']
GIZA = MONUMENTS['Great Pyramid of Giza']


def _azimuth(site, star, year, event='rise', **kwargs):
    return float(event_azimuths(site['lat'], star, [year], event,
                                elevation_m=site['elevation_m'], **kwargs)[0])


@pytest.mark.parametrize('star, event, year', [
    ('Sirius', 'rise', -4321), ('Aldebaran', 'set', -2000), ('Altair', 'rise', 1500),
])
def test_finds_the_year_whose_azimuth_it_was_given(star, event, year):
    target = _azimuth(STONEHENGE, star, year, event)
    result = find_axis_crossings(STONEHENGE['lat'], STONEHENGE['lon'], star, target, event,
                                 elevation_m=STONEHENGE['elevation_m'])
    assert year in [c['year'] for c in result['crossings']]


@pytest.mark.parametrize('star, event, year', [('Sirius', 'rise', -4321), ('Aldebaran', 'set', -2000)])
def test_reported_night_puts_the_star_on_the_horizon_at_the_axis(star, event, year):
    target = _azimuth(STONEHENGE, star, year, event)
    result = find_axis_crossings(STONEHENGE['lat'], STONEHENGE['lon'], star, target, event,
                                 elevation_m=STONEHENGE['elevation_m'])
    night = next(c['night'] for c in result['crossings'] if c['year'] == year)
    s = STARS[star]
    jd = _date_to_jd(year, night['month'], night['day'], night['hour'] - STONEHENGE['lon'] / 15.0)
    alt, az = star_altaz(s['ra'], s['dec'], s['pm_ra'], s['pm_dec'], s['dist'], jd,
                         STONEHENGE['lat'], STONEHENGE['lon'], s['rv'])
    assert abs(float(_observed(alt, STONEHENGE['elevation_m'], True))) < 0.01
    assert abs(float(az) - target) < 0.02
    assert night['hour'] < 0.5 or night['hour'] > 23.5


def test_setting_mirrors_rising():
    rise = _azimuth(GIZA, 'Sirius', -2780, 'rise')
    assert _azimuth(GIZA, 'Sirius', -2780, 'set') == pytest.approx(360.0 - rise)


def test_hills_move_a_southern_star_further_south():
    assert _azimuth(GIZA, 'Sirius', -2780, horizon=2.0) > _azimuth(GIZA, 'Sirius', -2780, horizon=0.0)


def test_refraction_shows_a_rising_star_sooner_and_further_east():
    assert _azimuth(GIZA, 'Sirius', -2780) < _azimuth(GIZA, 'Sirius', -2780, refract=False)


def test_a_circumpolar_star_never_crosses_the_horizon():
    result = find_axis_crossings(STONEHENGE['lat'], STONEHENGE['lon'], 'Polaris', 10.0, 'rise')
    assert result == {'found': False, 'reason': 'never_crosses_horizon'}


def test_reports_closest_approach_when_the_axis_is_never_reached():
    result = find_axis_crossings(GIZA['lat'], GIZA['lon'], 'Sirius', 60.0, 'rise')
    assert not result['found']
    assert result['closest']['offset'] > 0 and np.isfinite(result['closest']['azimuth'])
