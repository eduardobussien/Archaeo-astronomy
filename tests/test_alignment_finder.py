import json
from pathlib import Path

import erfa
import numpy as np
import pytest

from alignment import _date_to_jd, _observed
from alignment_finder import (_sun_centre_altitude, event_azimuths, find_axis_crossings,
                              find_sun_alignment, horizon_azimuth, sun_events_in_year)
from data import MONUMENTS, STARS
from precession import J2000_JD, mean_obliquity, star_altaz
from refraction import apparent_altitude
from solar_system import sun_altaz

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


# ── The Sun ──

def test_obliquity_matches_iau_2006_at_j2000():
    assert float(mean_obliquity(J2000_JD)) * 3600 == pytest.approx(84381.406, abs=0.01)
    century = np.diff(mean_obliquity(J2000_JD + np.array([-18262.5, 18262.5])))[0] * 3600
    assert century == pytest.approx(-46.84, abs=0.1)


@pytest.mark.parametrize('site, year', [('Stonehenge', -2499), ('Stonehenge', -8999),
                                        ('Newgrange', -3199), ('Machu Picchu', 1450)])
def test_solstice_formula_agrees_with_the_years_sunrises(site, year):
    # Two independent routes: the obliquity in a closed formula, and every
    # sunrise of the year found from the solar ephemeris.
    m = MONUMENTS[site]
    _, az = sun_events_in_year(m['lat'], m['lon'], year, 'rise', elevation_m=m['elevation_m'])
    eps = float(mean_obliquity(_date_to_jd(year, 7, 2, 0.0)))
    h0 = _sun_centre_altitude(0.0, 'upper', m['elevation_m'], True)
    assert np.nanmin(az) == pytest.approx(float(horizon_azimuth(eps, m['lat'], h0)), abs=0.01)
    assert np.nanmax(az) == pytest.approx(float(horizon_azimuth(-eps, m['lat'], h0)), abs=0.01)


def test_reported_sunrises_put_the_first_gleam_on_the_axis():
    m = STONEHENGE
    result = find_sun_alignment(m['lat'], m['lon'], m['orientation_az'], 'rise', -2499,
                                elevation_m=m['elevation_m'])
    assert len(result['dates']) == 2
    for d in result['dates']:
        jd = _date_to_jd(-2499, d['month'], d['day'], d['hour'] - m['lon'] / 15.0)
        alt, az = sun_altaz(jd, m['lat'], m['lon'])
        upper_limb = float(apparent_altitude(alt + 0.2666, m['elevation_m']))  # refract the limb itself
        assert abs(upper_limb) < 0.01
        assert abs(float(az) - m['orientation_az']) < 0.5     # nearest whole day


# Sunrise and sunset in UT from the US Naval Observatory for 2025, from 59 N to 46 S; the times
# are rounded to the minute. USNO takes 34' of refraction at the horizon where this app's
# Saemundsson formula gives 34.4', worth a few seconds.
USNO = json.loads((Path(__file__).parent / 'data' / 'usno_sun_2025.json').read_text())


@pytest.mark.parametrize('event', ['rise', 'set'])
@pytest.mark.parametrize('site', sorted(USNO['sites']))
def test_sunrise_and_sunset_times_match_the_us_naval_observatory(site, event):
    s = USNO['sites'][site]
    jd, _ = sun_events_in_year(s['lat'], s['lon'], USNO['year'], event)
    jd = jd[np.isfinite(jd)]
    _, month, day, fraction = erfa.jd2cal(jd, 0.0)
    for date, times in s['rise_set'].items():
        m, d = map(int, date.split('-'))
        hh, mm = map(int, times[0 if event == 'rise' else 1].split(':'))
        (k,) = np.flatnonzero((month == m) & (day == d))
        # 30 s of rounding, plus a margin. Refracting the Sun's centre instead of its limb
        # put sunrise up to 74 s early here.
        assert abs(fraction[k] * 86400.0 - (hh * 3600 + mm * 60)) < 45.0, date


def test_june_sunrise_is_north_of_east_in_the_southern_hemisphere_too():
    m = MONUMENTS['Machu Picchu']
    result = find_sun_alignment(m['lat'], m['lon'], 65.0, 'rise', 1450, elevation_m=m['elevation_m'])
    assert result['june']['month'] == 6 and result['june']['azimuth'] < 90
    assert result['december']['month'] == 12 and result['december']['azimuth'] > 90


def test_finds_the_year_a_solstice_sunrise_fell_on_its_azimuth():
    m = STONEHENGE
    eps = float(mean_obliquity(_date_to_jd(-4000, 7, 2, 0.0)))
    target = float(horizon_azimuth(eps, m['lat'], _sun_centre_altitude(0.0, 'upper', m['elevation_m'], True)))
    result = find_sun_alignment(m['lat'], m['lon'], target, 'rise', -2499, elevation_m=m['elevation_m'])
    assert result['solstice']['which'] == 'june'
    assert -4000 in [c['year'] for c in result['solstice']['crossings']]
