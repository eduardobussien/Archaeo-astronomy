import numpy as np
import pytest

from alignment import _date_to_jd, calculate_alignments, find_heliacal_rising
from alignment_finder import _sun_centre_altitude, find_axis_crossings, sun_events_in_year
from data import MONUMENTS
from horizon import MAX_POINTS, HorizonProfile, as_profile
from solar_system import sun_altaz

STONEHENGE = MONUMENTS['Stonehenge']
GIZA = MONUMENTS['Great Pyramid of Giza']


def test_interpolates_linearly_between_points():
    p = HorizonProfile([0, 90, 180, 270], [1.0, 3.0, 2.0, 0.0])
    assert p(45) == pytest.approx(2.0)
    assert p(135) == pytest.approx(2.5)
    assert np.allclose(p([90, 180]), [3.0, 2.0])


def test_wraps_through_north():
    p = HorizonProfile([10, 350], [4.0, 2.0])
    assert p(0) == pytest.approx(3.0)
    assert p(360) == pytest.approx(3.0)
    assert p(355) == pytest.approx(2.5)
    assert p(-5) == pytest.approx(2.5)              # the same direction as 355
    assert p(180) == pytest.approx(3.0)             # halfway round the long way


def test_a_single_point_is_a_level_horizon():
    p = HorizonProfile([123], [1.5])
    assert np.allclose(p([0, 90, 359.9]), 1.5)


def test_reads_survey_and_stellarium_style_text():
    text = """azimuth,altitude,distance_km
    # surveyed 2026-05-01
    0, 1.0, 12
    90\t3.0\t4
    180 2.0      # behind the ridge
    270;0.5
    """
    p = HorizonProfile.from_text(text)
    assert list(p.azimuths) == [0, 90, 180, 270]
    assert list(p.altitudes) == [1.0, 3.0, 2.0, 0.5]


def test_compact_form_round_trips():
    p = HorizonProfile([270, 0, 90.5], [0.5, 1.25, 3.0])
    assert p.compact() == '0:1.25,90.5:3,270:0.5'
    q = HorizonProfile.from_compact(p.compact())
    assert np.array_equal(q.azimuths, p.azimuths) and np.array_equal(q.altitudes, p.altitudes)


@pytest.mark.parametrize('text, message', [
    ('', 'no points'),
    ('# only a comment', 'no points'),
    ('0 1\n90', 'needs an azimuth and an altitude'),
    ('0 1\nnorth 2', 'Line 2'),
    ('0 1\n90 50', 'between -2 and 45'),
    ('0 1\n400 2', 'between 0 and 360'),
    ('0 1\n360 2', 'only once'),
    ('0 1\nnan 2', 'must be numbers'),
])
def test_rejects_bad_profiles_with_a_clear_reason(text, message):
    with pytest.raises(ValueError, match=message):
        HorizonProfile.from_text(text)


def test_rejects_too_many_points_and_bad_compact_text():
    with pytest.raises(ValueError, match=f'1 to {MAX_POINTS}'):
        HorizonProfile(np.linspace(0, 359.9, MAX_POINTS + 1), np.zeros(MAX_POINTS + 1))
    with pytest.raises(ValueError, match='azimuth:altitude'):
        HorizonProfile.from_compact('0:1,90')


def test_a_number_is_a_level_horizon_and_the_outline_closes():
    level = as_profile(1.5)
    assert level.is_level and level(77) == 1.5
    assert as_profile(level) is level
    az, alt = HorizonProfile([0, 90.5, 180], [1, 3, 2]).outline()
    assert az[0] == 0 and az[-1] == 360 and alt[0] == alt[-1]
    assert 90.5 in az and alt[list(az).index(90.5)] == 3     # every surveyed point is on it


# ── The sky and the searches follow the skyline ──

def test_visible_means_above_the_skyline():
    profile = HorizonProfile([0, 90, 180, 270], [0, 30, 0, 10])
    results = calculate_alignments(STONEHENGE['lat'], STONEHENGE['lon'], -2499, 6, 21, 4.0,
                                   horizon=profile)
    for body in [*results['stars'].values(), *results['planets'].values()]:
        assert body['visible'] == (body['altitude'] > profile(body['azimuth']))
    assert any(0 < s['altitude'] < profile(s['azimuth']) for s in results['stars'].values())


def _transitions(lat, lon, midnight, profile, rising):
    """Every moment, to the second, the Sun's upper limb crosses the skyline in a half day."""
    jd = midnight + (np.arange(0, 43200) + (0 if rising else 43200)) / 86400.0
    alt, az = sun_altaz(jd, lat, lon)
    clear = alt >= _sun_centre_altitude(profile(az), 'upper', 0.0, True)
    return jd[1:][clear[1:] != clear[:-1]]


def test_sunrise_and_sunset_over_a_jagged_skyline_match_a_second_by_second_scan():
    # Teeth 4 degrees high every 2 degrees: the Sun slips in and out of view as it
    # rises or sets, so sunrise must be the first gleam and sunset the last.
    sawtooth = HorizonProfile(np.arange(360), np.where(np.arange(360) % 2, 4.0, 0.0))
    lat, lon = STONEHENGE['lat'], STONEHENGE['lon']
    first = _date_to_jd(2025, 1, 1, -lon / 15.0)
    several_crossings = False
    for event in ('rise', 'set'):
        jd, _ = sun_events_in_year(lat, lon, 2025, event, sawtooth, elevation_m=0.0)
        for day in (0, 80, 172, 265):
            moments = _transitions(lat, lon, first + day, sawtooth, event == 'rise')
            expected = moments[0] if event == 'rise' else moments[-1]
            assert abs(jd[day] - expected) * 86400.0 < 2.0
            several_crossings |= len(moments) > 1
    assert several_crossings


def test_sunrise_depends_only_on_the_skyline_where_the_sun_rises():
    lat, lon = STONEHENGE['lat'], STONEHENGE['lon']
    western_peak = HorizonProfile([0, 250, 270, 290], [1.5, 1.5, 20.0, 1.5])
    with_peak, _ = sun_events_in_year(lat, lon, 2025, 'rise', western_peak)
    level, _ = sun_events_in_year(lat, lon, 2025, 'rise', 1.5)
    assert np.nanmax(np.abs(with_peak - level)) * 86400.0 < 1.0


def test_a_star_meets_the_axis_where_the_skyline_crosses_it():
    lat, lon = 42.5, 12.3
    profile = HorizonProfile([0, 100, 140], [0.0, 1.0, 3.0])   # 2.0 at the 120 degree axis
    with_profile = find_axis_crossings(lat, lon, 'Sirius', 120.0, horizon=profile)
    level = find_axis_crossings(lat, lon, 'Sirius', 120.0, horizon=2.0)
    assert with_profile['crossings'] == level['crossings'] and level['crossings']


def test_a_ridge_where_sirius_rises_delays_its_heliacal_rising():
    # Over a level horizon Sirius is first seen at 4.5 deg (it needs 3.8 to shine
    # through the haze) and stands about 1 deg higher at each dawn, so a 10 deg
    # ridge holds it back about a week.
    args = GIZA['lat'], GIZA['lon'], -2780, 'Sirius'
    level = find_heliacal_rising(*args, elevation_m=GIZA['elevation_m'])
    ridge = HorizonProfile([0, 95, 100, 130, 135], [0, 0, 10, 10, 0])
    behind = find_heliacal_rising(*args, elevation_m=GIZA['elevation_m'], horizon=ridge)
    assert level['found'] and behind['found']
    assert 100 < behind['star_azimuth'] < 130 and behind['horizon_altitude'] == 10.0
    assert behind['star_altitude'] > 10.0
    days_later = (_date_to_jd(behind['year'], behind['month'], behind['day'], 0)
                  - _date_to_jd(level['year'], level['month'], level['day'], 0))
    assert 5 <= days_later <= 9
