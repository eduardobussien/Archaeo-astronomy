import csv

import numpy as np
import pytest

from app import app

CLIENT = app.test_client()


def _get(path, **params):
    return CLIENT.get(path, query_string=params)


def test_stars_default_request_succeeds():
    r = _get('/api/stars')
    assert r.status_code == 200
    body = r.get_json()
    assert body['meta']['year'] == -2499
    assert len(body['stars']) == 60
    assert set(body['planets']) == {'Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn'}


def test_stars_response_includes_the_ecliptic():
    params = {'site': 'Stonehenge', 'year': '-2499', 'month': '6', 'day': '21', 'hour': '4'}
    stars = _get('/api/stars', **params).get_json()
    ecliptic = _get('/api/ecliptic', **params).get_json()
    assert len(stars['ecliptic']) == 73
    assert stars['ecliptic'] == ecliptic['points']


def test_stars_response_reports_delta_t_and_its_uncertainty():
    ancient = _get('/api/stars', year='-2499').get_json()['meta']['delta_t']
    assert not ancient['measured']
    assert ancient['seconds'] == pytest.approx(59700, rel=0.01)       # about 16.6 hours
    low, high = ancient['moon_uncertainty_deg']
    assert low < high and high == pytest.approx(0.93, abs=0.01)     # the README table, 2500 BC
    modern = _get('/api/stars', year='2000').get_json()['meta']['delta_t']
    assert modern['measured'] and modern['uncertainty_s'][1] < 5.0   # 0.8 t^2 gives 2.6 s


def test_site_lookup_is_exact_and_case_insensitive():
    r = _get('/api/stars', site='stonehenge')
    assert r.status_code == 200
    assert r.get_json()['monument']['name'] == 'Stonehenge'


@pytest.mark.parametrize('path', ['/api/stars', '/api/ecliptic', '/api/heliacal'])
def test_unknown_site_is_404_on_every_endpoint(path):
    r = _get(path, site='a')
    assert r.status_code == 404
    assert 'error' in r.get_json()


@pytest.mark.parametrize('params', [
    {'lat': 'abc'},
    {'lat': '120'},
    {'lat': 'nan'},
    {'lon': '-181'},
    {'lat': ''},
    {'year': 'xyz'},
    {'year': '2.5'},
    {'year': '-30000'},
    {'month': '13'},
    {'month': '2', 'day': '31'},
    {'year': '2023', 'month': '2', 'day': '29'},
    {'hour': '25'},
])
@pytest.mark.parametrize('path', ['/api/stars', '/api/ecliptic'])
def test_invalid_parameters_return_json_400(path, params):
    r = _get(path, **params)
    assert r.status_code == 400
    assert r.is_json and r.get_json()['error']


@pytest.mark.parametrize('year', ['2024', '0', '-400'])
def test_february_29_accepted_in_leap_years(year):
    assert _get('/api/stars', year=year, month='2', day='29').status_code == 200


@pytest.mark.parametrize('params', [
    {'star': 'Nibiru'},
    {'extinction': '5'},
    {'extinction': 'clear'},
    {'year': '-30000'},
])
def test_heliacal_rejects_invalid_parameters(params):
    r = _get('/api/heliacal', **params)
    assert r.status_code == 400
    assert r.get_json()['error']


def test_refraction_lifts_objects_near_the_horizon():
    params = {'site': 'Great Pyramid of Giza', 'year': '-2499', 'month': '3', 'day': '20', 'hour': '22'}
    apparent = _get('/api/stars', **params).get_json()
    geometric = _get('/api/stars', refraction='0', **params).get_json()
    assert apparent['meta']['refraction'] and not geometric['meta']['refraction']
    for name, star in apparent['stars'].items():
        lift = star['altitude'] - geometric['stars'][name]['altitude']
        assert 0.0 <= lift <= 0.66
        assert star['azimuth'] == geometric['stars'][name]['azimuth']


def test_high_sites_refract_less():
    params = {'year': '-2499', 'month': '3', 'day': '20', 'hour': '22', 'lat': '-16.5544', 'lon': '-68.6742'}
    sea = _get('/api/stars', elevation='0', **params).get_json()
    high = _get('/api/stars', elevation='3850', **params).get_json()
    flat = _get('/api/stars', refraction='0', **params).get_json()
    name = min(flat['stars'], key=lambda n: abs(flat['stars'][n]['altitude'] - 5))
    lift_sea = sea['stars'][name]['altitude'] - flat['stars'][name]['altitude']
    lift_high = high['stars'][name]['altitude'] - flat['stars'][name]['altitude']
    assert abs(lift_high / lift_sea - 0.620) < 0.005


@pytest.mark.parametrize('params', [{'refraction': 'yes'}, {'elevation': '10000'}])
def test_refraction_parameters_are_validated(params):
    assert _get('/api/stars', **params).status_code == 400


def test_calendars_describe_the_same_moment():
    common = {'site': 'Great Pyramid of Giza', 'year': '-2780', 'hour': '4'}
    julian = _get('/api/stars', month='7', day='19', calendar='julian', **common).get_json()
    gregorian = _get('/api/stars', month='6', day='26', **common).get_json()
    assert julian['meta']['jd'] == gregorian['meta']['jd']
    assert julian['meta']['dates']['gregorian'] == {'year': -2780, 'month': 6, 'day': 26}
    assert gregorian['meta']['dates']['julian'] == {'year': -2780, 'month': 7, 'day': 19}


def test_calendar_is_validated_and_sets_leap_years():
    assert _get('/api/stars', calendar='mayan').status_code == 400
    assert _get('/api/stars', year='1900', month='2', day='29', calendar='julian').status_code == 200
    assert _get('/api/stars', year='1900', month='2', day='29').status_code == 400


def test_heliacal_date_follows_the_requested_calendar():
    params = {'site': 'Great Pyramid of Giza', 'year': '-2780', 'star': 'Sirius'}
    gregorian = _get('/api/heliacal', **params).get_json()
    julian = _get('/api/heliacal', calendar='julian', **params).get_json()
    assert julian['calendar'] == 'julian' and julian['month'] == 7
    assert (julian['day'] - gregorian['day']) % 30 == 23


def test_alignment_uses_the_monument_axis_by_default():
    r = _get('/api/alignment', site='Göbekli Tepe', star='Sirius')
    assert r.status_code == 200
    body = r.get_json()
    assert body['azimuth'] == 165.0 and body['found']
    night = body['crossings'][0]['night']
    assert {'month', 'day', 'hour', 'azimuth'} <= set(night)


def test_sun_alignment_lists_dates_and_solstice_range():
    body = _get('/api/alignment', site='Stonehenge', star='Sun', year='-2499').get_json()
    assert body['star'] == 'Sun' and body['limb'] == 'upper'
    assert body['azimuth'] == 49.95 and body['horizon'] == 0.6      # the surveyed axis and horizon
    assert body['june']['azimuth'] < body['azimuth'] < body['december']['azimuth']
    assert body['solstice']['which'] == 'june'


@pytest.mark.parametrize('params', [
    {'site': 'Stonehenge', 'star': 'Sun', 'limb': 'top'},
    {'site': 'Stonehenge', 'azimuth': '200', 'event': 'rise'},
    {'site': 'Stonehenge', 'azimuth': '60', 'event': 'set'},
    {'site': 'Great Pyramid of Giza'},
    {'lat': '30', 'lon': '31'},
    {'site': 'Stonehenge', 'event': 'culminate'},
    {'site': 'Stonehenge', 'star': 'Nibiru'},
    {'site': 'Stonehenge', 'horizon': '50'},
    {'site': 'Stonehenge', 'horizon': '0:1,90'},
    {'site': 'Stonehenge', 'from': '-2000', 'to': '-3000'},
])
def test_alignment_rejects_invalid_requests(params):
    r = _get('/api/alignment', **params)
    assert r.status_code == 400
    assert r.get_json()['error']


def _csv_table(response):
    """(header lines as a dict, rows as dicts) from a CSV download."""
    lines = response.get_data(as_text=True).splitlines()
    header = dict(l[2:].split(': ', 1) for l in lines if l.startswith('# '))
    return header, list(csv.DictReader(l for l in lines if not l.startswith('#')))


def test_stars_csv_matches_the_json_and_says_how_it_was_made():
    params = {'site': 'Stonehenge', 'year': '-2499', 'month': '6', 'day': '21', 'hour': '4'}
    r = _get('/api/stars', format='csv', **params)
    assert r.status_code == 200 and r.mimetype == 'text/csv'
    assert 'attachment' in r.headers['Content-Disposition']
    header, rows = _csv_table(r)
    assert header['site'] == 'Stonehenge' and header['local date (gregorian)'] == '2500 BC, June 21'
    assert 'format=csv' in header['request']
    body = _get('/api/stars', **params).get_json()
    assert len(rows) == len(body['stars']) + len(body['planets'])
    sirius = next(row for row in rows if row['name'] == 'Sirius')
    assert float(sirius['altitude_deg']) == pytest.approx(body['stars']['Sirius']['altitude'], abs=1e-4)


def test_alignment_csv_lists_the_same_crossings_as_the_json():
    params = {'lat': '42.5', 'lon': '12.3', 'elevation': '300', 'star': 'Sirius',
              'azimuth': '120', 'horizon': '1.5'}
    body = _get('/api/alignment', **params).get_json()
    header, rows = _csv_table(_get('/api/alignment', format='csv', **params))
    assert header['site'] == 'custom coordinates' and header['horizon altitude_deg'] == '1.5'
    crossings = [int(row['year']) for row in rows if row['section'] == 'crossing']
    assert crossings == [c['year'] for c in body['crossings']]
    assert sum(row['section'] == 'curve' for row in rows) > 100


def test_sun_alignment_csv_lists_the_dates_and_solstices():
    params = {'site': 'Pyramid of the Sun (Teotihuacan)', 'star': 'Sun', 'event': 'set',
              'year': '150', 'limb': 'centre'}
    header, rows = _csv_table(_get('/api/alignment', format='csv', **params))
    assert header['part of the Sun'] == 'centre'
    dates = [(int(row['month']), int(row['day'])) for row in rows if row['section'] == 'date']
    assert dates == [(4, 30), (8, 13)]                     # Šprajc (2000)
    assert {row['section'] for row in rows} >= {'june_solstice', 'december_solstice', 'solstice_curve'}


@pytest.mark.parametrize('path', ['/api/stars', '/api/alignment'])
def test_unknown_format_is_rejected(path):
    r = _get(path, site='Stonehenge', format='xlsx')
    assert r.status_code == 400 and r.get_json()['error']


PROFILE = '0:0,90:8,180:0,270:2'


def test_horizon_endpoint_reads_a_profile_file():
    r = CLIENT.post('/api/horizon', data='azimuth,altitude\n0,0\n90,8   # the ridge\n180\t0\n270 2\n')
    assert r.status_code == 200
    assert r.get_json() == {'profile': PROFILE, 'points': 4, 'lowest': 0.0, 'highest': 8.0}
    bad = CLIENT.post('/api/horizon', data='0 1\nnorth 2\n')
    assert bad.status_code == 400 and 'Line 2' in bad.get_json()['error']


def test_stars_are_visible_only_above_the_skyline_and_the_outline_is_returned():
    params = {'site': 'Stonehenge', 'year': '-2499', 'month': '6', 'day': '21', 'hour': '4'}
    assert _get('/api/stars', **params).get_json()['horizon'] is None
    body = _get('/api/stars', horizon=PROFILE, **params).get_json()
    outline = body['horizon']['outline']
    assert body['horizon']['profile'] == PROFILE and outline[0] == [0.0, 0.0] and outline[-1][0] == 360
    assert [90.0, 8.0] in outline
    for s in body['stars'].values():
        assert s['visible'] == (s['altitude'] > np.interp(s['azimuth'], [0, 90, 180, 270, 360],
                                                           [0, 8, 0, 2, 0]))
    header, _ = _csv_table(_get('/api/stars', horizon=PROFILE, format='csv', **params))
    assert header['horizon'].endswith(PROFILE)


def test_alignment_with_a_profile_uses_the_skyline_at_the_axis():
    params = {'site': 'Stonehenge', 'star': 'Sun', 'year': '-2499', 'horizon': '0:0,45:2,90:0'}
    body = _get('/api/alignment', **params).get_json()
    assert body['horizon'] == pytest.approx(2 - 2 * (49.95 - 45) / 45)     # 2 at 45, 0 at 90
    assert body['horizon_profile'] == '0:0,45:2,90:0'
    level = _get('/api/alignment', site='Stonehenge', star='Sun', year='-2499').get_json()
    assert level['horizon_profile'] is None and level['horizon'] == 0.6


def test_heliacal_finds_sirius():
    body = _get('/api/heliacal', site='Great Pyramid of Giza', year='-2780', star='Sirius').get_json()
    assert body['found'] and body['month'] == 6
