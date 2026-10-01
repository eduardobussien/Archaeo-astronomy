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


def test_heliacal_finds_sirius():
    body = _get('/api/heliacal', site='Great Pyramid of Giza', year='-2780', star='Sirius').get_json()
    assert body['found'] and body['month'] == 6
