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
    {'arc_vision': '5'},
    {'year': '-30000'},
])
def test_heliacal_rejects_invalid_parameters(params):
    r = _get('/api/heliacal', **params)
    assert r.status_code == 400
    assert r.get_json()['error']


def test_heliacal_finds_sirius():
    body = _get('/api/heliacal', site='Great Pyramid of Giza', year='-2780', star='Sirius').get_json()
    assert body['found'] and body['month'] == 6
