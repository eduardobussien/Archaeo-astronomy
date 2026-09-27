import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from flask import Flask, jsonify, request, render_template
from alignment import (calculate_alignments, calculate_ecliptic, days_in_month,
                       find_heliacal_rising, format_year)
from data import MONUMENTS, STARS
from precession import EPOCH_RANGE

app = Flask(__name__)

# Keeps every requested date, and the heliacal scan's previous-year start, inside the model.
YEAR_RANGE = (math.ceil(EPOCH_RANGE[0]) + 2, math.floor(EPOCH_RANGE[1]) - 2)
_SITE_NAMES = {name.lower(): name for name in MONUMENTS}


class ApiError(Exception):
    """A request the API cannot serve; reported to the client as JSON."""

    def __init__(self, message, status=400):
        super().__init__(message)
        self.status = status


@app.errorhandler(ApiError)
def handle_api_error(err):
    return jsonify({'error': str(err)}), err.status


@app.errorhandler(500)
def handle_server_error(err):
    if request.path.startswith('/api/'):
        return jsonify({'error': 'Internal server error'}), 500
    return err


def _number(name, default, cast, low, high):
    raw = request.args.get(name)
    if raw is None:
        return default
    kind = 'a whole number' if cast is int else 'a number'
    try:
        value = cast(raw)
    except ValueError:
        raise ApiError(f"'{name}' must be {kind}") from None
    if not (math.isfinite(value) and low <= value <= high):
        raise ApiError(f"'{name}' must be between {low} and {high}")
    return value


def _location():
    """(lat, lon, monument name or None) from ?site= or ?lat=&lon=."""
    site = request.args.get('site')
    if site:
        name = _SITE_NAMES.get(site.strip().lower())
        if name is None:
            raise ApiError("Unknown site; see /api/sites for the list", 404)
        return MONUMENTS[name]['lat'], MONUMENTS[name]['lon'], name
    return (_number('lat', 29.9792, float, -90.0, 90.0),
            _number('lon', 31.1342, float, -180.0, 180.0),
            None)


def _date_time():
    """(year, month, day, hour): astronomical year, proleptic Gregorian date, local hour."""
    year = _number('year', -2499, int, *YEAR_RANGE)
    month = _number('month', 3, int, 1, 12)
    day = _number('day', 20, int, 1, days_in_month(year, month))
    hour = _number('hour', 22.0, float, 0.0, 24.0)
    return year, month, day, hour


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/stars')
def stars():
    """
    Return star, Sun, Moon and planet positions for a location and historical date.

    Query parameters:
        site  - monument name from /api/sites (overrides lat/lon)
        lat   - latitude in degrees, -90 to 90 (default: Giza)
        lon   - longitude in degrees, -180 to 180 (default: Giza)
        year  - astronomical year: 0 = 1 BC, -2499 = 2500 BC (default: -2499)
        month - 1-12 (default: 3)
        day   - 1 to the length of the month (default: 20)
        hour  - local mean solar time, 0-24 (default: 22.0)

    Invalid parameters return HTTP 400 (404 for an unknown site) with {"error": message}.
    """
    lat, lon, site = _location()
    year, month, day, hour = _date_time()
    results = calculate_alignments(lat, lon, year, month, day, hour)

    monument_info = None
    if site:
        monument_info = {
            'name':           site,
            'orientation_az': MONUMENTS[site]['orientation_az'],
            'note':           MONUMENTS[site]['note'],
        }

    return jsonify({
        'meta': {
            'lat':     lat,
            'lon':     lon,
            'year':    year,
            'era':     'BC' if year <= 0 else 'AD',
            'month':   month,
            'day':     day,
            'hour':    hour,
            'hour_ut': round(hour - lon / 15.0, 4),
            'jd':      float(results['jd']),
            'lst':     float(results['lst']),
            'method':  results['method'],
        },
        'monument': monument_info,
        'stars': {
            name: {
                'altitude':      float(d['altitude']),
                'azimuth':       float(d['azimuth']),
                'visible':       bool(d['visible']),
                'magnitude':     STARS[name]['mag'],
                'constellation': STARS[name]['constellation'],
            }
            for name, d in results['stars'].items()
        },
        'planets': {
            name: {
                'altitude': float(d['altitude']),
                'azimuth':  float(d['azimuth']),
                'visible':  bool(d['visible']),
            }
            for name, d in results['planets'].items()
        },
    })


@app.route('/api/heliacal')
def heliacal():
    """
    Find the heliacal rising of a star in a given year and location.

    Query parameters: site or lat/lon (as /api/stars), year (astronomical),
    star (catalog name, default Sirius), arc_vision (Sun altitude at the
    moment of observation, -20 to 0 degrees, default -10).
    """
    lat, lon, _ = _location()
    year = _number('year', -2780, int, *YEAR_RANGE)
    arc_vision = _number('arc_vision', -10.0, float, -20.0, 0.0)
    star = request.args.get('star', 'Sirius')
    if star not in STARS:
        raise ApiError("Unknown star; it must be one of the catalog names")

    result = find_heliacal_rising(lat, lon, year, star, arc_vision)
    if result is None:
        return jsonify({
            'found': False,
            'message': f"No heliacal rising of {star} found at this location in {format_year(year)}.",
        })

    return jsonify({'found': True, **result})


@app.route('/api/ecliptic')
def ecliptic():
    """
    Return the ecliptic great circle (73 points) in local alt-az.

    Accepts the same query parameters as /api/stars.
    """
    lat, lon, _ = _location()
    year, month, day, hour = _date_time()
    return jsonify({'points': calculate_ecliptic(lat, lon, year, month, day, hour)})


@app.route('/api/sites')
def sites():
    """Return the list of all known monument sites."""
    return jsonify({
        name: {
            'lat':            data['lat'],
            'lon':            data['lon'],
            'orientation_az': data['orientation_az'],
            'note':           data['note'],
        }
        for name, data in MONUMENTS.items()
    })


if __name__ == '__main__':
    # The Werkzeug debugger can run arbitrary code, so it is opt-in: FLASK_DEBUG=1.
    app.run(debug=os.environ.get('FLASK_DEBUG') == '1')
