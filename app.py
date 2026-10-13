import datetime
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from flask import Flask, jsonify, request, render_template
from alignment import (CALENDARS, calculate_alignments, calculate_ecliptic, days_in_month,
                       find_heliacal_rising, format_year, local_calendar_dates)
from alignment_finder import EVENTS, SUN_LIMBS, find_axis_crossings, find_sun_alignment
from data import MONUMENTS, STARS
from precession import EPOCH_RANGE
from visibility import DEFAULT_EXTINCTION

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


def _flag(name, default):
    raw = request.args.get(name)
    if raw is None:
        return default
    if raw.lower() in ('1', 'true'):
        return True
    if raw.lower() in ('0', 'false'):
        return False
    raise ApiError(f"'{name}' must be 0 or 1")


def _location():
    """(lat, lon, elevation in metres, monument name or None) from ?site= or ?lat=&lon=&elevation=."""
    site = request.args.get('site')
    if site:
        name = _SITE_NAMES.get(site.strip().lower())
        if name is None:
            raise ApiError("Unknown site; see /api/sites for the list", 404)
        m = MONUMENTS[name]
        return m['lat'], m['lon'], m['elevation_m'], name
    return (_number('lat', 29.9792, float, -90.0, 90.0),
            _number('lon', 31.1342, float, -180.0, 180.0),
            _number('elevation', 0.0, float, -500.0, 9000.0),
            None)


def _calendar():
    calendar = request.args.get('calendar', 'gregorian').lower()
    if calendar not in CALENDARS:
        raise ApiError("'calendar' must be 'gregorian' or 'julian'")
    return calendar


def _date_time():
    """(year, month, day, hour, calendar): astronomical year, date in that calendar, local hour."""
    calendar = _calendar()
    year = _number('year', -2499, int, *YEAR_RANGE)
    month = _number('month', 3, int, 1, 12)
    day = _number('day', 20, int, 1, days_in_month(year, month, calendar))
    hour = _number('hour', 22.0, float, 0.0, 24.0)
    return year, month, day, hour, calendar


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/stars')
def stars():
    """
    Return star, Sun, Moon and planet positions and the ecliptic for a location and
    historical date, so one request draws the whole sky.

    Query parameters:
        site       - monument name from /api/sites (overrides lat/lon/elevation)
        lat        - latitude in degrees, -90 to 90 (default: Giza)
        lon        - longitude in degrees, -180 to 180 (default: Giza)
        elevation  - metres above sea level, -500 to 9000 (default: 0)
        year       - astronomical year: 0 = 1 BC, -2499 = 2500 BC (default: -2499)
        month      - 1-12 (default: 3)
        day        - 1 to the length of the month (default: 20)
        hour       - local mean solar time, 0-24 (default: 22.0)
        calendar   - 'gregorian' (default) or 'julian', both proleptic
        refraction - 1 for apparent altitudes (default), 0 for geometric

    meta.dates gives the same local date in both calendars; ecliptic is the same
    list /api/ecliptic returns.
    Invalid parameters return HTTP 400 (404 for an unknown site) with {"error": message}.
    """
    lat, lon, elevation, site = _location()
    year, month, day, hour, calendar = _date_time()
    refract = _flag('refraction', True)
    results = calculate_alignments(lat, lon, year, month, day, hour, elevation, refract, calendar)
    ecliptic_points = calculate_ecliptic(lat, lon, year, month, day, hour, elevation, refract, calendar)

    monument_info = None
    if site:
        monument_info = {
            'name':           site,
            'orientation_az': MONUMENTS[site]['orientation_az'],
            'note':           MONUMENTS[site]['note'],
        }

    return jsonify({
        'meta': {
            'lat':         lat,
            'lon':         lon,
            'elevation_m': elevation,
            'refraction':  refract,
            'year':    year,
            'era':     'BC' if year <= 0 else 'AD',
            'month':   month,
            'day':     day,
            'calendar': calendar,
            'dates':   local_calendar_dates(results['jd'], lon),
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
        'ecliptic': ecliptic_points,
    })


@app.route('/api/heliacal')
def heliacal():
    """
    Find the heliacal rising of a star in a given year and location.

    Query parameters: site or lat/lon/elevation, refraction and calendar (as
    /api/stars), year (astronomical, in that calendar), star (catalog name,
    default Sirius), extinction (atmospheric extinction in magnitudes per
    airmass, 0.1 to 0.6, default 0.27; higher means hazier air, so the star
    must climb higher to be seen). The date returned is in the same calendar.
    """
    lat, lon, elevation, _ = _location()
    calendar = _calendar()
    year = _number('year', -2780, int, *YEAR_RANGE)
    extinction = _number('extinction', DEFAULT_EXTINCTION, float, 0.1, 0.6)
    refract = _flag('refraction', True)
    star = request.args.get('star', 'Sirius')
    if star not in STARS:
        raise ApiError("Unknown star; it must be one of the catalog names")

    result = find_heliacal_rising(lat, lon, year, star, extinction, elevation, refract, calendar)
    result['calendar'] = calendar
    if not result['found']:
        when = format_year(year)
        result['message'] = {
            'always_visible': f"{star} is visible at every dawn here in {when}: it never sinks "
                              "out of sight, so it has no heliacal rising.",
            'never_visible': f"{star} never climbs high enough in a dark enough sky to be seen "
                             f"at dawn here in {when}.",
        }.get(result['reason'], f"No heliacal rising of {star} found here in {when}.")

    return jsonify(result)


@app.route('/api/alignment')
def alignment():
    """
    Find when a star, or the Sun, rose or set along an axis.

    Query parameters: site or lat/lon/elevation, refraction and calendar (as
    /api/stars); star (catalog name or 'Sun', default Sirius); event ('rise'
    or 'set', default 'rise'); azimuth (degrees from north through east;
    defaults to the monument's orientation, required with custom coordinates;
    0-180 for a rising, 180-360 for a setting); horizon (apparent altitude of
    the visible horizon, -2 to 10 degrees; defaults to the monument's measured
    horizon where one is known, else 0); from and to (astronomical years, default 12,000 BC to the current
    year). For the Sun also: year (astronomical, default -2499) to list the
    days of that year, and limb ('upper' for the first gleam, the default,
    'centre' or 'lower' for the whole disc on the horizon).
    """
    lat, lon, elevation, site = _location()
    calendar = _calendar()
    refract = _flag('refraction', True)
    star = request.args.get('star', 'Sirius')
    if star != 'Sun' and star not in STARS:
        raise ApiError("Unknown star; it must be 'Sun' or one of the catalog names")
    event = request.args.get('event', 'rise').lower()
    if event not in EVENTS:
        raise ApiError("'event' must be 'rise' or 'set'")

    default_azimuth = MONUMENTS[site]['orientation_az'] if site else None
    azimuth = _number('azimuth', default_azimuth, float, 0.0, 360.0)
    if azimuth is None:
        raise ApiError("'azimuth' is required with custom coordinates")
    if azimuth % 180.0 == 0.0:
        raise ApiError("An axis due north or due south meets stars at their highest or lowest "
                       "point in the sky, never as they rise or set; try the rising or setting "
                       "direction of the star you have in mind")
    low, high, side = (0.0, 180.0, 'eastern') if event == 'rise' else (180.0, 360.0, 'western')
    if not low < azimuth < high:
        raise ApiError(f"Stars and the Sun {event} on the {side} half of the horizon, so the azimuth "
                       f"must be between {low:.0f} and {high:.0f} degrees; got {azimuth:g}")

    site_horizon = MONUMENTS[site]['horizon_deg'] if site else None
    horizon = _number('horizon', 0.0 if site_horizon is None else site_horizon, float, -2.0, 10.0)
    first_year = _number('from', -11999, int, *YEAR_RANGE)
    last_year = _number('to', datetime.date.today().year, int, *YEAR_RANGE)
    if last_year <= first_year:
        raise ApiError("'to' must be later than 'from'")

    if star == 'Sun':
        return jsonify(_sun_alignment(lat, lon, elevation, azimuth, event, horizon, refract,
                                      first_year, last_year, calendar))

    result = find_axis_crossings(lat, lon, star, azimuth, event, horizon, elevation,
                                 refract, first_year, last_year, calendar)
    result.update({'star': star, 'event': event, 'azimuth': azimuth, 'horizon': horizon,
                   'calendar': calendar})
    if not result['found']:
        span = f"between {format_year(first_year)} and {format_year(last_year)}"
        if result.get('reason') == 'never_crosses_horizon':
            result['message'] = (f"{star} never {event}s over this horizon {span}: it stays "
                                 "either always above it or always below it.")
        else:
            c = result['closest']
            result['message'] = (f"{star} did not {event} on the {azimuth:g} degree axis {span}; "
                                 f"it came closest in {format_year(c['year'])}, "
                                 f"{abs(c['offset']):.1f} degrees away.")
    return jsonify(result)


def _sun_alignment(lat, lon, elevation, azimuth, event, horizon, refract, first_year, last_year,
                   calendar):
    year = _number('year', -2499, int, *YEAR_RANGE)
    limb = request.args.get('limb', 'upper').lower()
    if limb not in SUN_LIMBS:
        raise ApiError("'limb' must be 'upper', 'centre' or 'lower'")

    result = find_sun_alignment(lat, lon, azimuth, event, year, horizon, elevation, refract,
                                limb, first_year, last_year, calendar)
    result.update({'star': 'Sun', 'event': event, 'azimuth': azimuth, 'horizon': horizon,
                   'limb': limb, 'calendar': calendar})
    when = 'sunrise' if event == 'rise' else 'sunset'
    if result['june'] is None:
        result['message'] = f"There is no {when} here in {format_year(year)}."
    elif not result['dates']:
        result['message'] = (f"In {format_year(year)} {when} never fell on the {azimuth:g} degree "
                             f"axis: it ranged from {result['june']['azimuth']:.1f} degrees at the "
                             f"June solstice to {result['december']['azimuth']:.1f} degrees at the "
                             "December solstice.")
    return result


@app.route('/api/ecliptic')
def ecliptic():
    """
    Return the ecliptic great circle (73 points) in local alt-az.

    Accepts the same query parameters as /api/stars.
    """
    lat, lon, elevation, _ = _location()
    year, month, day, hour, calendar = _date_time()
    refract = _flag('refraction', True)
    return jsonify({'points': calculate_ecliptic(lat, lon, year, month, day, hour,
                                                 elevation, refract, calendar)})


@app.route('/api/sites')
def sites():
    """Return the list of all known monument sites."""
    return jsonify({
        name: {
            'lat':            data['lat'],
            'lon':            data['lon'],
            'elevation_m':    data['elevation_m'],
            'orientation_az': data['orientation_az'],
            'horizon_deg':    data['horizon_deg'],
            'source':         data['source'],
            'note':           data['note'],
        }
        for name, data in MONUMENTS.items()
    })


if __name__ == '__main__':
    # The Werkzeug debugger can run arbitrary code, so it is opt-in: FLASK_DEBUG=1.
    app.run(debug=os.environ.get('FLASK_DEBUG') == '1')
