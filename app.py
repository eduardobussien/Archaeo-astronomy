import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from flask import Flask, jsonify, request, render_template
from alignment import calculate_alignments, calculate_ecliptic, find_heliacal_rising
from data import MONUMENTS, STARS

app = Flask(__name__)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/stars')
def stars():
    """
    Return star positions for a given location and historical date.

    Query parameters:
        lat   - latitude in degrees (default: Giza)
        lon   - longitude in degrees (default: Giza)
        year  - astronomical year, negative = BC (default: -2500)
        month - 1-12 (default: 3)
        day   - 1-31 (default: 20)
        hour  - 0-23.99 (default: 22.0)
        site  - optional monument name (overrides lat/lon)

    Returns JSON with star positions and metadata.
    """
    lat   = float(request.args.get('lat',   29.9792))
    lon   = float(request.args.get('lon',   31.1342))
    year  = int(request.args.get('year',   -2500))
    month = int(request.args.get('month',  3))
    day   = int(request.args.get('day',    20))
    hour  = float(request.args.get('hour', 22.0))
    site  = request.args.get('site', None)

    monument_info = None
    if site:
        match = next((k for k in MONUMENTS if site.lower() in k.lower()), None)
        if match:
            lat = MONUMENTS[match]['lat']
            lon = MONUMENTS[match]['lon']
            monument_info = {
                'name':           match,
                'orientation_az': MONUMENTS[match]['orientation_az'],
                'note':           MONUMENTS[match]['note'],
            }
        else:
            return jsonify({'error': f"Site '{site}' not found"}), 404

    results = calculate_alignments(lat, lon, year, month, day, hour)

    era = 'BC' if year < 0 else 'AD'
    return jsonify({
        'meta': {
            'lat':    lat,
            'lon':    lon,
            'year':   year,
            'era':    era,
            'month':  month,
            'day':    day,
            'hour':    hour,
            'hour_ut': round(hour - lon / 15.0, 4),
            'jd':      float(results['jd']),
            'lst':    float(results['lst']),
            'method': results['method'],
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
            for name, d in results.get('planets', {}).items()
        },
    })


@app.route('/api/heliacal')
def heliacal():
    """
    Find the heliacal rising of a star in a given year and location.

    Query params: lat, lon, year, star, site (optional), arc_vision (default -10).
    """
    lat        = float(request.args.get('lat',        29.9792))
    lon        = float(request.args.get('lon',        31.1342))
    year       = int(request.args.get('year',        -2780))
    star       = request.args.get('star', 'Sirius')
    arc_vision = float(request.args.get('arc_vision', -10.0))
    site       = request.args.get('site', None)

    if site:
        match = next((k for k in MONUMENTS if site.lower() in k.lower()), None)
        if match:
            lat = MONUMENTS[match]['lat']
            lon = MONUMENTS[match]['lon']

    if star not in STARS:
        return jsonify({'found': False, 'message': f"Unknown star: '{star}'"}), 400

    result = find_heliacal_rising(lat, lon, year, star, arc_vision)
    if result is None:
        era = 'BC' if year < 0 else 'AD'
        return jsonify({
            'found': False,
            'message': f"No heliacal rising of {star} found at this location in {abs(year)} {era}.",
        })

    return jsonify({'found': True, **result})


@app.route('/api/ecliptic')
def ecliptic():
    """
    Return the ecliptic great circle (73 points) in local alt-az.

    Accepts the same query parameters as /api/stars (lat, lon, year, month, day, hour, site).
    """
    lat   = float(request.args.get('lat',   29.9792))
    lon   = float(request.args.get('lon',   31.1342))
    year  = int(request.args.get('year',   -2500))
    month = int(request.args.get('month',  3))
    day   = int(request.args.get('day',    20))
    hour  = float(request.args.get('hour', 22.0))
    site  = request.args.get('site', None)

    if site:
        match = next((k for k in MONUMENTS if site.lower() in k.lower()), None)
        if match:
            lat = MONUMENTS[match]['lat']
            lon = MONUMENTS[match]['lon']

    points = calculate_ecliptic(lat, lon, year, month, day, hour)
    return jsonify({'points': points})


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
    app.run(debug=True)
