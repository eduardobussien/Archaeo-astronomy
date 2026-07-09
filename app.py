import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from flask import Flask, jsonify, request
from alignment import calculate_alignments
from data import MONUMENTS

app = Flask(__name__)


@app.route('/')
def index():
    return jsonify({'status': 'ok', 'message': 'Archaeo-Astronomy API is running'})


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
            'hour':   hour,
            'jd':     float(results['jd']),
            'lst':    float(results['lst']),
            'method': results['method'],
        },
        'monument': monument_info,
        'stars': {
            name: {
                'altitude': float(d['altitude']),
                'azimuth':  float(d['azimuth']),
                'visible':  bool(d['visible']),
            }
            for name, d in results['stars'].items()
        },
    })


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
