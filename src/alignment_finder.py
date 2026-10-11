"""
When did a star rise or set along a monument's axis?

A star of declination d meets a horizon at geometric altitude h0 at azimuth

    cos A = (sin d - sin phi sin h0) / (cos phi cos h0)

for an observer at latitude phi: it rises at A (measured from north through
east) and sets at 360 - A. h0 is the visible horizon (0 deg, or the height of
distant hills) lowered by refraction, since a star first shows while still
geometrically below it. Precession and proper motion change d over the
centuries, so the rising point drifts along the horizon; the search evaluates
it for every year and reports the years in which it crosses the axis.
"""
import numpy as np

from alignment import _date_to_jd, _jd_to_date, _observed
from data import STARS
from precession import star_altaz, star_declination
from refraction import true_altitude

EVENTS = ('rise', 'set')
_SIDEREAL_SHIFT_HOURS_PER_DAY = 24.0 / 365.2422
_CURVE_STEP_YEARS = 50


def horizon_azimuth(declination, lat, horizon_altitude):
    """Rising azimuth in degrees; NaN where the star never crosses that horizon."""
    d, phi, h = map(np.radians, (declination, lat, horizon_altitude))
    cos_a = (np.sin(d) - np.sin(phi) * np.sin(h)) / (np.cos(phi) * np.cos(h))
    with np.errstate(invalid='ignore'):
        return np.degrees(np.arccos(np.where(np.abs(cos_a) <= 1.0, cos_a, np.nan)))


def _star_args(name):
    s = STARS[name]
    return s['ra'], s['dec'], s['pm_ra'], s['pm_dec'], s['dist']


def event_azimuths(lat, star_name, years, event='rise', horizon=0.0, elevation_m=0.0, refract=True):
    """Azimuth (degrees) at which the star rises or sets in each year (mid-year epoch)."""
    jd = np.array([_date_to_jd(int(y), 7, 2, 0.0) for y in years])
    dec = star_declination(*_star_args(star_name), STARS[star_name]['rv'], jd)
    h0 = true_altitude(horizon, elevation_m) if refract else horizon
    rise = horizon_azimuth(dec, lat, h0)
    return rise if event == 'rise' else 360.0 - rise


def night_of_event(lat, lon, year, star_name, event='rise', horizon=0.0, elevation_m=0.0,
                   refract=True, calendar='gregorian'):
    """
    A date in the year when the star rises (or sets) close to local midnight,
    and the local time it does so: {'month', 'day', 'hour', 'azimuth'}.

    The event comes about 3.94 minutes earlier each day, so the day on which
    it falls at midnight follows from its time on 1 January.
    """
    s = STARS[star_name]
    first_midnight = _date_to_jd(year, 1, 1, -lon / 15.0, calendar)

    def crossing_near(midnight_jd):
        jd = midnight_jd + np.linspace(-0.5, 0.5, 1441)
        alt, az = star_altaz(*_star_args(star_name), jd, lat, lon, s['rv'])
        alt = _observed(alt, elevation_m, refract) - horizon
        above = alt >= 0
        edges = np.flatnonzero(above[1:] != above[:-1])
        edges = [i for i in edges if above[i + 1] == (event == 'rise')]
        if not edges:
            return None
        i = min(edges, key=lambda k: abs(jd[k] - midnight_jd))
        f = alt[i] / (alt[i] - alt[i + 1])
        azimuth = az[i] + f * ((az[i + 1] - az[i] + 180.0) % 360.0 - 180.0)
        return jd[i] + f * (jd[i + 1] - jd[i]), azimuth % 360.0

    first = crossing_near(first_midnight)
    if first is None:
        return None
    hours_after_midnight = ((first[0] - first_midnight) * 24.0) % 24.0
    day = round(hours_after_midnight / _SIDEREAL_SHIFT_HOURS_PER_DAY) % 365
    found = crossing_near(first_midnight + day)
    if found is None:
        return None
    jd_event, azimuth = found
    local = jd_event + lon / 360.0
    _, month, day_of_month = _jd_to_date(local, calendar)
    return {
        'month': month,
        'day': day_of_month,
        'hour': round(float((local + 0.5) % 1.0) * 24.0, 3),
        'azimuth': round(float(azimuth), 3),
    }


def find_axis_crossings(lat, lon, star_name, azimuth, event='rise', horizon=0.0, elevation_m=0.0,
                        refract=True, first_year=-11999, last_year=2026, calendar='gregorian'):
    """
    Years (astronomical numbering) in which the star rose or set at the given
    azimuth, each with the drift of its rising point and a night to see it.

    Returns {'found', 'crossings', 'closest', 'latest', 'curve'} or, if the
    star never meets that horizon in the period, {'found': False,
    'reason': 'never_crosses_horizon'}.
    """
    years = np.arange(first_year, last_year + 1)
    az = event_azimuths(lat, star_name, years, event, horizon, elevation_m, refract)
    offset = az - azimuth
    valid = np.isfinite(offset)
    if not valid.any():
        return {'found': False, 'reason': 'never_crosses_horizon'}

    with np.errstate(invalid='ignore'):
        changes = valid[:-1] & valid[1:] & (np.sign(offset[:-1]) != np.sign(offset[1:]))
    crossings = []
    for i in np.flatnonzero(changes):
        f = offset[i] / (offset[i] - offset[i + 1])
        year = int(round(years[i] + f))
        crossings.append({
            'year': year,
            'drift_deg_per_century': round(float(az[i + 1] - az[i]) * 100.0, 3),
            'night': night_of_event(lat, lon, year, star_name, event, horizon,
                                    elevation_m, refract, calendar),
        })

    k = int(np.nanargmin(np.abs(offset)))
    latest = int(np.flatnonzero(valid)[-1])
    curve = [{'year': int(y), 'azimuth': None if np.isnan(a) else round(float(a), 3)}
             for y, a in zip(years[::_CURVE_STEP_YEARS], az[::_CURVE_STEP_YEARS])]
    return {
        'found': bool(crossings),
        'crossings': crossings,
        'closest': {'year': int(years[k]), 'azimuth': round(float(az[k]), 3),
                    'offset': round(float(offset[k]), 3)},
        'latest': {'year': int(years[latest]), 'azimuth': round(float(az[latest]), 3)},
        'curve': curve,
    }
