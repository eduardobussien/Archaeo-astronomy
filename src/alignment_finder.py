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
from precession import mean_obliquity, star_altaz, star_declination
from refraction import true_altitude
from solar_system import sun_altaz

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
    return {**_local_moment(jd_event, lon, calendar), 'azimuth': round(float(azimuth), 3)}


def _local_moment(jd, lon, calendar):
    """{'month', 'day', 'hour'} in local mean time for a UT Julian Date."""
    local = jd + lon / 360.0
    _, month, day = _jd_to_date(local, calendar)
    return {'month': month, 'day': day, 'hour': round(float((local + 0.5) % 1.0) * 24.0, 3)}


def _axis_summary(years, az, azimuth):
    """
    Where a yearly azimuth series crosses the axis, with its closest approach,
    its last defined value and a thinned curve for plotting; None if the
    series is undefined throughout.
    """
    offset = az - azimuth
    valid = np.isfinite(offset)
    if not valid.any():
        return None
    with np.errstate(invalid='ignore'):
        changes = valid[:-1] & valid[1:] & (np.sign(offset[:-1]) != np.sign(offset[1:]))
    crossings = []
    for i in np.flatnonzero(changes):
        f = offset[i] / (offset[i] - offset[i + 1])
        crossings.append({
            'year': int(round(years[i] + f)),
            'drift_deg_per_century': round(float(az[i + 1] - az[i]) * 100.0, 3),
        })
    k = int(np.nanargmin(np.abs(offset)))
    latest = int(np.flatnonzero(valid)[-1])
    return {
        'crossings': crossings,
        'closest': {'year': int(years[k]), 'azimuth': round(float(az[k]), 3),
                    'offset': round(float(offset[k]), 3)},
        'latest': {'year': int(years[latest]), 'azimuth': round(float(az[latest]), 3)},
        'curve': [{'year': int(y), 'azimuth': None if np.isnan(a) else round(float(a), 3)}
                  for y, a in zip(years[::_CURVE_STEP_YEARS], az[::_CURVE_STEP_YEARS])],
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
    summary = _axis_summary(years, event_azimuths(lat, star_name, years, event, horizon,
                                                  elevation_m, refract), azimuth)
    if summary is None:
        return {'found': False, 'reason': 'never_crosses_horizon'}
    for c in summary['crossings']:
        c['night'] = night_of_event(lat, lon, c['year'], star_name, event, horizon,
                                    elevation_m, refract, calendar)
    return {'found': bool(summary['crossings']), **summary}


# ── The Sun ──
# Its rising point swings through the year between the solstices rather than
# drifting over the centuries, so it gets two searches: the dates in one year
# when sunrise (or sunset) falls on the axis, and, since the solstice extremes
# shift as the obliquity changes, the years in which a solstice sunrise did.

SUN_LIMBS = ('upper', 'centre', 'lower')
_SUN_SEMIDIAMETER = 0.2666   # mean apparent radius, degrees
_LIMB_OFFSET = {'upper': -_SUN_SEMIDIAMETER, 'centre': 0.0, 'lower': _SUN_SEMIDIAMETER}


def _sun_centre_altitude(horizon, limb, elevation_m, refract):
    """
    Geometric altitude of the Sun's centre when the chosen limb touches the
    visible horizon. Refraction applies where the limb is seen, at the horizon;
    the centre lies one semidiameter from the limb's geometric position. (Refracting
    the centre instead, 16' lower where refraction is stronger, put first gleam
    0.055 deg too low and sunrise about 20 s early.)
    """
    limb_altitude = true_altitude(horizon, elevation_m) if refract else horizon
    return limb_altitude + _LIMB_OFFSET[limb]


def sun_events_in_year(lat, lon, year, event='rise', horizon=0.0, elevation_m=0.0, refract=True,
                       limb='upper', calendar='gregorian'):
    """
    Every sunrise (or sunset) of a year from the full solar ephemeris, by
    bisection over all days at once: (UT Julian Dates, azimuths in degrees),
    NaN on days without one (polar day or night). limb is which part of the
    disc touches the horizon: 'upper' is the first (or last) gleam.
    """
    h0 = _sun_centre_altitude(horizon, limb, elevation_m, refract)
    first = _date_to_jd(year, 1, 1, -lon / 15.0, calendar)
    midnights = first + np.arange(round(_date_to_jd(year + 1, 1, 1, -lon / 15.0, calendar) - first))
    rising = event == 'rise'
    lo = np.full_like(midnights, 0.0 if rising else 12.0)
    hi = lo + 12.0
    alt_lo, alt_hi = (sun_altaz(midnights + h / 24.0, lat, lon)[0] for h in (lo, hi))
    has_event = (alt_lo < h0) & (alt_hi > h0) if rising else (alt_lo > h0) & (alt_hi < h0)
    for _ in range(20):
        mid = (lo + hi) / 2.0
        later = (sun_altaz(midnights + mid / 24.0, lat, lon)[0] < h0) == rising
        lo = np.where(later, mid, lo)
        hi = np.where(later, hi, mid)
    jd = midnights + (lo + hi) / 48.0
    az = sun_altaz(jd, lat, lon)[1]
    return np.where(has_event, jd, np.nan), np.where(has_event, az, np.nan)


def find_sun_alignment(lat, lon, azimuth, event='rise', year=-2499, horizon=0.0, elevation_m=0.0,
                       refract=True, limb='upper', first_year=-11999, last_year=2026,
                       calendar='gregorian'):
    """
    When sunrise (or sunset) fell on the axis.

    'dates': the days of the given year on which it did, each with its local
    time; 'june' and 'december': where the solstice sunrises fell that year;
    'solstice': which solstice lies on the axis's side of east (or west) and
    the years in which it fell exactly on the axis, as find_axis_crossings
    reports for stars, each with that year's solstice sunrise to show.
    """
    jd, az = sun_events_in_year(lat, lon, year, event, horizon, elevation_m, refract, limb, calendar)
    offset = az - azimuth
    with np.errstate(invalid='ignore'):
        changes = (np.isfinite(offset[:-1]) & np.isfinite(offset[1:])
                   & (np.sign(offset[:-1]) != np.sign(offset[1:])))
    dates = []
    for i in np.flatnonzero(changes):
        k = i if abs(offset[i]) <= abs(offset[i + 1]) else i + 1
        dates.append({**_local_moment(jd[k], lon, calendar), 'azimuth': round(float(az[k]), 3)})

    def extreme(k):
        return {**_local_moment(jd[k], lon, calendar), 'azimuth': round(float(az[k]), 3)}

    if np.isfinite(az).any():
        northmost, southmost = (int(np.nanargmin(az)), int(np.nanargmax(az))) if event == 'rise' \
            else (int(np.nanargmax(az)), int(np.nanargmin(az)))
        june, december = extreme(northmost), extreme(southmost)
    else:
        june = december = None

    north = azimuth < 90.0 if event == 'rise' else azimuth > 270.0
    years = np.arange(first_year, last_year + 1)
    epochs = np.array([_date_to_jd(int(y), 7, 2, 0.0) for y in years])
    declination = mean_obliquity(epochs) * (1.0 if north else -1.0)
    rise = horizon_azimuth(declination, lat, _sun_centre_altitude(horizon, limb, elevation_m, refract))
    summary = _axis_summary(years, rise if event == 'rise' else 360.0 - rise, azimuth)
    if summary is not None:
        for c in summary['crossings']:
            c_jd, c_az = sun_events_in_year(lat, lon, c['year'], event, horizon, elevation_m,
                                            refract, limb, calendar)
            k = int(np.nanargmin(c_az) if north == (event == 'rise') else np.nanargmax(c_az))
            c['night'] = {**_local_moment(c_jd[k], lon, calendar), 'azimuth': round(float(c_az[k]), 3)}

    return {
        'found': bool(dates) or bool(summary and summary['crossings']),
        'year': year,
        'dates': dates,
        'june': june,
        'december': december,
        'solstice': {'which': 'june' if north else 'december', **(summary or {'crossings': []})},
    }
