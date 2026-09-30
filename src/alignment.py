import argparse
import math
import warnings
import numpy as np
import erfa

from data import STARS, MONUMENTS
from precession import star_altaz, ecliptic_altaz, local_mean_sidereal_time
from refraction import apparent_altitude
from solar_system import body_altaz, sun_altitude
from visibility import DEFAULT_EXTINCTION, visibility_thresholds

warnings.simplefilter('ignore', category=erfa.ErfaWarning)

_STAR_NAMES = list(STARS)
_CATALOG = {
    key: np.array([STARS[name][key] for name in _STAR_NAMES])
    for key in ('ra', 'dec', 'pm_ra', 'pm_dec', 'dist', 'rv')
}


def _date_to_jd(year, month, day, hour):
    """
    Convert a proleptic Gregorian calendar date to Julian Day Number.

    Uses the Meeus algorithm (Astronomical Algorithms, ch. 7).
    Python's // (floor division) matches Meeus's INT() for negative years.
    Year 0 = 1 BC, year -1 = 2 BC, etc. (astronomical year numbering).
    """
    y, m = year, month
    d = day + hour / 24.0
    if m <= 2:
        y -= 1
        m += 12
    a = y // 100
    b = 2 - a + a // 4
    return math.floor(365.25 * (y + 4716)) + math.floor(30.6001 * (m + 1)) + d + b - 1524.5


def days_in_month(year, month):
    """Length of a proleptic Gregorian month (astronomical year numbering)."""
    if month == 2:
        return 29 if (year % 4 == 0 and year % 100 != 0) or year % 400 == 0 else 28
    return 30 if month in (4, 6, 9, 11) else 31


def format_year(year):
    """Label an astronomical year (0 = 1 BC, -1 = 2 BC) in historical form."""
    return f"{1 - year} BC" if year <= 0 else f"{year} AD"


_DAYS_PER_400_YEARS = 146097


def _jd_to_date(jd):
    """
    Proleptic Gregorian (year, month, day) containing a Julian Date.

    Inverse of _date_to_jd (Meeus ch. 7). Negative Julian Dates are first
    shifted forward by whole 400-year Gregorian cycles, over which the
    calendar repeats exactly.
    """
    cycles = max(0, math.ceil(-jd / _DAYS_PER_400_YEARS) + 1)
    z = math.floor(jd + 0.5) + cycles * _DAYS_PER_400_YEARS
    alpha = math.floor((z - 1867216.25) / 36524.25)
    b = z + 1 + alpha - math.floor(alpha / 4) + 1524
    c = math.floor((b - 122.1) / 365.25)
    e = math.floor((b - math.floor(365.25 * c)) / 30.6001)
    day = b - math.floor(365.25 * c) - math.floor(30.6001 * e)
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715
    return year - 400 * cycles, month, day


def _observed(altitude, elevation_m, refract):
    """Geometric altitude as an observer sees it, with refraction if requested."""
    return apparent_altitude(altitude, elevation_m) if refract else np.asarray(altitude, dtype=float)


def calculate_alignments(lat, lon, year, month, day, hour, elevation_m=0.0, refract=True):
    """
    Compute altitude and azimuth for every star, the Sun, Moon and planets.

    Stars use the long-term precession engine in precession.py; the Sun, Moon
    and planets come from solar_system.py and share the same Earth rotation.
    Altitudes include atmospheric refraction for the site's elevation unless
    refract is False.

    Returns a dict with keys:
        jd      - Julian Date (UT1)
        lst     - Local mean sidereal time in degrees
        stars   - {name: {altitude, azimuth, visible}} for each star
        planets - {name: {altitude, azimuth, visible}} for the Sun, Moon and planets
        method  - name of the star position model
    """
    # Convert local mean solar time to UT (the user picks local time; lon/15 is the offset)
    hour_ut = hour - lon / 15.0
    jd = _date_to_jd(year, month, day, hour_ut)

    alt, az = star_altaz(_CATALOG['ra'], _CATALOG['dec'], _CATALOG['pm_ra'],
                         _CATALOG['pm_dec'], _CATALOG['dist'], jd, lat, lon, _CATALOG['rv'])
    alt = _observed(alt, elevation_m, refract)
    stars_out = {
        name: {'altitude': float(a), 'azimuth': float(z), 'visible': bool(a > 0)}
        for name, a, z in zip(_STAR_NAMES, alt, az)
    }

    planets = {}
    for name, (a, z) in body_altaz(jd, lat, lon).items():
        a = float(_observed(a, elevation_m, refract))
        planets[name] = {'altitude': a, 'azimuth': z, 'visible': a > 0}

    return {
        'jd': jd,
        'lst': float(local_mean_sidereal_time(jd, lon)),
        'stars': stars_out,
        'planets': planets,
        'method': 'Vondrák 2011',
    }


def _dawn_hours(midnights, lat, lon, sun_limit):
    """
    Local hour (0-14) at which the Sun rises through sun_limit on each day,
    by bisection over all days at once. NaN where it never crosses (polar
    day or night).
    """
    lo = np.zeros_like(midnights)
    hi = np.full_like(midnights, 14.0)
    crosses = ((sun_altitude(midnights, lat, lon) < sun_limit)
               & (sun_altitude(midnights + hi / 24.0, lat, lon) > sun_limit))
    for _ in range(20):
        mid = (lo + hi) / 2.0
        below = sun_altitude(midnights + mid / 24.0, lat, lon) < sun_limit
        lo = np.where(below, mid, lo)
        hi = np.where(below, hi, mid)
    return np.where(crosses, (lo + hi) / 2.0, np.nan)


def find_heliacal_rising(lat, lon, year, star_name, extinction=DEFAULT_EXTINCTION,
                         elevation_m=0.0, refract=True):
    """
    Find the first dawn in the target year on which star_name can be seen,
    after a dawn on which it could not.

    A dawn counts when, at the last moment the sky is dark enough for a star
    of this magnitude (the Sun at its limit from visibility.py), the star
    already stands at the altitude it needs to shine through the extinction
    k (magnitudes per airmass). The star's altitude includes refraction
    unless refract is False; the Sun's is geometric.

    The scan starts on 1 October of the previous year so that a star already
    visible on 1 January is not reported as rising that day.

    Returns {'found': True, date and geometry...}, or {'found': False,
    'reason': ...} where reason is 'always_visible' (seen every dawn, as for a
    star that never sets), 'never_visible' (never high enough in a dark
    enough sky) or 'not_this_year'.
    """
    s = STARS[star_name]
    star_limit, sun_limit = (float(v) for v in visibility_thresholds(s['mag'], extinction))
    first_midnight = _date_to_jd(year - 1, 10, 1, -lon / 15.0)
    n_days = round(_date_to_jd(year + 1, 1, 1, -lon / 15.0) - first_midnight)
    first_target_day = round(_date_to_jd(year, 1, 1, -lon / 15.0) - first_midnight)
    midnights = first_midnight + np.arange(n_days)

    dawn = _dawn_hours(midnights, lat, lon, sun_limit)
    has_dawn = ~np.isnan(dawn)
    dawn_jd = midnights + np.where(has_dawn, dawn, 0.0) / 24.0
    star_alt, star_az = star_altaz(s['ra'], s['dec'], s['pm_ra'], s['pm_dec'], s['dist'],
                                   dawn_jd, lat, lon, s['rv'])
    star_alt = _observed(star_alt, elevation_m, refract)
    visible = has_dawn & (star_alt >= star_limit)

    rising = np.flatnonzero(visible[1:] & ~visible[:-1]) + 1
    rising = rising[rising >= first_target_day]
    if rising.size == 0:
        dark_dawns = visible[first_target_day:][has_dawn[first_target_day:]]
        if dark_dawns.size and dark_dawns.all():
            reason = 'always_visible'
        elif not dark_dawns.any():
            reason = 'never_visible'
        else:
            reason = 'not_this_year'
        return {'found': False, 'reason': reason}

    k = rising[0]
    y, m, d = _jd_to_date(midnights[k] + 0.5)
    return {
        'found': True,
        'year':  y,
        'month': m,
        'day':   d,
        'star_altitude':   round(float(star_alt[k]), 2),
        'star_azimuth':    round(float(star_az[k]), 2),
        'sun_altitude':    round(float(sun_altitude(dawn_jd[k], lat, lon)), 2),
        'dawn_hour_local': round(float(dawn[k]), 2),
        'required_star_altitude': round(star_limit, 2),
        'extinction':      extinction,
    }


def calculate_ecliptic(lat, lon, year, month, day, hour, elevation_m=0.0, refract=True):
    """
    Return 73 points (0°..360° ecliptic longitude, step 5°) projected onto
    the local alt-az frame.  The 73rd point closes the loop back to 0°.

    Uses the mean ecliptic of date from the long-term precession model, with
    refraction applied as in calculate_alignments.
    """
    jd = _date_to_jd(year, month, day, hour - lon / 15.0)
    longitudes = np.arange(73) * 5.0
    alt, az = ecliptic_altaz(jd, lat, lon, longitudes)
    alt = _observed(alt, elevation_m, refract)
    return [
        {'longitude': float(lam), 'altitude': float(a), 'azimuth': float(z)}
        for lam, a, z in zip(longitudes, alt, az)
    ]


def check_alignments(star_results, orientation_az, threshold_deg=2.0):
    """
    Return stars whose azimuth falls within threshold_deg of orientation_az.

    orientation_az - compass bearing (0=N, 90=E, 180=S, 270=W) in degrees
    threshold_deg  - angular tolerance for a match

    The angular difference wraps correctly (e.g. 359° vs 1° = 2° apart).
    """
    matches = []
    for name, data in star_results['stars'].items():
        diff = abs(data['azimuth'] - orientation_az)
        if diff > 180:
            diff = 360 - diff
        if diff <= threshold_deg:
            matches.append({
                'star':           name,
                'azimuth':        round(data['azimuth'], 2),
                'orientation_az': round(orientation_az, 2),
                'diff_deg':       round(diff, 2),
                'altitude':       round(data['altitude'], 2),
                'visible':        data['visible'],
            })
    return sorted(matches, key=lambda x: x['diff_deg'])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _print_star_table(results):
    engine = results['method']
    print(f"\n=== Archaeo-Astronomy Alignment [{engine}] ===")
    print(f"Location : {results['_lat']:.4f}, {results['_lon']:.4f}")
    print(f"Date     : {format_year(results['_year'])}, {results['_month']:02d}-{results['_day']:02d}")
    print(f"Time     : {results['_hour']:05.2f}h  |  JD: {results['jd']:.2f}  |  LST: {results['lst']:.2f}")
    print("-" * 62)
    print(f"{'STAR':<12s} | {'CONST':<14s} | {'ALT':>8s} | {'AZ':>8s} | STATUS")
    print("-" * 62)
    for name, d in results['stars'].items():
        status = "VISIBLE" if d['visible'] else "hidden"
        const  = STARS[name]['constellation']
        print(f"{name:<12s} | {const:<14s} | {d['altitude']:6.2f}°   | {d['azimuth']:6.2f}°   | {status}")
    print()


def _print_alignment_results(matches, monument_name, threshold_deg):
    print(f"\n--- Alignment check vs '{monument_name}' (±{threshold_deg}°) ---")
    if not matches:
        print("  No stars within the threshold at this date/time.\n")
        return
    for m in matches:
        vis = "VISIBLE" if m['visible'] else "hidden"
        print(
            f"  {m['star']:<12s}  az={m['azimuth']:6.2f}°  "
            f"(Δ={m['diff_deg']:.2f}°)  alt={m['altitude']:+.2f}°  [{vis}]"
        )
    print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Archaeo-Astronomy Alignment Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  python src/alignment.py --monument Stonehenge --year -2499   (2500 BC)\n"
               "  python src/alignment.py --list-monuments",
    )
    parser.add_argument("--lat",    type=float, default=29.9792)
    parser.add_argument("--lon",    type=float, default=31.1342)
    parser.add_argument("--year",   type=int,   default=-2499,
                        help="Astronomical year: 0 = 1 BC, -2499 = 2500 BC")
    parser.add_argument("--month",  type=int,   default=3)
    parser.add_argument("--day",    type=int,   default=20)
    parser.add_argument("--hour",   type=float, default=0.0,
                        help="Hour in 24h format (decimals OK, e.g. 22.5 = 22:30)")
    parser.add_argument("--monument", type=str, default=None,
                        help="Monument name (overrides --lat/--lon). Use --list-monuments to see options.")
    parser.add_argument("--threshold", type=float, default=2.0,
                        help="Azimuth tolerance in degrees for alignment check (default: 2.0)")
    parser.add_argument("--list-monuments", action="store_true",
                        help="Print all available monuments and exit")

    args = parser.parse_args()

    if args.list_monuments:
        print("\nAvailable monuments:")
        for name, m in MONUMENTS.items():
            print(f"  \"{name}\"")
            print(f"    lat={m['lat']}, lon={m['lon']}, orientation={m['orientation_az']}°")
            print(f"    {m['note']}\n")
        raise SystemExit(0)

    lat, lon, elevation = args.lat, args.lon, 0.0
    monument_name = None

    if args.monument:
        # Find monument by case-insensitive substring match
        key = args.monument.lower()
        matches_found = [k for k in MONUMENTS if key in k.lower()]
        if not matches_found:
            print(f"Monument '{args.monument}' not found. Use --list-monuments to see options.")
            raise SystemExit(1)
        monument_name = matches_found[0]
        lat = MONUMENTS[monument_name]['lat']
        lon = MONUMENTS[monument_name]['lon']
        elevation = MONUMENTS[monument_name]['elevation_m']
        print(f"\n[Using monument: {monument_name}  lat={lat}, lon={lon}, elevation={elevation} m]")

    results = calculate_alignments(lat, lon, args.year, args.month, args.day, args.hour, elevation)
    results.update({
        '_lat': lat, '_lon': lon,
        '_year': args.year, '_month': args.month,
        '_day': args.day,   '_hour': args.hour,
    })

    _print_star_table(results)

    if monument_name:
        orientation_az = MONUMENTS[monument_name]['orientation_az']
        alignment_matches = check_alignments(results, orientation_az, args.threshold)
        _print_alignment_results(alignment_matches, monument_name, args.threshold)
