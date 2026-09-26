import argparse
import math
import warnings
import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_body, solar_system_ephemeris
import astropy.units as u
from astropy.utils.exceptions import AstropyWarning
import erfa

from data import STARS, MONUMENTS
from precession import star_altaz, ecliptic_altaz, local_mean_sidereal_time

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=erfa.ErfaWarning)

_STAR_NAMES = list(STARS)
_CATALOG = {
    key: np.array([STARS[name][key] for name in _STAR_NAMES])
    for key in ('ra', 'dec', 'pm_ra', 'pm_dec', 'dist')
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


def format_year(year):
    """Label an astronomical year (0 = 1 BC, -1 = 2 BC) in historical form."""
    return f"{1 - year} BC" if year <= 0 else f"{year} AD"


def get_observation_time(year, month, day, hour):
    """
    Build an astropy Time for any historical date.

    Uses Julian Day Number (avoids the ISO parser, which rejects BC years).
    UT1 scale is used so astropy treats the JD directly as Earth rotation
    time; delta_ut1_utc is pre-set to 0 so no IERS lookup is attempted.
    """
    jd = _date_to_jd(year, month, day, hour)
    t = Time(jd, format='jd', scale='ut1')
    t.delta_ut1_utc = 0.0   # bypass the IERS table; ΔUT1=0 is the best
    return t                 # estimate for ancient dates anyway


def _gmst_degrees(jd):
    """
    Greenwich Mean Sidereal Time in degrees via the Meeus formula.
    Valid for any Julian Date (no ERFA dependency).
    """
    D = jd - 2451545.0
    T = D / 36525.0
    gmst = (280.46061837 + 360.98564736629 * D
            + 0.000387933 * T**2 - T**3 / 38710000.0)
    return gmst % 360.0


_PLANET_NAMES = ['sun', 'moon', 'mercury', 'venus', 'mars', 'jupiter', 'saturn']
solar_system_ephemeris.set('builtin')   # set once at import; safe across threads


def _calculate_planets(obs_time, loc, altaz_frame):
    """Return altitude/azimuth for the Sun, Moon, and five naked-eye planets."""
    planets = {}
    for name in _PLANET_NAMES:
        try:
            body  = get_body(name, obs_time, location=loc)
            altaz = body.transform_to(altaz_frame)
            planets[name.capitalize()] = {
                'altitude': float(altaz.alt.degree),
                'azimuth':  float(altaz.az.degree),
                'visible':  bool(altaz.alt.degree > 0),
            }
        except Exception:
            pass  # skip bodies outside the builtin ephemeris range
    return planets


def calculate_alignments(lat, lon, year, month, day, hour):
    """
    Compute altitude and azimuth for every star in the catalog.

    Stars use the long-term precession engine in precession.py at every epoch.
    Planets use astropy's builtin ephemeris and are omitted before ~4800 BC,
    where ERFA's calendar routines refuse the date.

    Returns a dict with keys:
        jd      - Julian Date (UT1)
        lst     - Local mean sidereal time in degrees
        stars   - {name: {altitude, azimuth, visible}} for each star
        planets - {name: {altitude, azimuth, visible}}, possibly empty
        method  - name of the star position model
    """
    # Convert local mean solar time to UT (the user picks local time; lon/15 is the offset)
    hour_ut = hour - lon / 15.0
    jd = _date_to_jd(year, month, day, hour_ut)

    alt, az = star_altaz(_CATALOG['ra'], _CATALOG['dec'], _CATALOG['pm_ra'],
                         _CATALOG['pm_dec'], _CATALOG['dist'], jd, lat, lon)
    stars_out = {
        name: {'altitude': float(a), 'azimuth': float(z), 'visible': bool(a > 0)}
        for name, a, z in zip(_STAR_NAMES, alt, az)
    }

    try:
        obs_time = get_observation_time(year, month, day, hour_ut)
        loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
        planets = _calculate_planets(obs_time, loc, AltAz(obstime=obs_time, location=loc))
    except erfa.ErfaError:
        planets = {}

    return {
        'jd': jd,
        'lst': float(local_mean_sidereal_time(jd, lon)),
        'stars': stars_out,
        'planets': planets,
        'method': 'Vondrák 2011',
    }


def _approx_sun_altitude(jd, lat, lon):
    """
    Low-precision Sun altitude (±1° within ±3000 yr of J2000).
    Uses Meeus Astronomical Algorithms ch. 25.  jd must be UT-based.
    """
    T = (jd - 2451545.0) / 36525.0
    L0 = (280.46646 + 36000.76983 * T + 0.0003032 * T**2) % 360
    M  = math.radians((357.52911 + 35999.05029 * T - 0.0001537 * T**2) % 360)
    C  = ((1.914602 - 0.004817 * T - 0.000014 * T**2) * math.sin(M)
          + (0.019993 - 0.000101 * T) * math.sin(2 * M)
          + 0.000289 * math.sin(3 * M))
    sun_lon = math.radians((L0 + C) % 360)
    eps     = math.radians(23.439291 - 0.013004 * T)
    sun_ra  = math.atan2(math.cos(eps) * math.sin(sun_lon), math.cos(sun_lon))
    sun_dec = math.asin(max(-1.0, min(1.0, math.sin(eps) * math.sin(sun_lon))))
    lst     = math.radians((_gmst_degrees(jd) + lon) % 360)
    ha      = lst - sun_ra
    sin_alt = (math.sin(math.radians(lat)) * math.sin(sun_dec)
               + math.cos(math.radians(lat)) * math.cos(sun_dec) * math.cos(ha))
    return math.degrees(math.asin(max(-1.0, min(1.0, sin_alt))))


def _find_dawn_hour(lat, lon, year, month, day, arc_vision=-10.0):
    """
    Binary search for the local hour (0–14) when the Sun crosses arc_vision
    on the way up (morning twilight).  Returns None for polar day/night.
    """
    def sun_alt(h):
        return _approx_sun_altitude(_date_to_jd(year, month, day, h - lon / 15.0), lat, lon)

    alt_lo, alt_hi = sun_alt(0.0), sun_alt(14.0)
    lo, hi = 0.0, 14.0
    # Need arc_vision to be strictly between the two endpoints
    if not ((alt_lo < arc_vision < alt_hi) or (alt_hi < arc_vision < alt_lo)):
        return None
    for _ in range(20):
        mid = (lo + hi) / 2.0
        alt_mid = sun_alt(mid)
        if (alt_lo - arc_vision) * (alt_mid - arc_vision) <= 0:
            hi, alt_hi = mid, alt_mid
        else:
            lo, alt_lo = mid, alt_mid
    return (lo + hi) / 2.0


def find_heliacal_rising(lat, lon, year, star_name, arc_vision=-10.0):
    """
    Scan day-by-day to find when star_name first becomes visible at
    astronomical twilight (Sun altitude = arc_vision) in the target year.

    Starts scanning from Oct 1 of (year-1) so the initial visibility state
    is correctly seeded before Jan 1 of the target year.

    Returns a dict with the date and geometry, or None if not found.
    """
    if star_name not in STARS:
        return None

    star_info = STARS[star_name]
    DAYS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    scan_year, scan_month, scan_day = year - 1, 10, 1
    prev_visible = None

    for _ in range(460):   # 92 days warm-up + 366 target year + buffer
        dawn_hour = _find_dawn_hour(lat, lon, scan_year, scan_month, scan_day, arc_vision)
        star_alt = star_az = 0.0
        visible = False

        if dawn_hour is not None:
            jd = _date_to_jd(scan_year, scan_month, scan_day, dawn_hour - lon / 15.0)
            star_alt, star_az = (float(v) for v in star_altaz(
                star_info['ra'], star_info['dec'], star_info['pm_ra'],
                star_info['pm_dec'], star_info['dist'], jd, lat, lon))
            visible = star_alt > 0.5

        if scan_year == year and visible and prev_visible is False:
            sun_alt = _approx_sun_altitude(
                _date_to_jd(scan_year, scan_month, scan_day, dawn_hour - lon / 15.0),
                lat, lon,
            )
            return {
                'year':  scan_year,
                'month': scan_month,
                'day':   scan_day,
                'star_altitude':   round(star_alt, 2),
                'star_azimuth':    round(star_az, 2),
                'sun_altitude':    round(sun_alt, 2),
                'dawn_hour_local': round(dawn_hour, 2),
            }

        prev_visible = visible

        # Advance one day
        scan_day += 1
        if scan_day > DAYS[scan_month - 1]:
            scan_day = 1
            scan_month += 1
            if scan_month > 12:
                scan_month = 1
                scan_year += 1
        if scan_year > year:
            break

    return None


def calculate_ecliptic(lat, lon, year, month, day, hour):
    """
    Return 73 points (0°..360° ecliptic longitude, step 5°) projected onto
    the local alt-az frame.  The 73rd point closes the loop back to 0°.

    Uses the mean ecliptic of date from the long-term precession model.
    """
    jd = _date_to_jd(year, month, day, hour - lon / 15.0)
    longitudes = np.arange(73) * 5.0
    alt, az = ecliptic_altaz(jd, lat, lon, longitudes)
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
    engine = results.get('method', 'astropy')
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

    lat, lon = args.lat, args.lon
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
        print(f"\n[Using monument: {monument_name}  lat={lat}, lon={lon}]")

    results = calculate_alignments(lat, lon, args.year, args.month, args.day, args.hour)
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
