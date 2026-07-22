import argparse
import math
import warnings
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body, solar_system_ephemeris
import astropy.units as u
from astropy.utils.exceptions import AstropyWarning
import erfa

from data import STARS, MONUMENTS

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=erfa.ErfaWarning)


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


# ---------------------------------------------------------------------------
# Manual calculation fallback (for dates outside ERFA's ~4800 BC limit)
# ---------------------------------------------------------------------------

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


def _apply_precession_lieske(ra_deg, dec_deg, T):
    """
    Apply IAU 1976 precession (Lieske 1979) to J2000 equatorial coordinates.
    T = Julian centuries from J2000.0.  Valid to roughly ±10 000 years.
    """
    zeta  = (0.6406161 + 0.0000839 * T + 0.0000050 * T**2) * T
    z     = (0.6406161 + 0.0003041 * T + 0.0000051 * T**2) * T
    theta = (0.5567530 - 0.0001185 * T - 0.0000116 * T**2) * T

    zeta_r  = math.radians(zeta)
    z_r     = math.radians(z)
    theta_r = math.radians(theta)
    ra_r    = math.radians(ra_deg)
    dec_r   = math.radians(dec_deg)

    A = math.cos(dec_r) * math.cos(ra_r + zeta_r)
    B = (math.cos(theta_r) * math.cos(dec_r) * math.sin(ra_r + zeta_r)
         - math.sin(theta_r) * math.sin(dec_r))
    C = (math.sin(theta_r) * math.cos(dec_r) * math.sin(ra_r + zeta_r)
         + math.cos(theta_r) * math.sin(dec_r))

    ra_prec  = (math.degrees(math.atan2(B, A)) + z) % 360.0
    dec_prec = math.degrees(math.asin(max(-1.0, min(1.0, C))))
    return ra_prec, dec_prec


def _ha_dec_to_altaz(ha_deg, dec_deg, lat_deg):
    """Convert Hour Angle and Declination to (altitude, azimuth) in degrees."""
    ha  = math.radians(ha_deg)
    dec = math.radians(dec_deg)
    lat = math.radians(lat_deg)

    sin_alt = math.sin(lat) * math.sin(dec) + math.cos(lat) * math.cos(dec) * math.cos(ha)
    sin_alt = max(-1.0, min(1.0, sin_alt))
    alt = math.degrees(math.asin(sin_alt))

    cos_alt = math.cos(math.radians(alt))
    if cos_alt < 1e-10:
        return alt, 0.0

    cos_az = (math.sin(dec) - math.sin(lat) * sin_alt) / (math.cos(lat) * cos_alt)
    cos_az = max(-1.0, min(1.0, cos_az))
    az = math.degrees(math.acos(cos_az))
    if math.sin(ha) >= 0:
        az = 360.0 - az

    return alt, az


def _calculate_manual(lat, lon, jd):
    """
    Compute star positions using polynomial GMST + Lieske precession.
    Used automatically for dates before ~4800 BC where ERFA errors out.
    """
    T   = (jd - 2451545.0) / 36525.0
    dt_years = T * 100.0
    lst = (_gmst_degrees(jd) + lon) % 360.0

    stars_out = {}
    for name, data in STARS.items():
        # Proper motion: linear angular displacement
        ra  = data['ra']  + (data['pm_ra']  / 3.6e6 / math.cos(math.radians(data['dec']))) * dt_years
        dec = data['dec'] + (data['pm_dec'] / 3.6e6) * dt_years

        ra_prec, dec_prec = _apply_precession_lieske(ra, dec, T)

        ha = (lst - ra_prec) % 360.0
        alt, az = _ha_dec_to_altaz(ha, dec_prec, lat)

        stars_out[name] = {'altitude': alt, 'azimuth': az, 'visible': alt > 0}

    return stars_out, lst


_PLANET_NAMES = ['sun', 'moon', 'mercury', 'venus', 'mars', 'jupiter', 'saturn']

def _calculate_planets(obs_time, loc, altaz_frame):
    """Return altitude/azimuth for the Sun, Moon, and five naked-eye planets."""
    planets = {}
    try:
        with solar_system_ephemeris.set('builtin'):
            for name in _PLANET_NAMES:
                body  = get_body(name, obs_time, location=loc)
                altaz = body.transform_to(altaz_frame)
                planets[name.capitalize()] = {
                    'altitude': float(altaz.alt.degree),
                    'azimuth':  float(altaz.az.degree),
                    'visible':  bool(altaz.alt.degree > 0),
                }
    except Exception:
        pass  # silently skip if date is outside the builtin ephemeris range
    return planets


def calculate_alignments(lat, lon, year, month, day, hour):
    """
    Compute altitude and azimuth for every star in the catalog.

    Uses astropy/ERFA for dates after ~4800 BC (full precession + proper motion).
    Falls back to Meeus polynomial formulas for deeper historical dates.

    Returns a dict with keys:
        jd      - Julian Date
        lst     - Local Sidereal Time in degrees
        stars   - {name: {altitude, azimuth, visible}} for each star
        method  - 'astropy' or 'manual' (indicates which engine was used)
    """
    jd = _date_to_jd(year, month, day, hour)

    try:
        obs_time  = get_observation_time(year, month, day, hour)
        loc       = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
        altaz_frame = AltAz(obstime=obs_time, location=loc)
        lst = obs_time.sidereal_time('mean', longitude=loc.lon).degree

        stars_out = {}
        for name, data in STARS.items():
            star = SkyCoord(
                ra=data['ra'] * u.deg,
                dec=data['dec'] * u.deg,
                distance=data['dist'] * u.pc,
                pm_ra_cosdec=data['pm_ra'] * u.mas / u.yr,
                pm_dec=data['pm_dec'] * u.mas / u.yr,
                obstime=Time('J2000'),
            )
            star_at_epoch = star.apply_space_motion(new_obstime=obs_time)
            altaz = star_at_epoch.transform_to(altaz_frame)
            stars_out[name] = {
                'altitude': altaz.alt.degree,
                'azimuth':  altaz.az.degree,
                'visible':  altaz.alt.degree > 0,
            }

        planets = _calculate_planets(obs_time, loc, altaz_frame)
        return {'jd': jd, 'lst': lst, 'stars': stars_out, 'planets': planets, 'method': 'astropy'}

    except erfa.ErfaError:
        stars_out, lst = _calculate_manual(lat, lon, jd)
        return {'jd': jd, 'lst': lst, 'stars': stars_out, 'planets': {}, 'method': 'manual'}


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
    era = "BC" if results['_year'] < 0 else "AD"
    engine = results.get('method', 'astropy')
    print(f"\n=== Archaeo-Astronomy Alignment [{engine}] ===")
    print(f"Location : {results['_lat']:.4f}, {results['_lon']:.4f}")
    print(f"Date     : {abs(results['_year'])} {era}-{results['_month']:02d}-{results['_day']:02d}")
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
        epilog="Example:\n  python src/alignment.py --monument Stonehenge --year -2500\n"
               "  python src/alignment.py --list-monuments",
    )
    parser.add_argument("--lat",    type=float, default=29.9792)
    parser.add_argument("--lon",    type=float, default=31.1342)
    parser.add_argument("--year",   type=int,   default=-2500,
                        help="Historical year (negative = BC)")
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
