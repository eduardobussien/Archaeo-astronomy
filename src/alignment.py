import argparse
import warnings
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from astropy.utils.exceptions import AstropyWarning

warnings.simplefilter('ignore', category=AstropyWarning)

# High-Precision Star Data (Epoch J2000)
# dist: approximate distance in parsecs
# pm_ra: proper motion in Right Ascension (mas/yr)
# pm_dec: proper motion in Declination (mas/yr)
STARS = {
    'Sirius':     {'ra': 101.287, 'dec': -16.716, 'dist': 2.64, 'pm_ra': -546.01, 'pm_dec': -1223.07},
    'Betelgeuse': {'ra': 88.793,  'dec': 7.407,   'dist': 168,  'pm_ra': 27.54,   'pm_dec': 11.30},
    'Rigel':      {'ra': 78.634,  'dec': -8.202,  'dist': 264,  'pm_ra': 1.87,    'pm_dec': -0.56},
    'Alnitak':    {'ra': 85.190,  'dec': -1.942,  'dist': 226,  'pm_ra': 3.99,    'pm_dec': 2.54},
    'Alnilam':    {'ra': 84.053,  'dec': -1.202,  'dist': 606,  'pm_ra': 1.44,    'pm_dec': -0.73},
    'Mintaka':    {'ra': 83.002,  'dec': -0.299,  'dist': 280,  'pm_ra': 0.18,    'pm_dec': -0.58},
}

def get_observation_time(year, month, day, hour):
    """Safely construct a historical time object, handling BC years."""
    h = int(hour)
    m = int((hour - h) * 60)
    s = ((hour - h) * 60 - m) * 60

    sign = "-" if year < 0 else ""
    time_str = f"{sign}{abs(year):04d}-{month:02d}-{day:02d} {h:02d}:{m:02d}:{s:05.2f}"
    return Time(time_str)

def calculate_alignments(lat, lon, year, month, day, hour):
    obs_time = get_observation_time(year, month, day, hour)
    loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg)
    
    altaz_frame = AltAz(obstime=obs_time, location=loc)
    lst = obs_time.sidereal_time('mean', longitude=loc.lon).degree
    
    results = {'jd': obs_time.jd, 'lst': lst, 'stars': {}}
    
    for name, data in STARS.items():
        star = SkyCoord(
            ra=data['ra'] * u.deg, 
            dec=data['dec'] * u.deg,
            distance=data['dist'] * u.pc,
            pm_ra_cosdec=data['pm_ra'] * u.mas/u.yr,
            pm_dec=data['pm_dec'] * u.mas/u.yr,
            obstime=Time('J2000')
        )
        
        star_historical = star.apply_space_motion(new_obstime=obs_time)
        
        altaz = star_historical.transform_to(altaz_frame)
        
        results['stars'][name] = {
            'altitude': altaz.alt.degree,
            'azimuth': altaz.az.degree,
            'visible': altaz.alt.degree > 0
        }
        
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Archaeo-Astronomy Alignment Tool (High Precision)")
    parser.add_argument("--lat", type=float, default=29.9792, help="Latitude")
    parser.add_argument("--lon", type=float, default=31.1342, help="Longitude")
    parser.add_argument("--year", type=int, default=-10499, help="Historical year (negative for BC)")
    parser.add_argument("--month", type=int, default=3, help="Month")
    parser.add_argument("--day", type=int, default=20, help="Day")
    parser.add_argument("--hour", type=float, default=0.0, help="Hour (Military)")

    args = parser.parse_args()
    r = calculate_alignments(args.lat, args.lon, args.year, args.month, args.day, args.hour)
    
    era = "BC" if args.year < 0 else "AD"
    print(f"\n=== Archaeo-Astronomy Alignment (Astropy Precision) ===")
    print(f"Location: {args.lat}°, {args.lon}°")
    print(f"Date:     {abs(args.year)} {era}-{args.month:02d}-{args.day:02d}")
    print(f"Time:     {args.hour:05.2f} (Military)")
    print(f"Metrics:  JD: {r['jd']:.2f} | LST: {r['lst']:.2f}°")
    print("-" * 55)
    print(f"{'STAR':<12s} | {'ALTITUDE':<10s} | {'AZIMUTH':<10s} | VISIBILITY")
    print("-" * 55)
    
    for name, data in r['stars'].items():
        vis = "VISIBLE" if data['visible'] else "HIDDEN"
        print(f"{name:12s} | {data['altitude']:7.2f}°   | {data['azimuth']:7.2f}°   | {vis}")
    print("\n")