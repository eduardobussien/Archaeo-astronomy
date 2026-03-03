" Archaeo-Astronomy Alignment Tool - Calculate stellar alignments for any location and historical date. By eduardobussien"

import numpy as np

J2000_JD = 2451545.0
DAYS_IN_MONTH = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]

STARS = {
    'Sirius': {'ra': 101.287, 'dec': -16.716},
    'Betelgeuse': {'ra': 88.793, 'dec': 7.407},
    'Rigel': {'ra': 78.634, 'dec': -8.202},
    'Alnitak': {'ra': 85.190, 'dec': -1.942},
    'Alnilam': {'ra': 84.053, 'dec': -1.202},
    'Mintaka': {'ra': 83.002, 'dec': -0.299},
}


def date_to_jd(year, month, day, hour):
    years_from_j2000 = year - 2000
    days_from_years = years_from_j2000 * 365.25
    day_of_year = DAYS_IN_MONTH[month - 1] + day
    day_fraction = hour / 24.0
    return J2000_JD + days_from_years + day_of_year + day_fraction


def calculate_lst(jd, longitude):
    days_since_j2000 = jd - J2000_JD
    centuries = days_since_j2000 / 36525.0
    gmst = 280.46061837 + 360.98564736629 * days_since_j2000 + 0.000387933 * centuries**2
    gmst = gmst % 360
    return (gmst + longitude) % 360


def star_altaz(ra, dec, lst, latitude):
    ha = lst - ra
    if ha < -180:
        ha += 360
    if ha > 180:
        ha -= 360
    
    ha_rad = np.radians(ha)
    dec_rad = np.radians(dec)
    lat_rad = np.radians(latitude)
    
    sin_alt = np.sin(dec_rad) * np.sin(lat_rad) + np.cos(dec_rad) * np.cos(lat_rad) * np.cos(ha_rad)
    altitude = np.degrees(np.arcsin(sin_alt))
    
    cos_az = (np.sin(dec_rad) - np.sin(lat_rad) * sin_alt) / (np.cos(lat_rad) * np.cos(np.arcsin(sin_alt)))
    cos_az = np.clip(cos_az, -1, 1)
    azimuth = np.degrees(np.arccos(cos_az))
    
    if np.sin(ha_rad) > 0:
        azimuth = 360 - azimuth
    
    return altitude, azimuth


def calculate_alignments(lat, lon, year, month, day, hour):
    jd = date_to_jd(year, month, day, hour)
    lst = calculate_lst(jd, lon)
    
    results = {'jd': jd, 'lst': lst, 'stars': {}}
    
    for name, coords in STARS.items():
        alt, az = star_altaz(coords['ra'], coords['dec'], lst, lat)
        results['stars'][name] = {
            'altitude': alt,
            'azimuth': az,
            'visible': alt > 0
        }
    
    return results


if __name__ == "__main__":
    r = calculate_alignments(29.9792, 31.1342, -10499, 3, 20, 0)
    
    print(f"10,500 BC - Great Pyramid - Midnight")
    print(f"JD: {r['jd']:.2f} | LST: {r['lst']:.2f}°\n")
    
    for name, data in r['stars'].items():
        vis = "VISIBLE" if data['visible'] else "HIDDEN"
        print(f"{name:12s} | Alt: {data['altitude']:6.2f}° | Az: {data['azimuth']:6.2f}° | {vis}")