"""
Archaeo-Astronomy Alignment Tool. By EduardoBussien. 02/17/2026
"""

from astropy.coordinates import EarthLocation, SkyCoord, ICRS
from astropy.time import Time
import astropy.units as u
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# Great Pyramid of Giza!!
PYRAMID_LAT = 29.9792 * u.deg
PYRAMID_LON = 31.1342 * u.deg
PYRAMID_HEIGHT = 0 * u.m
pyramid_location = EarthLocation(lat=PYRAMID_LAT, lon=PYRAMID_LON, height=PYRAMID_HEIGHT)


TARGET_YEAR = -10499
# To calculate Julian Date
days_per_year = 365.25
years_from_j2000 = TARGET_YEAR - 2000
jd_offset = years_from_j2000 * days_per_year
j2000_jd = 2451545.0
spring_equinox_offset = 79 + (0/24) 
target_jd = j2000_jd + jd_offset + spring_equinox_offset

print(f"Target observation: Spring Equinox, 10,500 BC (Midnight)")
print(f"Julian Date: {target_jd:.2f}")

# Calculate Local Sidereal Time manually (Angle of Earth's rotation at this moment)
days_since_j2000 = target_jd - j2000_jd
centuries_since_j2000 = days_since_j2000 / 36525.0

# Greenwich Mean Sidereal Time at 0h UT
gmst = 280.46061837 + 360.98564736629 * days_since_j2000 + 0.000387933 * centuries_since_j2000**2
gmst = gmst % 360  # Keep in 0-360 range

# Local Sidereal Time = GMST + longitude
lst = (gmst + PYRAMID_LON.value) % 360

print(f"Local Sidereal Time: {lst:.2f}°")

print("\n" + "="*60)
print("STELLAR ALIGNMENTS - MARCH 20, 10,500 BC")
print("="*60)
print(f"Location: Great Pyramid ({PYRAMID_LAT.value:.4f}°N, {PYRAMID_LON.value:.4f}°E)\n")

# Define stars
famous_stars = {
    'Sirius (Brightest Star)': {'ra': 101.287, 'dec': -16.716},
    'Betelgeuse (Orion)': {'ra': 88.793, 'dec': 7.407},
    'Rigel (Orion)': {'ra': 78.634, 'dec': -8.202},
    'Aldebaran (Taurus)': {'ra': 68.980, 'dec': 16.509},
    'Procyon': {'ra': 114.825, 'dec': 5.225},
    'Alnitak (Belt Left)': {'ra': 85.190, 'dec': -1.942},
    'Alnilam (Belt Center)': {'ra': 84.053, 'dec': -1.202},
    'Mintaka (Belt Right)': {'ra': 83.002, 'dec': -0.299},
}

# Calculate altitude and azimuth for each star
for star_name, coords in famous_stars.items():
    ra = coords['ra']
    dec = coords['dec']
    
    # Hour Angle = LST - RA
    ha = (lst - ra)
    if ha < 0:
        ha += 360
    if ha > 180:
        ha -= 360
    
    # Convert to radians for calculation
    ha_rad = np.radians(ha)
    dec_rad = np.radians(dec)
    lat_rad = np.radians(PYRAMID_LAT.value)
    
    # Calculate altitude
    sin_alt = np.sin(dec_rad) * np.sin(lat_rad) + np.cos(dec_rad) * np.cos(lat_rad) * np.cos(ha_rad)
    altitude = np.degrees(np.arcsin(sin_alt))
    
    # Calculate azimuth
    cos_az = (np.sin(dec_rad) - np.sin(lat_rad) * sin_alt) / (np.cos(lat_rad) * np.cos(np.arcsin(sin_alt)))
    cos_az = np.clip(cos_az, -1, 1)  # Prevent numerical errors
    azimuth = np.degrees(np.arccos(cos_az))
    
    if np.sin(ha_rad) > 0:
        azimuth = 360 - azimuth
    
    # Visibility
    if altitude > 0:
        visibility = "+VISIBLE"
    else:
        visibility = "-Below horizon"
    
    print(f"{star_name:30s} | Alt: {altitude:7.2f}° | Az: {azimuth:6.2f}° | {visibility}")