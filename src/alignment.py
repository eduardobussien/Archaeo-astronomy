"""
Archaeo-Astronomy Alignment Tool. By EduardoBussien. 02/17/2026

"""

from astropy.coordinates import EarthLocation, AltAz, get_sun, SkyCoord
from astropy.time import Time
import astropy.units as u
from datetime import datetime
import numpy as np

# Great Pyramid of Giza
PYRAMID_LAT = 29.9792 * u.deg  # degrees North
PYRAMID_LON = 31.1342 * u.deg  # degrees East
PYRAMID_HEIGHT = 0 * u.m       # meters above sea level 

# EarthLocation object
pyramid_location = EarthLocation(lat=PYRAMID_LAT, lon=PYRAMID_LON, height=PYRAMID_HEIGHT)

# Target date: Spring Equinox 2500 BC
TARGET_YEAR = -2499  
TARGET_MONTH = 3
TARGET_DAY = 20

target_date = Time({
    'year': TARGET_YEAR,
    'month': TARGET_MONTH, 
    'day': TARGET_DAY,
    'hour': 6,
    'minute': 0,
    'second': 0
}, format='ymdhms', scale='ut1')

print(f"Target observation date: {target_date.iso}")
print(f"Julian Date: {target_date.jd}")