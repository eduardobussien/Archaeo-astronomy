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

# Sunrise time
sun_position = get_sun(target_date)

# AltAz coordinate frame
altaz_frame = AltAz(obstime=target_date, location=pyramid_location)
sun_altaz = sun_position.transform_to(altaz_frame)

print(f"\nAt 6:00 AM local time:")
print(f"Sun altitude: {sun_altaz.alt:.2f}")
print(f"Sun azimuth: {sun_altaz.az:.2f}")

print("\nFinding exact sunrise time...")

for minutes_before in range(0, 180, 15):  
    total_minutes = 360 - minutes_before  
    hours = total_minutes // 60
    minutes = total_minutes % 60
    
    test_time = Time({
        'year': TARGET_YEAR,
        'month': TARGET_MONTH,
        'day': TARGET_DAY,
        'hour': int(hours),
        'minute': int(minutes),
        'second': 0
    }, format='ymdhms', scale='ut1')
    
    sun_pos = get_sun(test_time)
    altaz = AltAz(obstime=test_time, location=pyramid_location)
    sun_alt = sun_pos.transform_to(altaz).alt.degree
    
    # Print when cross horizon
    if sun_alt < 0.5 and sun_alt > -0.5:
        print(f"Sunrise at approximately: {hours:02d}:{minutes:02d}")
        print(f"Sun altitude at this time: {sun_alt:.2f}°")
        break