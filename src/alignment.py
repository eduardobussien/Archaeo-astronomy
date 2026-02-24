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

# Sunrise time ! 
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


print("\n" + "="*60)
print("STELLAR ALIGNMENTS AT SUNRISE")
print("="*60)

sunrise_time = Time({
    'year': TARGET_YEAR,
    'month': TARGET_MONTH,
    'day': TARGET_DAY,
    'hour': 16,
    'minute': 0,
    'second': 0
}, format='ymdhms', scale='ut1')

# Define stars and constellations (just trying out for now)
famous_stars = {
    'Sirius (Brightest Star)': SkyCoord(ra='06h45m08.9s', dec='-16d42m58s', frame='icrs'),
    'Betelgeuse (Orion)': SkyCoord(ra='05h55m10.3s', dec='+07d24m25s', frame='icrs'),
    'Rigel (Orion)': SkyCoord(ra='05h14m32.3s', dec='-08d12m06s', frame='icrs'),
    'Aldebaran (Taurus)': SkyCoord(ra='04h35m55.2s', dec='+16d30m33s', frame='icrs'),
    'Procyon': SkyCoord(ra='07h39m18.1s', dec='+05d13m30s', frame='icrs'),
}

# Create AltAz frame for sunrise time
sunrise_altaz = AltAz(obstime=sunrise_time, location=pyramid_location)

print(f"\nAt sunrise ({sunrise_time.iso}):")
print(f"Location: Great Pyramid ({PYRAMID_LAT.value:.4f}°N, {PYRAMID_LON.value:.4f}°E)\n")

for star_name, star_coord in famous_stars.items():
    star_altaz = star_coord.transform_to(sunrise_altaz)
    
    altitude = star_altaz.alt.degree
    azimuth = star_altaz.az.degree
    
    # Determine if star is visible
    if altitude > 0:
        visibility = "+VISIBLE"
    else:
        visibility = "-Below horizon"
    
    print(f"{star_name:30s} | Alt: {altitude:7.2f}° | Az: {azimuth:6.2f}° | {visibility}")