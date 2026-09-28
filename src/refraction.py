"""
Atmospheric refraction: the lift that makes objects appear higher than they are.

Saemundsson (1986): R = 1.02' / tan(h + 10.3 / (h + 5.11)) for a geometric
altitude h in degrees, at 1010 hPa and 10 C (about 34' at the horizon). It is
scaled by the air pressure at the site's elevation from the International
Standard Atmosphere. Below -1 deg, where an object is already hidden and the
formula diverges, the -1 deg value is held.
"""
import numpy as np

_MIN_ALTITUDE = -1.0


def pressure_ratio(elevation_m):
    """Air pressure at elevation_m relative to sea level (standard atmosphere)."""
    return (1.0 - 2.25577e-5 * np.asarray(elevation_m, dtype=float)) ** 5.25588


def refraction(true_altitude, elevation_m=0.0):
    """Refraction in degrees for geometric altitudes in degrees."""
    h = np.maximum(np.asarray(true_altitude, dtype=float), _MIN_ALTITUDE)
    arcmin = 1.02 / np.tan(np.radians(h + 10.3 / (h + 5.11)))
    return np.maximum(arcmin, 0.0) / 60.0 * pressure_ratio(elevation_m)


def apparent_altitude(true_altitude, elevation_m=0.0):
    """Altitude in degrees as seen through the atmosphere."""
    return np.asarray(true_altitude, dtype=float) + refraction(true_altitude, elevation_m)
