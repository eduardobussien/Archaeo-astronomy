"""
Delta T (TT - UT1), the accumulated slowing of Earth's rotation.

Polynomials from Espenak & Meeus (2006), Five Millennium Canon of Solar
Eclipses, which follow Morrison & Stephenson (2004) and reduce to their
long-term parabola -20 + 32 u^2 outside 500 BC to 2150 AD.
"""
import erfa
import numpy as np

# (first year, last year, reference year, scale, coefficients in (year - ref) / scale, constant first)
_SEGMENTS = (
    (-np.inf, -500, 1820, 100, (-20.0, 0.0, 32.0)),
    (-500, 500, 0, 100, (10583.6, -1014.41, 33.78311, -5.952053, -0.1798452,
                         0.022174192, 0.0090316521)),
    (500, 1600, 1000, 100, (1574.2, -556.01, 71.23472, 0.319781, -0.8503463,
                            -0.005050998, 0.0083572073)),
    (1600, 1700, 1600, 1, (120.0, -0.9808, -0.01532, 1 / 7129)),
    (1700, 1800, 1700, 1, (8.83, 0.1603, -0.0059285, 0.00013336, -1 / 1174000)),
    (1800, 1860, 1800, 1, (13.72, -0.332447, 0.0068612, 0.0041116, -0.00037436,
                           0.0000121272, -0.0000001699, 0.000000000875)),
    (1860, 1900, 1860, 1, (7.62, 0.5737, -0.251754, 0.01680668, -0.0004473624, 1 / 233174)),
    (1900, 1920, 1900, 1, (-2.79, 1.494119, -0.0598939, 0.0061966, -0.000197)),
    (1920, 1941, 1920, 1, (21.20, 0.84493, -0.076100, 0.0020936)),
    (1941, 1961, 1950, 1, (29.07, 0.407, -1 / 233, 1 / 2547)),
    (1961, 1986, 1975, 1, (45.45, 1.067, -1 / 260, -1 / 718)),
    (1986, 2005, 2000, 1, (63.86, 0.3345, -0.060374, 0.0017275, 0.000651814, 0.00002373599)),
    (2005, 2050, 2000, 1, (62.92, 0.32217, 0.005589)),
    (2050, 2150, 1820, 100, (-205.724, 56.28, 32.0)),
    (2150, np.inf, 1820, 100, (-20.0, 0.0, 32.0)),
)


def delta_t(year):
    """Delta T in seconds for a decimal year (astronomical numbering)."""
    year = np.asarray(year, dtype=float)
    result = np.empty_like(year)
    for first, last, ref, scale, coeffs in _SEGMENTS:
        mask = (year >= first) & (year < last)
        result[mask] = np.polynomial.polynomial.polyval((year[mask] - ref) / scale, coeffs)
    return result if result.ndim else float(result)


def ut1_to_tt(jd_ut1):
    """Terrestrial Time Julian Date for a UT1 Julian Date."""
    return jd_ut1 + delta_t(erfa.epj(jd_ut1, 0.0)) / 86400.0
