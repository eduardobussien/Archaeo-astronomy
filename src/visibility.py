"""
When a star first becomes visible in morning twilight.

Reijs's fits to Schaefer's (1993) visibility model (archaeocosmology.org,
"Extinction angle and heliacal events"). At first sighting a star of visual
magnitude m stands at apparent altitude

    h_star = 0.1204 m^2 + 0.7941 m + 5.59063 + (23.276 k - 7.1958)

while the Sun is at geometric altitude

    h_sun = -0.1504 m^2 - 1.646 m - 8.1015 + (0.602 - 2.1165 k)

where k is the total atmospheric extinction in magnitudes per airmass. Reijs
quotes one-sigma uncertainties of 1.4 deg (Sun) and 1.9 deg (star), which
amount to a few days in the date of a heliacal rising. For Sirius this gives
the Sun 6 deg down and the star 4 deg up (about 10 deg apart), matching
Schaefer (2000) and IMCCE's criteria for Egypt.
"""
import numpy as np

# Schaefer's (2000) estimate for ancient Memphis; about 0.35 suits humid coasts.
DEFAULT_EXTINCTION = 0.27


def visibility_thresholds(magnitude, extinction=DEFAULT_EXTINCTION):
    """(minimum apparent star altitude, maximum Sun altitude) in degrees at first sighting."""
    m = np.asarray(magnitude, dtype=float)
    star = 0.1204 * m**2 + 0.7941 * m + 5.59063 + (23.276 * extinction - 7.1958)
    sun = -0.1504 * m**2 - 1.646 * m - 8.1015 + (0.602 - 2.1165 * extinction)
    return star, sun
