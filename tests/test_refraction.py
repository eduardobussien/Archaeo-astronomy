import math

import erfa
import numpy as np
import pytest

from refraction import apparent_altitude, pressure_ratio, refraction


def test_horizontal_refraction_is_about_34_arcmin():
    true_altitudes = np.linspace(-1.0, 0.0, 100001)
    at_horizon = true_altitudes[np.argmin(np.abs(apparent_altitude(true_altitudes)))]
    assert abs(-at_horizon * 60.0 - 34.4) < 0.3


@pytest.mark.parametrize('altitude', [10, 15, 20, 30, 45, 60, 80])
def test_matches_erfa_refraction_above_10_degrees(altitude):
    a, b = erfa.refco(1013.25, 10.0, 0.5, 0.55)
    tan_z = math.tan(math.radians(90.0 - altitude))
    erfa_arcmin = math.degrees(a * tan_z + b * tan_z ** 3) * 60.0
    assert abs(refraction(altitude) * 60.0 - erfa_arcmin) < 0.15


def test_vanishes_at_zenith():
    assert refraction(90.0) * 3600.0 < 0.1


def test_apparent_altitude_preserves_order():
    true_altitudes = np.linspace(-10.0, 90.0, 20001)
    assert np.all(np.diff(apparent_altitude(true_altitudes)) > 0)


@pytest.mark.parametrize('elevation_m, ratio', [(0, 1.0), (2300, 0.756), (3850, 0.620)])
def test_pressure_follows_standard_atmosphere(elevation_m, ratio):
    assert abs(pressure_ratio(elevation_m) - ratio) < 0.002
