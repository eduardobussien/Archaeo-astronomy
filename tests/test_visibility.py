import numpy as np
import pytest

from alignment import _date_to_jd, find_heliacal_rising
from data import STARS
from visibility import DEFAULT_EXTINCTION, visibility_thresholds

GIZA_LAT, GIZA_LON, GIZA_ELEVATION = 29.9792, 31.1342, 60
SIRIUS_MAG = STARS['Sirius']['mag']


def test_sirius_thresholds_match_egyptological_criteria():
    # Schaefer (2000): Sirius about 6 deg up with the Sun about 5 deg down, some 10 deg
    # apart (9 to 11 with atmospheric transparency). IMCCE: Sun at -7 deg, Sirius at 2-3 deg.
    star, sun = visibility_thresholds(SIRIUS_MAG)
    assert 9.0 <= star - sun <= 11.0
    assert -7.0 <= sun <= -5.0
    assert 2.0 <= star <= 6.0


def test_fainter_stars_need_a_darker_sky_and_more_altitude():
    magnitudes = np.linspace(-1.5, 3.7, 50)
    star, sun = visibility_thresholds(magnitudes)
    assert np.all(np.diff(sun) < 0)
    assert np.all(np.diff(star) > 0)


def test_haze_raises_the_altitude_a_star_needs():
    clear, _ = visibility_thresholds(1.0, 0.20)
    hazy, _ = visibility_thresholds(1.0, 0.35)
    assert hazy > clear


def test_haze_delays_the_heliacal_rising_of_sirius():
    def rising_jd(extinction):
        r = find_heliacal_rising(GIZA_LAT, GIZA_LON, -2780, 'Sirius', extinction, GIZA_ELEVATION)
        return _date_to_jd(r['year'], r['month'], r['day'], 0.0)
    assert rising_jd(0.20) < rising_jd(DEFAULT_EXTINCTION) < rising_jd(0.35)


def test_pole_star_of_2781_bc_has_no_heliacal_rising():
    result = find_heliacal_rising(GIZA_LAT, GIZA_LON, -2780, 'Thuban')
    assert result == {'found': False, 'reason': 'always_visible'}


def test_canopus_is_lost_in_haze_from_giza():
    result = find_heliacal_rising(GIZA_LAT, GIZA_LON, -2780, 'Canopus', 0.35, GIZA_ELEVATION)
    assert result == {'found': False, 'reason': 'never_visible'}
