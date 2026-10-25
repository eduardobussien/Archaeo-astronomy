import numpy as np
import pytest

from horizon import MAX_POINTS, HorizonProfile


def test_interpolates_linearly_between_points():
    p = HorizonProfile([0, 90, 180, 270], [1.0, 3.0, 2.0, 0.0])
    assert p(45) == pytest.approx(2.0)
    assert p(135) == pytest.approx(2.5)
    assert np.allclose(p([90, 180]), [3.0, 2.0])


def test_wraps_through_north():
    p = HorizonProfile([10, 350], [4.0, 2.0])
    assert p(0) == pytest.approx(3.0)
    assert p(360) == pytest.approx(3.0)
    assert p(355) == pytest.approx(2.5)
    assert p(-5) == pytest.approx(2.5)              # the same direction as 355
    assert p(180) == pytest.approx(3.0)             # halfway round the long way


def test_a_single_point_is_a_level_horizon():
    p = HorizonProfile([123], [1.5])
    assert np.allclose(p([0, 90, 359.9]), 1.5)


def test_reads_survey_and_stellarium_style_text():
    text = """azimuth,altitude,distance_km
    # surveyed 2026-05-01
    0, 1.0, 12
    90\t3.0\t4
    180 2.0      # behind the ridge
    270;0.5
    """
    p = HorizonProfile.from_text(text)
    assert list(p.azimuths) == [0, 90, 180, 270]
    assert list(p.altitudes) == [1.0, 3.0, 2.0, 0.5]


def test_compact_form_round_trips():
    p = HorizonProfile([270, 0, 90.5], [0.5, 1.25, 3.0])
    assert p.compact() == '0:1.25,90.5:3,270:0.5'
    q = HorizonProfile.from_compact(p.compact())
    assert np.array_equal(q.azimuths, p.azimuths) and np.array_equal(q.altitudes, p.altitudes)


@pytest.mark.parametrize('text, message', [
    ('', 'no points'),
    ('# only a comment', 'no points'),
    ('0 1\n90', 'needs an azimuth and an altitude'),
    ('0 1\nnorth 2', 'Line 2'),
    ('0 1\n90 50', 'between -2 and 45'),
    ('0 1\n400 2', 'between 0 and 360'),
    ('0 1\n360 2', 'only once'),
    ('0 1\nnan 2', 'must be numbers'),
])
def test_rejects_bad_profiles_with_a_clear_reason(text, message):
    with pytest.raises(ValueError, match=message):
        HorizonProfile.from_text(text)


def test_rejects_too_many_points_and_bad_compact_text():
    with pytest.raises(ValueError, match=f'1 to {MAX_POINTS}'):
        HorizonProfile(np.linspace(0, 359.9, MAX_POINTS + 1), np.zeros(MAX_POINTS + 1))
    with pytest.raises(ValueError, match='azimuth:altitude'):
        HorizonProfile.from_compact('0:1,90')
