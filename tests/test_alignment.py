import pytest

from alignment import _date_to_jd, _jd_to_date, format_year


@pytest.mark.parametrize('year, label', [
    (1, '1 AD'),
    (0, '1 BC'),
    (-1, '2 BC'),
    (-2499, '2500 BC'),
    (-10499, '10500 BC'),
])
def test_format_year_uses_astronomical_numbering(year, label):
    assert format_year(year) == label


@pytest.mark.parametrize('year', [-12000, -10499, -4713, -4712, -2780, -1, 0, 1, 1582, 1900, 2000, 2024])
@pytest.mark.parametrize('month, day', [(1, 1), (2, 28), (2, 29), (3, 1), (12, 31)])
@pytest.mark.parametrize('hour', [0.0, 12.0, 23.99])
def test_jd_to_date_inverts_date_to_jd(year, month, day, hour):
    is_leap = (year % 4 == 0 and year % 100 != 0) or year % 400 == 0
    if (month, day) == (2, 29) and not is_leap:
        pytest.skip('not a leap year')
    assert _jd_to_date(_date_to_jd(year, month, day, hour)) == (year, month, day)
