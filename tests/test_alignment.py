import pytest

from alignment import _date_to_jd, _jd_to_date, days_in_month, format_year


@pytest.mark.parametrize('year, label', [
    (1, '1 AD'),
    (0, '1 BC'),
    (-1, '2 BC'),
    (-2499, '2500 BC'),
    (-10499, '10500 BC'),
])
def test_format_year_uses_astronomical_numbering(year, label):
    assert format_year(year) == label


@pytest.mark.parametrize('calendar', ['gregorian', 'julian'])
@pytest.mark.parametrize('year', [-12000, -10499, -4713, -4712, -2780, -1, 0, 1, 1582, 1900, 2000, 2024])
@pytest.mark.parametrize('month, day', [(1, 1), (2, 28), (2, 29), (3, 1), (12, 31)])
@pytest.mark.parametrize('hour', [0.0, 12.0, 23.99])
def test_jd_to_date_inverts_date_to_jd(year, month, day, hour, calendar):
    if (month, day) == (2, 29) and days_in_month(year, 2, calendar) == 28:
        pytest.skip('not a leap year')
    assert _jd_to_date(_date_to_jd(year, month, day, hour, calendar), calendar) == (year, month, day)


def test_julian_calendar_anchors():
    assert _date_to_jd(-4712, 1, 1, 12.0, 'julian') == 0.0             # origin of the Julian Date
    assert _date_to_jd(1582, 10, 4, 0.0, 'julian') + 1 == _date_to_jd(1582, 10, 15, 0.0)
    assert _jd_to_date(_date_to_jd(2000, 1, 1, 0.0, 'julian')) == (2000, 1, 14)
    # The Sothic anchor: 19 July (Julian) of 2781 BC is 26 June (Gregorian).
    assert _date_to_jd(-2780, 7, 19, 0.0, 'julian') == _date_to_jd(-2780, 6, 26, 0.0)


def test_julian_leap_years_have_no_century_rule():
    assert days_in_month(1900, 2, 'julian') == 29
    assert days_in_month(1900, 2, 'gregorian') == 28
    assert days_in_month(-1, 2, 'julian') == 28
    assert days_in_month(0, 2, 'julian') == 29


@pytest.mark.parametrize('year', [-12000, -2500, -400, -100, -1, 0, 1, 1900, 2000, 2023, 2024])
def test_days_in_month_agrees_with_julian_dates(year):
    for month in range(1, 13):
        next_year, next_month = (year + 1, 1) if month == 12 else (year, month + 1)
        length = _date_to_jd(next_year, next_month, 1, 0.0) - _date_to_jd(year, month, 1, 0.0)
        assert days_in_month(year, month) == round(length)
