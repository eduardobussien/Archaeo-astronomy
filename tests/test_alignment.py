import pytest

from alignment import format_year


@pytest.mark.parametrize('year, label', [
    (1, '1 AD'),
    (0, '1 BC'),
    (-1, '2 BC'),
    (-2499, '2500 BC'),
    (-10499, '10500 BC'),
])
def test_format_year_uses_astronomical_numbering(year, label):
    assert format_year(year) == label
