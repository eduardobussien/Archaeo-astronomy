import astropy.units as u
import erfa
import numpy as np
import pytest
from astropy.coordinates import AltAz, EarthLocation, GeocentricTrueEcliptic, SkyCoord
from astropy.time import Time

from alignment import _date_to_jd, find_heliacal_rising
from data import STARS
from precession import cio_locator, ecliptic_altaz, star_altaz

GIZA_LAT, GIZA_LON = 29.9792, 31.1342
ARCMIN = 1.0 / 60.0
MAS = np.pi / 648_000_000.0


def _separation_deg(alt1, az1, alt2, az2):
    alt1, az1, alt2, az2 = map(np.radians, (alt1, az1, alt2, az2))
    cos_sep = (np.sin(alt1) * np.sin(alt2)
               + np.cos(alt1) * np.cos(alt2) * np.cos(az1 - az2))
    return np.degrees(np.arccos(np.clip(cos_sep, -1.0, 1.0)))


def _star(name, jd):
    s = STARS[name]
    return star_altaz(s['ra'], s['dec'], s['pm_ra'], s['pm_dec'], s['dist'],
                      jd, GIZA_LAT, GIZA_LON)


def _pole_separation_deg(name, year):
    alt, az = _star(name, _date_to_jd(year, 1, 1, 0.0))
    return _separation_deg(alt, az, GIZA_LAT, 0.0)


def _astropy_frame(jd):
    t = Time(jd, format='jd', scale='ut1')
    t.delta_ut1_utc = 0.0
    loc = EarthLocation(lat=GIZA_LAT * u.deg, lon=GIZA_LON * u.deg)
    return t, AltAz(obstime=t, location=loc)


def test_stars_match_astropy_near_j2000():
    names = ['Sirius', 'Alnitak', 'Deneb', 'Canopus', 'Polaris', 'Arcturus']
    jd = _date_to_jd(2000, 3, 20, 20.0)
    t, frame = _astropy_frame(jd)
    cat = {k: np.array([STARS[n][k] for n in names]) for k in ('ra', 'dec', 'pm_ra', 'pm_dec', 'dist')}
    ref = SkyCoord(ra=cat['ra'] * u.deg, dec=cat['dec'] * u.deg, distance=cat['dist'] * u.pc,
                   pm_ra_cosdec=cat['pm_ra'] * u.mas / u.yr, pm_dec=cat['pm_dec'] * u.mas / u.yr,
                   obstime=Time('J2000')).apply_space_motion(new_obstime=t).transform_to(frame)

    alt, az = star_altaz(cat['ra'], cat['dec'], cat['pm_ra'], cat['pm_dec'], cat['dist'],
                         jd, GIZA_LAT, GIZA_LON)

    assert _separation_deg(alt, az, ref.alt.deg, ref.az.deg).max() < ARCMIN


def test_ecliptic_matches_astropy_near_j2000():
    jd = _date_to_jd(2000, 3, 20, 20.0)
    t, frame = _astropy_frame(jd)
    lam = np.arange(0.0, 360.0, 15.0)
    ref = SkyCoord(lon=lam * u.deg, lat=np.zeros_like(lam) * u.deg,
                   frame=GeocentricTrueEcliptic(equinox=t)).transform_to(frame)

    alt, az = ecliptic_altaz(jd, GIZA_LAT, GIZA_LON, lam)

    assert _separation_deg(alt, az, ref.alt.deg, ref.az.deg).max() < ARCMIN


@pytest.mark.parametrize('year', [1800, 1900, 1950, 2050, 2100, 2200])
def test_cio_locator_matches_iau2006_mean_pole(year):
    jd = 2451545.0 + (year - 2000) * 365.25
    x, y = erfa.bpn2xy(erfa.pmat06(jd, 0.0))
    assert abs(cio_locator(erfa.epj(jd, 0.0)) - erfa.s06(jd, 0.0, x, y)) < 20 * MAS


@pytest.mark.parametrize('star, year, max_sep_deg', [
    ('Polaris', 2000, 1.0),
    ('Thuban', -2830, 0.5),
    ('Vega', -12000, 5.0),
])
def test_historical_pole_stars(star, year, max_sep_deg):
    assert _pole_separation_deg(star, year) < max_sep_deg


def test_orion_belt_culminates_lowest_near_10500_bc():
    day = np.linspace(0.0, 1.0, 2881)[:, None]
    culminations = {}
    for year in (-11500, -10500, -9500):
        alt, _ = _star('Alnitak', _date_to_jd(year, 1, 1, 0.0) + day)
        culminations[year] = alt.max()
    assert min(culminations, key=culminations.get) == -10500


def test_sirius_heliacal_rising_matches_sothic_anchor():
    # Julian 19 July, 2780 BC is proleptic Gregorian 26 June.
    result = find_heliacal_rising(GIZA_LAT, GIZA_LON, -2780, 'Sirius')
    assert result['month'] == 6
    assert abs(result['day'] - 26) <= 2
    # Sirius at dec -21.5 deg rises near azimuth 115 deg at Giza's latitude.
    assert 110.0 < result['star_azimuth'] < 122.0


def test_rejects_epochs_outside_model_range():
    with pytest.raises(ValueError):
        _star('Sirius', _date_to_jd(-30000, 1, 1, 0.0))
