import astropy.units as u
import erfa
import numpy as np
import pytest
from astropy.coordinates import AltAz, EarthLocation, GeocentricTrueEcliptic, SkyCoord
from astropy.time import Time

from alignment import _date_to_jd, find_heliacal_rising
from data import STARS
from precession import J2000_JD, cio_locator, ecliptic_altaz, space_motion, star_altaz

GIZA_LAT, GIZA_LON = 29.9792, 31.1342
ARCMIN = 1.0 / 60.0
MAS = np.pi / 648_000_000.0
KM_S_TO_PC_PER_YEAR = 1.0227121650537077e-6


def _separation_deg(alt1, az1, alt2, az2):
    alt1, az1, alt2, az2 = map(np.radians, (alt1, az1, alt2, az2))
    cos_sep = (np.sin(alt1) * np.sin(alt2)
               + np.cos(alt1) * np.cos(alt2) * np.cos(az1 - az2))
    return np.degrees(np.arccos(np.clip(cos_sep, -1.0, 1.0)))


def _star(name, jd):
    s = STARS[name]
    return star_altaz(s['ra'], s['dec'], s['pm_ra'], s['pm_dec'], s['dist'],
                      jd, GIZA_LAT, GIZA_LON, s['rv'])


def _moved(name, jd, rv):
    s = STARS[name]
    ra, dec = space_motion(s['ra'], s['dec'], s['pm_ra'], s['pm_dec'], s['dist'], rv, jd)
    return np.degrees(erfa.seps(np.radians(s['ra']), np.radians(s['dec']), ra, dec))


def _pole_separation_deg(name, year):
    alt, az = _star(name, _date_to_jd(year, 1, 1, 0.0))
    return _separation_deg(alt, az, GIZA_LAT, 0.0)


def _astropy_frame(jd):
    t = Time(jd, format='jd', scale='ut1')
    t.delta_ut1_utc = 0.0
    loc = EarthLocation(lat=GIZA_LAT * u.deg, lon=GIZA_LON * u.deg)
    return t, AltAz(obstime=t, location=loc)


# astropy warns that ancient dates have no leap seconds, no measured polar motion and fall
# outside its aberration formula's tuned range; each matters only at the arcsecond level.
@pytest.mark.filterwarnings('ignore::erfa.ErfaWarning')
@pytest.mark.filterwarnings('ignore:Tried to get polar motions:astropy.utils.exceptions.AstropyWarning')
@pytest.mark.parametrize('year', [2000, 1000, 1, -999])
def test_all_stars_match_astropy_back_to_1000_bc(year):
    # astropy is an independent implementation (IAU 2006 precession, its own frames and space
    # motion). Its CIO locator comes from the IAU 2006 series, which drifts from its own
    # definition before about 1000 BC (2.4' at 2000 BC, 38' at 4000 BC; see the next test),
    # so the comparison stops here. Nutation and aberration, which astropy applies and this
    # engine does not, account for most of the remaining 0.2' to 0.7'.
    names = list(STARS)
    jd = _date_to_jd(year, 3, 20, 20.0)
    t, frame = _astropy_frame(jd)
    cat = {k: np.array([STARS[n][k] for n in names])
           for k in ('ra', 'dec', 'pm_ra', 'pm_dec', 'dist', 'rv')}
    ref = SkyCoord(ra=cat['ra'] * u.deg, dec=cat['dec'] * u.deg, distance=cat['dist'] * u.pc,
                   pm_ra_cosdec=cat['pm_ra'] * u.mas / u.yr, pm_dec=cat['pm_dec'] * u.mas / u.yr,
                   radial_velocity=cat['rv'] * u.km / u.s,
                   obstime=Time('J2000')).apply_space_motion(new_obstime=t).transform_to(frame)

    alt, az = star_altaz(cat['ra'], cat['dec'], cat['pm_ra'], cat['pm_dec'], cat['dist'],
                         jd, GIZA_LAT, GIZA_LON, cat['rv'])

    assert _separation_deg(alt, az, ref.alt.deg, ref.az.deg).max() < ARCMIN


def test_cio_locator_matches_its_definition_integrated_along_the_iau2006_pole():
    # s is defined by ds/dt = -(X dY/dt - Y dX/dt) / (1 + Z) along the path of the pole.
    # Integrating that with IAU 2006's own mean pole, back to 4000 BC, lands on this
    # engine's value; the IAU 2006 series for s (erfa.s06, used by astropy) is a
    # polynomial fitted near J2000 and is 38' off its own definition by then.
    jd = J2000_JD - np.arange(0.0, 6000 * 365.25 + 1.0, 365.25 / 4)
    x, y = erfa.bpn2xy(erfa.pmat06(jd, 0.0))
    rate = -(x * np.gradient(y, jd) - y * np.gradient(x, jd)) / (1.0 + np.sqrt(1.0 - x * x - y * y))
    s = erfa.s06(jd[0], 0.0, x[0], y[0]) + np.concatenate(
        [[0.0], np.cumsum((rate[1:] + rate[:-1]) / 2.0 * np.diff(jd))])
    for years_back in (1000, 2000, 3000, 4000, 5000, 6000):
        i = years_back * 4
        assert abs(cio_locator(erfa.epj(jd[i], 0.0)) - s[i]) < 0.1 * ARCMIN * np.pi / 180.0


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


@pytest.mark.parametrize('name', sorted(STARS))
def test_space_motion_is_a_straight_line_through_space(name):
    s = STARS[name]
    jd = _date_to_jd(-10499, 3, 20, 0.0)
    ra, dec = np.radians(s['ra']), np.radians(s['dec'])
    toward = erfa.s2c(ra, dec)
    east = np.array([-np.sin(ra), np.cos(ra), 0.0])
    north = np.array([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)])
    velocity = ((s['pm_ra'] * east + s['pm_dec'] * north) * MAS * s['dist']
                + s['rv'] * KM_S_TO_PC_PER_YEAR * toward)
    expected = erfa.c2s(s['dist'] * toward + velocity * (jd - J2000_JD) / 365.25)

    ours = space_motion(s['ra'], s['dec'], s['pm_ra'], s['pm_dec'], s['dist'], s['rv'], jd)
    assert np.degrees(erfa.seps(*ours, *expected)) * 3600 < 1.0


def test_radial_velocity_changes_apparent_speed_across_the_sky():
    # Approaching stars were farther away in the past, so they crossed the sky more slowly.
    jd = _date_to_jd(-10499, 3, 20, 0.0)
    assert STARS['Alpha Centauri']['rv'] < 0 < STARS['Aldebaran']['rv']
    assert _moved('Alpha Centauri', jd, STARS['Alpha Centauri']['rv']) < _moved('Alpha Centauri', jd, 0.0) - 2.0
    assert _moved('Aldebaran', jd, STARS['Aldebaran']['rv']) > _moved('Aldebaran', jd, 0.0)


def test_rejects_epochs_outside_model_range():
    with pytest.raises(ValueError):
        _star('Sirius', _date_to_jd(-30000, 1, 1, 0.0))
