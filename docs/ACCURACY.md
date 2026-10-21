# Accuracy

How far to trust each answer the app gives, measured against independent references. Every number on this page is held by a test in `tests/`, which runs on every push, so it cannot quietly get worse. Where each reference comes from is in [SOURCES.md](SOURCES.md).

The app covers 12,000 BC to today. Errors grow with distance from the present, so each table gives them by era.

---

## At a glance

| Question | How accurate | What limits it |
|---|---|---|
| Where was a star at a given local time? | Within 0.7 arcminutes back to 1000 BC (checked against astropy); the same method holds to 12,000 BC | The precession model (Vondrák 2011), designed for +-200,000 years |
| In which year did a star rise on an axis? | Each reported night puts the star on the axis to within 0.02 degrees | The survey: each degree of horizon height moves the answer by decades to centuries |
| On which days did the Sun rise on an axis? | The published dates for Teotihuacan and El Castillo exactly, Baalbek within 1 day | Before 4000 BC, the Earth's rotation: up to about a day by 10,000 BC |
| What time was sunrise or sunset? | Within 33 seconds of the US Naval Observatory today | Weather changes refraction near the horizon by a few arcminutes |
| In which era was a solstice sunrise on an axis? | Within 0.01 degrees of the year's actual sunrises | Only the tilt of the Earth's axis; no timing involved |
| Where was the Moon? | Within 0.2 degrees back to 3000 BC; 1.75 degrees at 10,000 BC | Before 720 BC, the Earth's rotation: about 1 degree at 2500 BC, several degrees before 5000 BC |
| Where were the planets? | Within 0.2 degrees back to 3000 BC; up to 1.5 degrees at 10,000 BC | The planetary theory, tuned to recent millennia |
| When was a star's heliacal rising? | Sirius matches Schaefer (2000) and IMCCE | The visibility model: its uncertainty is a few days |

---

## Stars

The app moves each star through space from its catalogue position, then applies long-term precession and the Earth's rotation.

**Against astropy**, an independent implementation, all 60 stars at Giza (`tests/test_precession.py`):

| Era | Largest difference |
|---|---|
| AD 2000 | 0.44 arcminutes |
| AD 1000 | 0.47 |
| AD 1 | 0.39 |
| 1000 BC | 0.67 |

Most of the difference is nutation and aberration, which astropy applies and this app leaves out (each under 21 arcseconds).

**Before 1000 BC** astropy stops being a reference. It takes the CIO locator `s`, which ties the Earth's rotation to the sky, from the IAU 2006 polynomial series, and that series drifts from its own definition: 2.4 arcminutes at 2000 BC, 11 at 3000 BC, 38 at 4000 BC. This app integrates `s` from the definition instead, and integrating the definition along IAU 2006's own pole reproduces the app's value to within 0.06 arcminutes at 4000 BC (`tests/test_precession.py`). The README has the table.

**Over the whole range**, the tests check known historical facts:

- Thuban within 0.5 degrees of the pole around 2830 BC.
- Vega near the pole around 12,000 BC.
- Orion's Belt culminating lowest near 10,500 BC.
- The Sothic heliacal rising of Sirius in 2781 BC.

**Star motions** use SIMBAD radial velocities; without them Alpha Centauri would be 2.2 degrees out of place at 10,500 BC.

---

## The Sun

**Against JPL's DE441 ephemeris** on the first of every month (`tests/test_solar_system.py`):

| Era | Largest error | Largest error in declination |
|---|---|---|
| 9998 BC | 0.30 degrees | 0.084 degrees |
| 9001 BC | 0.23 | 0.063 |
| 7001 BC | 0.13 | 0.036 |
| 5001 BC | 0.06 | 0.018 |
| 3001 BC | 0.02 | 0.007 |
| 1001 BC | 0.005 | 0.001 |
| AD 1 to today | 0.002 | 0.0002 |

Declination sets where the Sun rises, and it stays within 0.1 degrees in every era, below the precision of any monument survey. The Sun's disc is 0.53 degrees wide. The solstice declinations match DE441 to 0.02 degrees in 9000 BC and 2500 BC.

**Sunrise and sunset times against the US Naval Observatory** for 2025, at ten sites from 59 degrees N to 46 degrees S on eight dates (`tests/test_alignment_finder.py`): all 160 times agree to within 33 seconds. The Observatory rounds to the minute, so a perfect match can be 30 seconds off.

---

## The Moon and planets

**Against JPL's DE441** on the first of every month, largest error in degrees (`tests/test_solar_system.py`):

| Era | Moon | Mercury | Venus | Mars | Jupiter | Saturn |
|---|---|---|---|---|---|---|
| 9998 BC | 1.75 | 0.47 | 0.51 | 0.18 | 1.23 | 1.47 |
| 7001 BC | 1.02 | 0.26 | 0.24 | 0.11 | 0.76 | 0.91 |
| 5001 BC | 0.44 | 0.11 | 0.10 | 0.04 | 0.10 | 0.70 |
| 3001 BC | 0.17 | 0.03 | 0.06 | 0.04 | 0.07 | 0.19 |
| 1001 BC | 0.04 | 0.01 | 0.02 | 0.01 | 0.03 | 0.11 |
| AD 1 to today | 0.02 | 0.003 | 0.003 | 0.009 | 0.015 | 0.055 |

The Moon's disc is about 0.5 degrees wide. The analytic theories are tuned to recent millennia, so before about 5000 BC the Moon, Jupiter and Saturn drift by 1 to 2 degrees.

---

## The Earth's rotation (Delta T)

The app's clock is local solar time, which follows the Earth's rotation, as the stars do. Placing the Sun, Moon and planets at a local time needs Delta T, which is measured only from 720 BC. Before that it is extrapolated, with this published uncertainty (`tests/test_solar_system.py` checks both formulas against NASA's table):

| Era | Uncertainty in Delta T | The Sun moved by | The Moon moved by |
|---|---|---|---|
| 1000 BC | 10 minutes | 0.007 degrees | 0.1 degrees |
| 2500 BC | 25 minutes to 1.7 hours | 0.02 to 0.07 | 0.2 to 0.9 |
| 4000 BC | 45 minutes to 4.5 hours | 0.03 to 0.19 | 0.4 to 2.5 |
| 6000 BC | 1.4 to 10 hours | 0.06 to 0.4 | 0.8 to 5.7 |
| 10,000 BC | 3 to 29 hours | 0.13 to 1.2 | 1.7 to 16 |

The lower figures are Morrison & Stephenson (2004), the higher Huber (2000); rows before 4000 BC extend the formulas and are only indicative. Not affected: anything about the stars, solstice azimuths and solstice eras. The `/api/stars` response reports this for every date, and the app shows it beside the Moon before 720 BC.

---

## Atmospheric refraction

Saemundsson's formula agrees with ERFA's refraction model to within 0.15 arcminutes above 10 degrees altitude (`tests/test_refraction.py`), and is scaled for each site's elevation. Standard conditions (10 C) are assumed. Real weather changes refraction near the horizon by a few arcminutes, so no sunrise or star rising can be pinned down more finely than that.

---

## Monuments: published results the app reproduces

With each site's surveyed axis and horizon (`tests/test_monuments.py`):

| Site | Published | This app |
|---|---|---|
| Teotihuacan, Pyramid of the Sun (Šprajc 2000) | sunsets on 30 April and 13 August | 30 April, 13 August |
| Chichén Itzá, El Castillo (Šprajc and Sánchez Nava 2015) | sunsets on 20 May and 24 July | 20 May, 24 July |
| Baalbek, Temple of Jupiter (Heritage 2021) | sunrises on about 30 April and 12 August | 1 May, 13 August |
| Mnajdra, lower temple (Micallef 2000) | equinox sunrise | 22 March, 23 September |
| Newgrange (Patrick 1974; O'Kelly 1967 timing) | sunlight in the chamber for azimuths 133.7 to 138.4 degrees | the timed 1967 beam ran from 134.2 to 137.6 degrees |
| Stonehenge (sarsen.org) | midsummer sunrise at 49.6 degrees in 2500 BC | 49.59 degrees |

---

## Heliacal rising

How dark the sky must be, and how high the star, follows Schaefer's visibility model as fitted by Reijs. For Sirius in typical air the Sun must be 6 degrees down and the star 4 degrees up, matching Schaefer (2000) and IMCCE (`tests/test_visibility.py`). The model's own uncertainty, about 1.5 to 2 degrees in each threshold, corresponds to a few days in the date.

---

## What is not checked

- **10,000 to 12,000 BC for the Sun, Moon and planets.** JPL's online service starts at 9999 BC. The errors there are expected to continue the trend in the tables, but that has not been measured.
- **Star positions before 1000 BC against an independent program.** The method is checked against its own definition back to 4000 BC, and against the historical facts above across the range, but no independent code reaches further.
- **Nutation and aberration** are left out. Each is under 21 arcseconds, far below anything an ancient alignment can resolve.
- **Local horizons.** The app takes the horizon height you give it. Where a site's horizon has not been surveyed (marked in the app), every rising and setting answer is only as good as that guess.
