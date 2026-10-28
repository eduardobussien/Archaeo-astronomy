# How the calculation engine works

This is the detailed companion to the README's "How it works". How far to trust each answer, era by era and with the reference it was checked against, is collected in [ACCURACY.md](ACCURACY.md); every source is listed in [SOURCES.md](SOURCES.md).

## Stars

Star positions use a single model at every epoch (`src/precession.py`):

- **Precession:** Vondrák, Capitaine & Wallace (2011), valid to +/-200,000 years. The standard IAU 2006 polynomials are fitted to a few centuries of observations and drift by a third of a degree by 10,500 BC.
- **Earth rotation:** hour angles come from the Earth Rotation Angle, which is linear in UT1, measured from the Celestial Intermediate Origin. The CIO locator `s` is integrated numerically from the long-term pole, because the IAU 2006 series for `s` diverges beyond a few millennia.
- **Space motion:** each star moves in a straight line through space (`erfa.pmsafe`) from its J2000 position, using its proper motion and its radial velocity from [SIMBAD](https://simbad.cds.unistra.fr/). The radial velocity matters for nearby fast stars: an approaching star was farther away in the past and crossed the sky more slowly. Without it Alpha Centauri would be 2.2 degrees out of place at 10,500 BC, and Altair, Sirius and Procyon 0.07 to 0.14 degrees.
- **Binary stars:** over millennia a binary moves with its centre of mass, not with the orbital wobble of one component. Alpha Centauri uses the barycentric motion of Kervella et al. (2016). Sirius keeps its Hipparcos values, which the orbit of Bond et al. (2017) confirms are barycentric: they predict Sirius B's Gaia proper motion to within about 20 mas/yr.

The engine computes geometric mean places: nutation and aberration (each under 21 arcseconds) are not applied.

The star engine is checked in `tests/` against astropy, an independent implementation, for all 60 stars: they agree to within 0.7 arcminutes at every era from AD 2000 back to 1000 BC, most of it from nutation and aberration, which astropy applies and this engine leaves out. It is also checked against known pole stars (Thuban around 2830 BC, Vega around 12,000 BC), against the lowest culmination of Orion's Belt near 10,500 BC, and against the Sothic heliacal rising of Sirius in 2781 BC.

Astropy cannot serve as a reference before about 1000 BC, and the reason matters to anyone computing ancient skies with it. Turning the Earth's rotation into star positions needs the CIO locator `s`, which the IAU 2006 standard (and so astropy, and other software using the standard IAU 2006 routines) takes from a polynomial series fitted near the present. `s` is defined by an integral along the path of the celestial pole, and integrating that definition with IAU 2006's own pole shows how far the series drifts (arcminutes):

| Era | `s` from its definition (IAU 2006 pole) | This app | IAU 2006 series |
|---|---|---|---|
| 1000 BC | -16.16 | -16.16 | -16.49 |
| 2000 BC | -37.92 | -37.93 | -40.29 |
| 3000 BC | -73.06 | -73.08 | -84.12 |
| 4000 BC | -124.00 | -124.06 | -162.68 |

At 3000 BC the series puts every star about 10 arcminutes out of place, a third of the Moon's width, and at 4000 BC more than half a degree. This app integrates `s` from its definition along the long-term pole instead, and a test holds it to the integral above to within 0.1 arcminutes, back to 4000 BC.

## Atmospheric refraction

The air bends light, so everything near the horizon appears higher than it geometrically is: about 34 arcminutes at the horizon, falling to zero overhead. That shifts when and where objects appear to rise, which is exactly what alignment claims depend on. Displayed altitudes use the Saemundsson (1986) formula, which agrees with ERFA's refraction model to within 0.15 arcminutes above 10 degrees. Refraction scales with air pressure, which is taken from each site's elevation using the standard atmosphere: at Tiwanaku (3,850 m) it is 38% weaker than at sea level. Temperature and weather also change refraction near the horizon by a few arcminutes; standard conditions (10 C) are assumed. The API returns geometric altitudes with `refraction=0`.

## Sun, Moon and planets

The Earth's rotation has slowed over the millennia, so clock time (UT1, which follows the rotation) drifts away from the uniform time (TT) that orbital theories need. The difference, Delta T, is about 16.6 hours at 2500 BC. `src/timescales.py` uses the Espenak & Meeus (2006) expressions based on Morrison & Stephenson (2004); ignoring it would put the Moon about 8 degrees out of place at 2500 BC.

Positions come from ERFA's analytic theories (`epv00` for the Earth, `plan94` for the planets, `moon98` for the Moon), are rotated into the local sky by the same engine as the stars, and include the Moon's parallax. Measured against JPL's DE441 ephemeris on the first of every month (at the same instant of uniform time), the largest error in each era, in degrees, is:

| Era | Sun | Moon | Mercury | Venus | Mars | Jupiter | Saturn |
|---|---|---|---|---|---|---|---|
| 9998 BC | 0.30 | 1.75 | 0.47 | 0.51 | 0.18 | 1.23 | 1.47 |
| 7001 BC | 0.13 | 1.02 | 0.26 | 0.24 | 0.11 | 0.76 | 0.91 |
| 5001 BC | 0.06 | 0.44 | 0.11 | 0.10 | 0.04 | 0.10 | 0.70 |
| 3001 BC | 0.02 | 0.17 | 0.03 | 0.06 | 0.04 | 0.07 | 0.19 |
| 1001 BC | 0.005 | 0.04 | 0.01 | 0.02 | 0.01 | 0.03 | 0.11 |
| AD 1 to today | 0.002 | 0.02 | 0.003 | 0.003 | 0.009 | 0.015 | 0.055 |

The Moon's disc is about 0.5 degrees wide. From about 3000 BC on, every body is within 0.2 degrees. Earlier, the analytic theories (tuned to recent millennia) let the Moon, Jupiter and Saturn drift by 1 to 2 degrees. The Sun's declination, which sets where it rises, stays within 0.1 degrees in every era. Horizons provides the Sun and planets only from 9999 BC, so the app's last two millennia (to 12,000 BC) are not checked against it. The tests check every era in the table, and those between, with JPL's values in `tests/data/de441_moon_planets.json`.

Delta T itself is measured only from 720 BC, the earliest dated eclipse records; before that it is extrapolated, and uncertain by hours (see [What the Earth's rotation does and does not affect](#what-the-earths-rotation-does-and-does-not-affect)).

Far from the present, `epv00`'s orbital plane drifts off the true ecliptic, so in 9000 BC its Sun reached only +24.02 / -23.97 degrees of declination where DE441 has +-24.21. The real Sun stays within a few arcseconds of the ecliptic, so the app keeps `epv00`'s position along the ecliptic but places the Sun on the long-term (Vondrák) ecliptic of date. That brings the solstice declinations to +-24.20 and the Sun's total error at a 9000 BC test date from 0.21 to 0.18 degrees, and leaves recent dates unchanged.

## What the Earth's rotation does and does not affect

The app's clock is local solar time, which follows the Earth's rotation: in any era, local noon is when the Sun crosses the meridian. The stars follow the same rotation, so where a star stands at a local time does not depend on Delta T at all. The Sun, Moon and planets move along their orbits in uniform time, and to place them at a local time the app adds Delta T. Before the earliest dated eclipse records (Babylonian, from 720 BC) Delta T is not measured but extrapolated, and its uncertainty grows quickly:

| Era | Uncertainty in Delta T | The Sun moved by | The Moon moved by |
|---|---|---|---|
| 1000 BC | 10 minutes | 0.007 degrees | 0.1 degrees |
| 2500 BC | 25 minutes to 1.7 hours | 0.02 to 0.07 | 0.2 to 0.9 |
| 4000 BC | 45 minutes to 4.5 hours | 0.03 to 0.19 | 0.4 to 2.5 |
| 6000 BC | 1.4 to 10 hours | 0.06 to 0.4 | 0.8 to 5.7 |
| 10,000 BC | 3 to 29 hours | 0.13 to 1.2 | 1.7 to 16 |

The lower values follow Morrison & Stephenson (2004), whose standard error is 0.8 t^2 seconds, with t in centuries from 1820. The higher ones follow Huber's (2000) model, which treats changes in the Earth's rotation as a random walk; [NASA tabulates both](https://eclipse.gsfc.nasa.gov/SEhelp/uncertainty2004.html) back to 4000 BC. The rows for 6000 and 10,000 BC extend the same formulas beyond the range they were published for, so they are only indicative.

What that means for the answers the app gives:

- **Not affected:** star positions, and where and when stars rose and set, in every era; the years the alignment finder gives for stars; the solstice sunrise azimuths, and the era in which a solstice sunrise fell on an axis, which depend only on the tilt of the Earth's axis.
- **Little affected:** the days of the year on which the Sun rose on an axis, and heliacal rising dates. Until 4000 BC the Sun's place in its yearly path is uncertain by at most a fifth of a degree (about five hours of its motion). By 10,000 BC the cautious estimate reaches about a degree, so a date could be a day out.
- **Affected:** the Moon at a given local time. Under the cautious estimate it is uncertain by about a degree in 2500 BC and by several degrees before 5000 BC, and in deep prehistory the night of a full moon cannot be fixed. The limits of the Moon's monthly swing north and south, which lunar standstill studies use, come from the geometry of its orbit rather than from timing, and are not affected. The planets move at most a sixth as fast as the Moon across the sky and are affected correspondingly less.

`/api/stars` reports this for every request in `meta.delta_t`: Delta T in seconds, its uncertainty under both estimates, how far that moves the Moon, and whether the date falls in the measured period (from 720 BC). Before 720 BC the web interface adds the cautious figure to the Moon's tooltip and info card.

The uncertainty in Delta T is larger than the error of the Moon theory itself (1.75 degrees at 10,000 BC, from the table above), so a more precise ephemeris would not make the prehistoric Moon much more certain.

## Alignment finder

Where a star rises on the horizon depends only on the site's latitude and the star's declination: with the star's centre at geometric altitude h0 as it appears, cos A = (sin d - sin phi sin h0) / (cos phi cos h0) gives the rising azimuth A, and it sets at 360 - A. Precession changes every star's declination over the centuries, so its rising point sweeps along the horizon: seen from Göbekli Tepe around 8800 BC, Sirius rose 1.5 degrees further east each century. `src/alignment_finder.py` evaluates this for every year from 12,000 BC to today in one vectorized pass and reports the years in which the rising (or setting) point crossed the axis. h0 is the visible horizon lowered by refraction. For an axis that points at hills the horizon height can be raised, and it matters: each degree of horizon moves the Göbekli Tepe answer by about 210 years and the Teotihuacan one by about 75.

For each crossing it also finds a night that year when the event happens near local midnight, so the sky can be shown with the star on the axis. The tests check the search by asking it for the year whose azimuth it was given, and check each night against the engine directly: the star stands on the horizon there to within 0.01 degrees, at the axis azimuth to within 0.02 degrees. For example, Sirius rose on the 165 degree axis of Göbekli Tepe in 8813 BC, within the site's period of use.

The Sun is different: its rising point swings between the solstices every year rather than drifting over centuries. For the Sun the finder lists the days of the chosen year when sunrise (or sunset) fell on the axis, found from every sunrise of that year, and the range between the two solstices. It also asks which era's solstice sunrise fell exactly on the axis, the method Norman Lockyer used to date Stonehenge in the early 1900s: the solstice declination equals the obliquity of the ecliptic, which shrank from 24.2 degrees around 8000 BC to 23.44 today, so the solstice points move along the horizon. The obliquity matches IAU 2006 at J2000 (84381.406 arcseconds, falling 46.8 per century). These two routes, a closed formula with the obliquity and the year's sunrises from the solar ephemeris, agree on the solstice azimuths to 0.01 degrees from 9000 BC to today, and the tests hold them to it. Which part of the Sun counts as sunrise is a choice: the first gleam (default), the centre or the whole disc. Refraction is applied to the edge of the disc where it meets the horizon, and the centre lies one solar radius from it. With that, the app's sunrise and sunset times for 2025 agree with the US Naval Observatory's at ten sites from 59 degrees N to 46 degrees S, at every season, to within their rounding to the minute (160 times, none more than 33 seconds apart). Refracting the centre instead had put sunrise up to 74 seconds early.

## Horizon profiles

A horizon profile is the skyline's altitude in every direction, interpolated linearly between surveyed points (`src/horizon.py`). With a profile, objects count as visible only above the skyline in their direction, and sunrise and sunset are the first and last gleam over it: the Sun can slip behind a peak and show again in a notch (a double sunset), and the last gleam is the one that counts. Before the Sun reaches the skyline's lowest point it is hidden everywhere, and once past the highest it is clear everywhere, so the search samples the time between and bisects the first (or last) change. The tests check it against a second-by-second scan over a sawtooth skyline, to within 2 seconds. A heliacal rising also waits until the star clears the skyline, and a star's alignment with an axis depends only on the skyline's height at the axis.

## Monument data

Every monument in `src/data.py` carries its axis azimuth, the height of the visible horizon along it where one is known, and the source of both. Where no survey supports an axis (Avebury, Carnac, Tiwanaku) the azimuth is marked approximate, in the data and in the app. The alignment finder uses each site's measured horizon by default. Every source, for the monuments and for the engine, is listed with its strength in [SOURCES.md](SOURCES.md).

Checking the original values against the literature found and fixed several errors: Stonehenge's axis was a flat-horizon sunset calculation (51.2 degrees) rather than the surveyed 49.95 degrees (Thom et al. 1974); Baalbek's azimuth was 12 degrees off and its latitude 0.2 degrees off; Carnac's coordinates pointed about 8 km away; and several notes described the wrong building or an unsupported claim.

With the corrected data the engine reproduces published results, which the tests in `tests/test_monuments.py` hold it to:

| Site | Published | This app |
|---|---|---|
| Teotihuacan, Pyramid of the Sun (Šprajc 2000) | sunsets on 30 April and 13 August | 30 April, 13 August |
| Chichén Itzá, El Castillo (Šprajc and Sánchez Nava 2015) | sunsets on 20 May and 24 July | 20 May, 24 July |
| Baalbek, Temple of Jupiter (Heritage 2021) | sunrises on about 30 April and 12 August | 1 May, 13 August |
| Mnajdra, lower temple (Micallef 2000) | equinox sunrise | 22 March, 23 September |
| Newgrange (Patrick 1974; O'Kelly 1967 timing) | light in the chamber for azimuths 133.7 to 138.4 degrees | the timed 1967 beam ran from 134.2 to 137.6 degrees |
| Stonehenge (sarsen.org, Vondrák obliquity) | midsummer sunrise 49.6 degrees in 2500 BC, disc centre, level horizon | 49.59 degrees |

## Heliacal rising

A star's heliacal rising is the first morning it can be glimpsed in the dawn twilight after weeks hidden in the Sun's glare. How dark the sky must be depends on the star's brightness, and how high the star must climb depends on how much the air dims it near the horizon (extinction, in magnitudes per airmass). `src/visibility.py` uses Reijs's fits to Schaefer's visibility model ([archaeocosmology.org](http://www.archaeocosmology.org/eng/extinction.htm)): for Sirius in typical air (0.27, Schaefer's estimate for ancient Memphis) the Sun must be 6 degrees down and the star 4 degrees up, matching Schaefer (2000) and [IMCCE](https://promenade.imcce.fr/en/pages6/724.html); a magnitude 2 star needs the Sun 12 degrees down. The web interface offers clear, typical and hazy air.

Every dawn from the previous October to the end of the target year is checked at once: a vectorized bisection finds when the Sun reaches the star's darkness limit each morning, and the first dawn in the target year on which the star is high enough, after one on which it was not, is reported. Stars that are seen every dawn (like Thuban, the pole star around 2800 BC) or never (Canopus from Giza in hazy air) are reported as such. The model's quoted uncertainty (about 1.5 to 2 degrees in each threshold) corresponds to a few days in the date.
