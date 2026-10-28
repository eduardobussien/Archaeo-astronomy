# Investigating your own site

The app works for any place on Earth, not only the monuments it lists. An answer is only as good as what goes in, though, and the astronomy is now the most accurate part of the chain (see [ACCURACY.md](ACCURACY.md)). This guide covers the inputs that matter most, the traps in measuring them, and how to report a result so others can check it.

In the web interface, choose **Custom coordinates** under Ancient Site, then enter the axis and horizon in the **Alignment Finder**. The axis is drawn on the sky, and **Copy Link** saves the whole investigation.

---

## 1. Position and elevation

Latitude and longitude to four decimal places (about 10 metres) are more than enough. Elevation only sets the air pressure for refraction: at 3,850 m refraction is 38% weaker than at sea level, so a rough value is fine.

## 2. The axis azimuth: true north, not magnetic

The app needs the bearing of the axis from **true** (geographic) north, measured clockwise: 90 is east, 180 south, 270 west. This is the input most often wrong, and by more than every error in the astronomy put together.

**A compass gives magnetic north.** The difference, the magnetic declination, depends on the place and changes from year to year; in parts of North America it is more than 10 degrees. Correct the reading with the declination for the place and the day you measured, from [NOAA's calculator](https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml):

```
true bearing = magnetic bearing + declination      (declination east positive, west negative)
```

A magnetic bearing of 60 degrees where the declination is 10 degrees east is a true bearing of 70. Iron in nearby rock (basalt, some granites), in fences or in a vehicle can deflect a compass by degrees as well, so a compass is a check, not a survey.

**Map grids are not true north either.** National grids and UTM differ from true north by the grid convergence, up to a few degrees away from the grid's central meridian. Imagery in Google Earth is aligned to true north, so a bearing measured there between two points far apart along the axis is a true bearing.

**The Sun gives true north directly.** This is how surveyors do it, and the app can do the calculation:

1. On site, measure the horizontal angle from the axis to the Sun's centre with a theodolite or total station, and note the time to within a few seconds in UTC (a phone's GPS clock is good enough). The Sun's azimuth changes by about a quarter of a degree a minute in mid-latitudes, faster when it is high.
2. Convert the time to local mean solar time, which the app uses: `hour = UTC + longitude / 15`, with longitude east positive. If that passes midnight, add or subtract 24 hours and move the date by a day.
3. Ask the app where the Sun was: `/api/stars?lat=...&lon=...&year=...&month=...&day=...&hour=...`, and read `planets.Sun.azimuth`. Today the app's Sun is accurate to a few thousandths of a degree.
4. The axis azimuth is the Sun's azimuth plus or minus the angle you measured.

**How precise does it need to be?** For a star, the Alignment Finder reports how fast its rising point drifts. The date uncertainty is the azimuth uncertainty divided by that drift: Sirius drifting 0.46 degrees per century on an axis known to plus or minus half a degree gives a date only to within about 110 years either way. For the Sun, compare the azimuths the finder lists for neighbouring days.

## 3. The horizon height

The horizon height is the altitude of the visible skyline along the axis, as you would see it, in degrees. A distant level sea horizon is about 0; hills and mountains raise it; a viewpoint high above a flat plain can make it slightly negative (the app accepts -2 to 10).

It matters a great deal. A star rises later, and further round, over a raised horizon: each degree of horizon moves the Göbekli Tepe answer by about 210 years and the Teotihuacan one by about 75. A site among mountains, such as Machu Picchu, gives the wrong dates entirely with a level horizon.

- **On site:** measure it with a clinometer or theodolite along the axis. The angle you see already includes the bending of light by the air, which is what the app expects.
- **From elevation data:** [HeyWhatsThat](https://www.heywhatsthat.com) computes horizon profiles from SRTM data, including the Earth's curvature and the bending of light over the terrain. Its grid is about 90 m outside the United States, so it suits distant mountain horizons, not a skyline a few hundred metres away, trees or buildings.

## 4. Which part of the Sun

The Sun is half a degree across, so "sunrise on the axis" needs a definition. The finder offers the first gleam of the upper edge (the default), the centre of the disc, or the whole disc clear of the horizon. Use the one that fits the claim you are testing, and match the convention of any study you compare with: Šprajc's Mesoamerican dates, for example, are for the centre of the disc.

## 5. Stars are not seen at the horizon

The finder reports when a star's centre crossed the horizon you give it. But near the horizon the air dims starlight strongly, so a star is first seen some way above the skyline: higher for fainter stars and for hazier air (Schaefer 1986, *Atmospheric extinction effects on stellar alignments*, JHA Archaeoastronomy Supplement 17, S32). For anything fainter than the brightest stars, consider raising the horizon height by the altitude at which the star would first have been seen, and say that you did.

## 6. Dates, calendars and times

- Historians give dates before 1582 in the Julian calendar; the app shows either calendar, and switching relabels the same moment. The Sothic rising of Sirius was on 19 July (Julian) in 2781 BC, 26 June in the Gregorian calendar.
- The web interface shows historical years (2500 BC). The API uses astronomical years, in which 0 is 1 BC and -2499 is 2500 BC.
- Hours are local mean solar time: noon is when the Sun crosses the meridian, in any era. There were no time zones in antiquity, so this is the meaningful clock.

## 7. How far to trust the answer

[ACCURACY.md](ACCURACY.md) gives the measured accuracy for every kind of answer, era by era. In short: the stars and the Sun are far more accurate than any survey of a monument; the Moon's position at a given time is uncertain by about a degree in 2500 BC and by several degrees before 5000 BC, because the Earth's rotation is not measured before 720 BC; and local weather changes refraction near the horizon by a few arcminutes.

## 8. A match is not evidence

Over 14,000 years, with 60 bright stars, the Sun and the Moon, almost any axis lines up with something. Searching one invented axis (120 degrees, horizon 1.5, at 42.5 N) for 3000 to 2000 BC alone finds five stars that rose on it: Kaus Australis, Menkar, Rigel, Shaula and Sirius. An alignment becomes interesting when:

- the hypothesis comes first: decide which object and why before searching, from the culture's own records or iconography where they exist;
- the search is limited to the years the site was built and used (`from` and `to` in the API), from its archaeological dating;
- the same orientation recurs at several sites of the same culture, which is far harder to produce by chance (the method of Šprajc's Mesoamerican surveys and Ruggles's work in Britain and Ireland);
- the uncertainties in the axis and horizon are carried through to the date, and near misses are reported too.

Ruggles (1999), *Astronomy in Prehistoric Britain and Ireland*, is the standard discussion of these pitfalls.

## 9. Using the API from Python

Run the app locally (see the README), then call it from Python with only the standard library. This finds every star that rose on an axis in a period, and the days the Sun did:

```python
import json
import urllib.error
import urllib.parse
import urllib.request

BASE = 'http://127.0.0.1:5000'


def api(path, **params):
    url = f'{BASE}{path}?{urllib.parse.urlencode(params)}'
    try:
        with urllib.request.urlopen(url) as response:
            return json.load(response)
    except urllib.error.HTTPError as err:          # 400 or 404: the reason is in the body
        raise ValueError(json.load(err)['error']) from None


site = {'lat': 42.5, 'lon': 12.3, 'elevation': 300}

# Every star that rose on a 120 degree axis over a 1.5 degree horizon, 3000 to 2000 BC
stars = api('/api/stars', **site)['stars']
for name in sorted(stars):
    result = api('/api/alignment', star=name, event='rise', azimuth=120, horizon=1.5,
                 **{'from': -2999, 'to': -1999}, **site)
    for c in result.get('crossings', []):
        print(f"{name}: {1 - c['year']} BC, drifting {c['drift_deg_per_century']} deg per century")

# The days the Sun rose on the axis in 2500 BC
sun = api('/api/alignment', star='Sun', event='rise', azimuth=120, horizon=1.5, year=-2499, **site)
print([(d['month'], d['day']) for d in sun['dates']])
```

This prints:

```
Kaus Australis: 2769 BC, drifting 0.776 deg per century
Menkar: 2821 BC, drifting -0.799 deg per century
Rigel: 2014 BC, drifting -0.703 deg per century
Shaula: 2542 BC, drifting 0.803 deg per century
Sirius: 2484 BC, drifting -0.455 deg per century
[(1, 17), (11, 22)]
```

Any request can return a CSV file instead with `format=csv`. With pandas, `pd.read_csv(url, comment='#')` reads it, skipping the lines at the top that record how it was made. The full list of endpoints and parameters is in [API.md](API.md).

## 10. Reporting a result

So that someone else can check it, report:

- the latitude, longitude and elevation;
- the axis azimuth, how it was measured (and the magnetic declination used, if any), and its uncertainty;
- the horizon height, how it was measured, and its uncertainty;
- the part of the Sun, or the star and whether the horizon was raised for its visibility;
- refraction on or off, and the calendar of any dates;
- the app's version (`git rev-parse --short HEAD` in the repository) and the link from **Copy Link**, or the CSV files, whose first lines record all of the above.
