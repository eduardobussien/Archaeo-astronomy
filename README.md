# Archaeo-Astronomy Sky Explorer

An interactive web app for reconstructing the ancient night sky at historical sites, accounting for stellar proper motion and axial precession over thousands of years.

![Sky chart at Giza, 2500 BC](docs/screenshot_hero.png)

---

## Features

- **Polar sky chart:** North-up, horizon-to-zenith radial view powered by Plotly.js
- **60 historically significant stars** with proper motion applied back to any epoch
- **13 ancient monument sites** (Giza, Stonehenge, Angkor Wat, Göbekli Tepe, Machu Picchu, and more) with known astronomical orientations
- **Custom coordinates:** enter any latitude/longitude for off-catalog sites
- **Magnitude limit slider:** filter stars by brightness (adjustable from 0 to 6.5)
- **Planet positions:** Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn in alt-az
- **Ecliptic overlay:** dotted great circle showing the path of the Sun/Moon/planets
- **Constellation lines:** stick figures connecting stars within each constellation
- **Star/planet labels:** toggle name labels on all objects
- **Time animation:** step or play forward/backward in 1h / 1d / 1mo / 1yr increments
- **Precession sweep:** watch the North Celestial Pole trace its 26,000-year arc; animate in 10yr / 100yr / 1kyr steps
- **Heliacal rising finder:** find the exact date a star first appears at dawn after a period of solar conjunction
- **URL permalink:** every state (site, date, toggles) is encoded in the URL; share or bookmark any view

![Heliacal rising result for Sirius at Giza](docs/screenshot_heliacal.png)

---

## How to run locally

**1. Clone the repo and create a virtual environment**

```bash
git clone https://github.com/your-username/Archaeo-astronomy.git
cd Archaeo-astronomy
python -m venv venv
```

**2. Activate the virtual environment**

On Windows:
```bash
venv\Scripts\activate
```

On macOS / Linux:
```bash
source venv/bin/activate
```

**3. Install dependencies**

```bash
pip install -r requirements.txt
```

**4. Run the Flask server**

```bash
python app.py
```

**5. Open in your browser**

Navigate to `http://127.0.0.1:5000`

---

## Project structure

```
app.py              Flask server and API endpoints
src/
  alignment.py      Site and date calculations (star table, ecliptic, heliacal rising)
  precession.py     Long-term precession, Earth rotation and proper motion engine
  solar_system.py   Sun, Moon and planets from ERFA's analytic theories
  timescales.py     Delta T (difference between uniform time and Earth rotation time)
  data.py           Star catalog (60 stars, J2000 positions + proper motion) and monument list
  visualize.py      Static chart export helpers
templates/
  index.html        Single-page app, all UI and Plotly.js rendering
tests/
  test_precession.py    Stars: astropy, pole stars, Orion's Belt, Sirius
  test_solar_system.py  Sun, Moon, planets: JPL DE441 and astropy; Delta T: NASA table
  test_alignment.py     Year labels and calendar conversion
requirements.txt
```

---

## Calculation engine

Star positions use a single model at every epoch (`src/precession.py`):

- **Precession:** Vondrák, Capitaine & Wallace (2011), valid to +/-200,000 years. The standard IAU 2006 polynomials are fitted to a few centuries of observations and drift by a third of a degree by 10,500 BC.
- **Earth rotation:** hour angles come from the Earth Rotation Angle, which is linear in UT1, measured from the Celestial Intermediate Origin. The CIO locator `s` is integrated numerically from the long-term pole, because the IAU 2006 series for `s` diverges beyond a few millennia.
- **Proper motion:** rigorous 3D space motion (`erfa.pmsafe`) from the J2000 catalog.

Positions are geometric mean places: nutation, aberration and refraction (all under half a degree) are not applied.

The model is checked in `tests/` against astropy near J2000 (under 1 arcmin), against known pole stars (Thuban around 2830 BC, Vega around 12,000 BC), against the lowest culmination of Orion's Belt near 10,500 BC, and against the Sothic heliacal rising of Sirius in 2781 BC.

### Sun, Moon and planets

The Earth's rotation has slowed over the millennia, so clock time (UT1, which follows the rotation) drifts away from the uniform time (TT) that orbital theories need. The difference, Delta T, is about 16.6 hours at 2500 BC. `src/timescales.py` uses the Espenak & Meeus (2006) expressions based on Morrison & Stephenson (2004); ignoring it would put the Moon about 8 degrees out of place at 2500 BC.

Positions come from ERFA's analytic theories (`epv00` for the Earth, `plan94` for the planets, `moon98` for the Moon), are rotated into the local sky by the same engine as the stars, and include the Moon's parallax. Measured against JPL's DE441 ephemeris at the same instant, the largest error across all seven bodies is:

| Epoch | Largest error |
|---|---|
| 1000 BC | 0.09 degrees |
| 3000 BC | 0.18 degrees |
| 4000 BC | 0.35 degrees |
| 9000 BC | 1.3 degrees (Moon) |

Delta T itself is uncertain for prehistoric dates: published long-term models differ by about two hours at 10,500 BC, enough to move the Moon by about a degree. The Sun and planets move slowly enough that this barely matters for them.

**Heliacal rising** is found by checking every dawn from the previous October to the end of the target year at once: a vectorized bisection finds when the Sun reaches the arc-vision threshold (default -10 degrees) each morning, and the first dawn in the target year on which the star is visible after a dawn on which it was not is reported.

To run the tests:

```bash
pip install pytest
python -m pytest
```

![Precession sweep, NCP traces its 26,000-year arc](docs/screenshot_precession.png)

---

## API endpoints

| Endpoint | Parameters | Description |
|---|---|---|
| `GET /api/stars` | `lat, lon, year, month, day, hour, site` | Star + planet alt-az positions |
| `GET /api/ecliptic` | same as `/api/stars` | 73 ecliptic great-circle points in alt-az |
| `GET /api/heliacal` | `lat, lon, year, star, arc_vision, site` | First heliacal rising date for a star |
| `GET /api/sites` | none | List of all monument sites with coordinates and orientation notes |

`year` uses astronomical convention: -2500 = 2501 BC, 0 = 1 BC, 1 = 1 AD. The web interface shows and accepts historical years instead (-2500 = 2500 BC, no year 0) and converts them, so permalinks carry the astronomical value.

---

## Dependencies

- **Flask:** web server
- **pyerfa:** IAU SOFA routines: the Vondrák 2011 long-term precession and the Sun, Moon and planet theories
- **astropy:** independent reference values in the tests
- **numpy:** vectorized star and time calculations
- **Plotly.js** (CDN): interactive polar chart in the browser
