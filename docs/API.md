# API reference

The Flask server (`python app.py`) serves the web interface and a JSON API on `http://127.0.0.1:5000`. [RESEARCH_GUIDE.md](RESEARCH_GUIDE.md) shows how to call it from Python with only the standard library.

| Endpoint | Parameters | Description |
|---|---|---|
| `GET /api/stars` | `lat, lon, elevation, year, month, day, hour, calendar, refraction, horizon, site, format` | Star, Sun, Moon and planet alt-az positions plus the ecliptic, everything one sky view needs |
| `GET /api/ecliptic` | same as `/api/stars` | Only the 73 ecliptic great-circle points in alt-az |
| `GET /api/heliacal` | `lat, lon, year, star, extinction, horizon, site` | First heliacal rising date for a star |
| `POST /api/horizon` | a horizon profile file as the request body | The profile in the compact form `horizon=` takes, with its point count and range |
| `GET /api/alignment` | `site` or `lat, lon, elevation`; `star, event, azimuth, horizon, from, to, calendar, refraction, format`; for `star=Sun` also `year, limb` | Years in which a star rose or set along an axis, a night to see each, and the drift curve; for the Sun, the days of the year on the axis, the solstice range and the solstice eras on the axis |
| `GET /api/sites` | none | List of all monument sites with coordinates and orientation notes |

## Dates and years

Dates are in the proleptic Gregorian calendar unless `calendar=julian` is given; `/api/stars` returns the same local date in both calendars in `meta.dates`, and `/api/heliacal` answers in the calendar requested. Historians give ancient dates in the Julian calendar (the Sothic rising of Sirius on 19 July), which by 2781 BC runs 23 days ahead of the Gregorian (26 June); the web interface can show either and switching relabels the same moment.

`year` uses astronomical convention: -2500 = 2501 BC, 0 = 1 BC, 1 = 1 AD. The web interface shows and accepts historical years instead (-2500 = 2500 BC, no year 0) and converts them, so permalinks carry the astronomical value. `hour` is local mean solar time.

## Sites, validation and refraction

`site` must be an exact monument name from `/api/sites` (case does not matter). Invalid input returns HTTP 400 (404 for an unknown site) with a JSON body such as `{"error": "'lat' must be between -90.0 and 90.0"}`. Dates are checked against the real month length, so 29 February is accepted only in leap years. With custom coordinates, `elevation` (metres, default 0) sets the air pressure for refraction; `refraction=0` returns geometric altitudes on `/api/stars`, `/api/ecliptic` and `/api/heliacal`.

## Horizon

`horizon` is the skyline: a number for a level horizon (default 0; on `/api/alignment` the monument's measured horizon along its axis), or a horizon profile as `azimuth:altitude` pairs, `horizon=0:1.5,90:3,180:0.5,270:2`, interpolated linearly between points. With a profile, objects count as visible only above the skyline in their direction, sunrise and sunset are the first and last gleam over it (so a double sunset, the Sun reappearing in a notch, ends at its last gleam), and a heliacal rising waits until the star clears it. A star's alignment depends only on the skyline at the axis. `POST /api/horizon` turns a profile file (one `azimuth altitude` point per line, as surveys and Stellarium horizon files hold them) into that form: `curl --data-binary @skyline.txt http://localhost:5000/api/horizon`.

## CSV

`format=csv` on `/api/stars` and `/api/alignment` returns the results as a CSV file instead of JSON: one row per star and body, or for an alignment one table whose `section` column says what each row is (a crossing, a date, a solstice, a point on the drift curve). Each file begins with `#` lines recording the site, date, settings and the exact request, so it can be regenerated; spreadsheets show them as text rows, and pandas skips them with `pd.read_csv(path, comment='#')`. The web interface offers the same files through its CSV buttons.
