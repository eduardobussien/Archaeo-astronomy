# Sources

Every number this app depends on, where it comes from, and the test that holds the code to it.

Sources are marked by strength:

- **Primary:** a peer-reviewed study, a survey, or an official ephemeris or catalogue.
- **Secondary:** a summary of someone else's work (an encyclopedia, a website, a popular article). It is used only where no primary value was available, and it is named so it can be replaced.

All tests are in `tests/` and run on every push (see the README). The accuracy each check found is collected in [ACCURACY.md](ACCURACY.md).

---

## 1. Calculation engine

| What | Source | Strength | Code | Checked by |
|---|---|---|---|---|
| Long-term precession of the equator and ecliptic | Vondrák, J., Capitaine, N. & Wallace, P. (2011), *New precession expressions, valid for long time intervals*, A&A 534, A22, [doi:10.1051/0004-6361/201117274](https://doi.org/10.1051/0004-6361/201117274) | Primary | `src/precession.py` (via ERFA `ltpb`, `ltpequ`, `ltpecl`, `ltecm`) | `tests/test_precession.py`: all 60 stars against astropy back to 1000 BC, the CIO locator against its definition back to 4000 BC, pole stars (Thuban, Vega) |
| Fundamental astronomy routines | IAU SOFA, through [ERFA / pyerfa](https://github.com/liberfa/pyerfa) | Primary | throughout | |
| Obliquity at J2000 (84381.406 arcsec) | IAU 2006 precession (Capitaine et al. 2003, A&A 412, 567) | Primary | `mean_obliquity` in `src/precession.py` | `tests/test_monuments.py`, `tests/test_alignment_finder.py` |
| Star positions and proper motions | Hipparcos new reduction: van Leeuwen, F. (2007), A&A 474, 653 | Primary | `STARS` in `src/data.py` | |
| Radial velocities of all 60 stars | [SIMBAD](https://simbad.cds.unistra.fr/) (Wenger et al. 2000, A&AS 143, 9) | Primary | `rv` fields in `src/data.py` | |
| Alpha Centauri barycentric motion | Kervella, P. et al. (2016), *Proper motions and orbits of Alpha Centauri*, A&A 594, A107 | Primary | `src/data.py` | |
| Sirius: Hipparcos values are barycentric | Bond, H. E. et al. (2017), *The Sirius system and its astrophysical puzzles*, ApJ 840, 70 | Primary | `src/data.py` | |
| Delta T (clock time vs uniform time) | Espenak, F. & Meeus, J. (2006), *Five Millennium Canon of Solar Eclipses*, NASA/TP-2006-214141, based on Morrison, L. V. & Stephenson, F. R. (2004), JHA 35, 327 | Primary | `src/timescales.py` | `tests/test_solar_system.py`: the NASA Delta T table |
| Delta T measured from 720 BC (earliest dated eclipse records) | Stephenson, F. R., Morrison, L. V. & Hohenkerk, C. Y. (2016), *Measurement of the Earth's rotation: 720 BC to AD 2015*, Proc. R. Soc. A 472, 20160404, [doi:10.1098/rspa.2016.0404](https://doi.org/10.1098/rspa.2016.0404) | Primary | `FIRST_MEASURED_YEAR` in `src/timescales.py` | |
| Uncertainty of Delta T (optimistic) | Morrison & Stephenson (2004), above: standard error 0.8 t^2 s, t in centuries from 1820 | Primary | `delta_t_uncertainty` in `src/timescales.py` | `tests/test_solar_system.py`: NASA's table |
| Uncertainty of Delta T (cautious) | Huber, P. J. (2000), *Modeling the length of day and extrapolating the rotation of the Earth*, in *Astronomical Amusements*, eds. Bonoli, De Meis & Panaino, Rome; both models tabulated by [NASA](https://eclipse.gsfc.nasa.gov/SEhelp/uncertainty2004.html) to 4000 BC | Primary (model), secondary (table) | `delta_t_uncertainty` in `src/timescales.py` | `tests/test_solar_system.py`: NASA's table |
| Sun, Moon and planets | ERFA `epv00`, `plan94`, `moon98` | Primary | `src/solar_system.py` | `tests/test_solar_system.py`: JPL DE441 |
| Reference ephemeris for the tests | Park, R. S. et al. (2021), *The JPL Planetary and Lunar Ephemerides DE440 and DE441*, AJ 161, 105; values from [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/) | Primary | test fixtures only (`tests/data/de441_moon_planets.json`) | `tests/test_solar_system.py` |
| Sunrise and sunset times (validation) | US Naval Observatory, [Astronomical Applications API](https://aa.usno.navy.mil/data/api), sunrise and sunset for 2025 at ten sites | Primary | test fixtures only (`tests/data/usno_sun_2025.json`) | `tests/test_alignment_finder.py` |
| Atmospheric refraction | Saemundsson, T. (1986), Sky & Telescope 72, 70; pressure from the standard atmosphere | Primary | `src/refraction.py` | `tests/test_refraction.py`: ERFA's refraction model |
| Heliacal visibility thresholds | Schaefer, B. E. (2000), *New methods and techniques for historical astronomy and archaeoastronomy*, Archaeoastronomy 15; fitted by V. Reijs, [archaeocosmology.org](http://www.archaeocosmology.org/eng/extinction.htm) | Primary (model), secondary (fit) | `src/visibility.py` | `tests/test_visibility.py`: Schaefer and [IMCCE](https://promenade.imcce.fr/en/pages6/724.html) for Sirius |
| Solstice dating from the obliquity | Lockyer, N. (1906), *Stonehenge and Other British Stone Monuments Astronomically Considered* | Primary (method) | `src/alignment_finder.py` | `tests/test_alignment_finder.py` |

---

## 2. Monuments

Each monument in `src/data.py` has an axis azimuth (`orientation_az`), the height of the visible horizon along it (`horizon_deg`, or none if unknown), and a `source`. An axis with no survey behind it is marked approximate, in the data and in the app.

### Great Pyramid of Giza

- **Axis:** 0 degrees (true north), within about 4 arcminutes.
- **Sources (primary):** Petrie, W. M. F. (1883), *The Pyramids and Temples of Gizeh*. Spence, K. (2000), *Ancient Egyptian chronology and the astronomical orientation of pyramids*, Nature 408, 320-324.

### Stonehenge

- **Axis:** 49.95 degrees true. **Horizon:** about 0.6 degrees.
- **Sources (primary):** Thom, A., Thom, A. S. & Thom, A. S. (1974), *Stonehenge*, JHA 5, 71-90. [Ruggles, C. (1997), *Astronomy and Stonehenge*, Proceedings of the British Academy 92, 203-229](https://www.thebritishacademy.ac.uk/documents/3920/92p203.pdf): best estimate 49.8 +/- 0.2 degrees, horizon between 0 and 1 degree.
- **Also consulted:** [MacKie, E., *A New Look at the Astronomy and Geometry of Stonehenge*, Culture and Cosmos 16](https://www.cultureandcosmos.org/pdfs/16/MacKie_INSAPVII_Stonehenge.pdf).
- **Cross-check (secondary):** [sarsen.org (2026), *81 degrees at Stonehenge*](https://www.sarsen.org/2026/10/81-degrees-at-stonehenge.html) gives 49.6 degrees for the midsummer sunrise in 2500 BC (disc centre, level horizon, no refraction). The app gives 49.59 (`tests/test_monuments.py`).

### Newgrange

- **Axis:** 136.05 degrees, the centre of the window through which sunlight reaches the chamber. **Horizon:** 1.3 degrees.
- **Sources (primary):** Patrick, J. (1974), *Midwinter sunrise at Newgrange*, Nature 249, 517-519: the window runs from 133 deg 42 min to 138 deg 24 min. O'Kelly, M. J. (1982), *Newgrange: Archaeology, Art and Legend*: on 21 December 1967 the beam entered at 08:58 GMT and was cut off at 09:15.
- **Also consulted (secondary):** [newgrange.com, winter solstice](https://www.newgrange.com/winter-solstice-newgrange.htm).
- **Check:** the Sun's azimuth during the timed 1967 beam (134.2 to 137.6 degrees) falls inside Patrick's window (`tests/test_monuments.py`). The horizon height is derived from the same timing.

### Pyramid of the Sun (Teotihuacan)

- **Axis:** 285 deg 28 min. **Horizon:** 1 deg 05 min, for an observer at the centre of the base.
- **Source (primary):** [Šprajc, I. (2000), *Astronomical alignments at Teotihuacan, Mexico*, Latin American Antiquity 11(4), 403-415](https://iaps.zrc-sazu.si/sites/default/files/Sprajc_AstronomicalAlignmentsTeotihuacan.pdf).
- **Site information (secondary):** [INAH, Visita Teotihuacán](https://visitateotihuacan.inah.gob.mx/en/explore/pyramid-of-the-sun/).
- **Check:** sunsets on 30 April and 13 August, exactly as Šprajc gives them (`tests/test_monuments.py`).

### El Castillo (Chichén Itzá)

- **Axis:** 291.72 degrees.
- **Source:** Šprajc, I. & Sánchez Nava, P. F. (2015), *Orientaciones astronómicas en la arquitectura maya de las tierras bajas*, INAH. The value here was taken from a secondary summary ([Satyori](https://satyori.com/ancient-mysteries/articles/chichen-itza-astronomical-alignments/)) and should be checked against the book.
- **Check:** sunsets on 20 May and 24 July, as published (`tests/test_monuments.py`).

### Angkor Wat

- **Axis:** 90.5 degrees, half a degree from true east-west.
- **Sources (primary):** [Mollerup, A. (2012), *Orientations of Khmer Temples*](https://cdn.angkordatabase.asia/libs/docs/publications/khmer-temples-and-the-astronomic-true-east/Mollerup2012-Orientations-Khmer-Temples.pdf). Stencel, R., Gifford, F. & Morón, E. (1976), *Astronomy and cosmology at Angkor Wat*, Science 193, 281-287: 0.75 degrees. [*Archaeoastronomy in the Khmer heartland* (arXiv:1604.05674)](https://arxiv.org/pdf/1604.05674).

### Göbekli Tepe

- **Axis:** 165 degrees (Enclosure C; 172 for D, 159 for B).
- **Sources (primary):** [Magli, G. (2016), *Sirius and the project of the megalithic enclosures at Göbekli Tepe*, Nexus Network Journal (arXiv:1307.8397)](https://arxiv.org/abs/1307.8397). [De Lorenzis, A. & Orofino, V. (2015), *New possible astronomic alignments at the megalithic site of Göbekli Tepe*](https://file.scirp.org/Html/5-1140033_53506.htm).
- **Caveat:** the pillars have not been surveyed with a theodolite, and the Sirius hypothesis is disputed.

### Mnajdra (Malta), lower temple

- **Axis:** 92 deg 30 min 11 s. **Horizon:** 3.5 degrees.
- **Source (primary):** [Micallef, P. I. (2000), *Alignments along the main axes at Mnajdra temples*, Xjenza 5(1)](https://www.xjenza.org/JOURNAL-XJENZA/2000-5-1/01-Alignments-Along-the-Main-Axes-at-Mnajdra-Temples.pdf): declination 0 on the axis. The horizon height is derived from the azimuth and that declination; another study gives about 4 degrees.
- **Check:** sunrise on the axis at the equinoxes (`tests/test_monuments.py`).

### Baalbek (Temple of Jupiter)

- **Axis:** 74.5 degrees. **Horizon:** 4 deg 16 min.
- **Sources (primary):** [*Archaeoastronomy of the temples of the Bekaa Valley*, Heritage 4 (2021)](https://doi.org/10.3390/heritage4030084). [*Archaeoastronomy and the chronology of the Temple of Jupiter at Baalbek* (arXiv:1606.05888)](https://arxiv.org/pdf/1606.05888): 75.5 degrees.
- **Coordinates (secondary):** [Wikipedia, Temple of Jupiter (Baalbek)](https://en.wikipedia.org/wiki/Temple_of_Jupiter_(Baalbek)).
- **Check:** sunrises within two days of 30 April and 12 August (`tests/test_monuments.py`).

### Machu Picchu

- **Axis:** 64.5 degrees, the Torreón window. **Horizon:** not known; set it before searching.
- **Source (primary):** Dearborn, D. S. P. & White, R. E. (1983), *The Torreón at Machu Picchu as an observatory*, Archaeoastronomy 5 (JHA 14).
- **Also consulted:** [*Astronomical observations at Machu Picchu: facts, hypothesis and wishful thinking* (ResearchGate)](https://www.researchgate.net/publication/365346202_Astronomical_Observations_at_Machu_Picchu_Facts_Hypothesis_and_Wishful_Thinking).

### Approximate axes (no survey found)

These azimuths are marked approximate in the data and in the app.

- **Avebury (50 degrees):** no surveyed axis; the Cove's midsummer sunrise claim is disputed (Ruggles 1999, *Astronomy in Prehistoric Britain and Ireland*). Consulted (secondary): [Stone Temple Gardening](https://www.stonetemplegardening.com/post/the-neolithic-world-and-the-astronomy-of-avebury-stone-circle-silbury-hill-west-kennet-long-barrow).
- **Carnac (83 degrees):** solstice and lunar claims (de Cleuziou 1887; Thom in the 1970s) are disputed (Ruggles 1999). Coordinates (secondary): [Wikipedia, Carnac stones](https://en.wikipedia.org/wiki/Carnac_stones).
- **Tiwanaku, Kalasasaya (89 degrees):** the equinox claim is not an established survey result. Consulted (secondary): [Benitez (2004), Archaeology field notes](https://interactive.archaeology.org/tiwanaku/fieldnotes/astronomy1.html).

---

## 3. What still needs a primary source

- El Castillo: confirm 291.72 degrees in Šprajc and Sánchez Nava (2015) itself.
- Avebury, Carnac and Tiwanaku: a surveyed axis.
- Machu Picchu: the horizon height along the Torreón window.
- Göbekli Tepe: a surveyed azimuth for each enclosure.
