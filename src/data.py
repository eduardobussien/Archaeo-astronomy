# All star and monument reference data.
#
# Star fields (J2000 epoch):
#   ra       - Right Ascension in degrees
#   dec      - Declination in degrees
#   dist     - Distance in parsecs
#   pm_ra    - Proper motion in RA (mas/yr, already includes cos(dec) factor)
#   pm_dec   - Proper motion in Dec (mas/yr)
#   mag      - Apparent visual magnitude (lower = brighter; Sirius is -1.46)
#   constellation - IAU constellation name

STARS = {
    # --- Orion ---
    'Betelgeuse': {
        'ra': 88.793,   'dec':  7.407,  'dist':  168.0,
        'pm_ra':   27.54,  'pm_dec':   11.30,
        'mag': 0.42,  'constellation': 'Orion',
    },
    'Rigel': {
        'ra': 78.634,   'dec': -8.202,  'dist':  264.0,
        'pm_ra':    1.87,  'pm_dec':   -0.56,
        'mag': 0.13,  'constellation': 'Orion',
    },
    'Alnitak': {
        'ra': 85.190,   'dec': -1.942,  'dist':  226.0,
        'pm_ra':    3.99,  'pm_dec':    2.54,
        'mag': 1.74,  'constellation': 'Orion',
    },
    'Alnilam': {
        'ra': 84.053,   'dec': -1.202,  'dist':  606.0,
        'pm_ra':    1.44,  'pm_dec':   -0.73,
        'mag': 1.69,  'constellation': 'Orion',
    },
    'Mintaka': {
        'ra': 83.002,   'dec': -0.299,  'dist':  280.0,
        'pm_ra':    0.18,  'pm_dec':   -0.58,
        'mag': 2.23,  'constellation': 'Orion',
    },

    # --- Canis Major ---
    'Sirius': {
        'ra': 101.287,  'dec': -16.716, 'dist':    2.64,
        'pm_ra': -546.01, 'pm_dec': -1223.07,
        'mag': -1.46, 'constellation': 'Canis Major',
    },

    # --- Carina ---
    'Canopus': {
        'ra': 95.988,   'dec': -52.696, 'dist':   95.9,
        'pm_ra':   19.93,  'pm_dec':   23.24,
        'mag': -0.72, 'constellation': 'Carina',
    },

    # --- Taurus ---
    'Aldebaran': {
        'ra': 68.980,   'dec': 16.509,  'dist':   20.0,
        'pm_ra':   62.78,  'pm_dec': -189.36,
        'mag': 0.86,  'constellation': 'Taurus',
    },
    'Alcyone': {                                     # Brightest of the Pleiades
        'ra': 56.871,   'dec': 24.105,  'dist':  116.0,
        'pm_ra':   19.35,  'pm_dec':  -43.11,
        'mag': 2.87,  'constellation': 'Taurus',
    },

    # --- Lyra ---
    'Vega': {                                        # Future north star (~13,700 AD)
        'ra': 279.235,  'dec': 38.784,  'dist':    7.68,
        'pm_ra':  200.94,  'pm_dec':  286.23,
        'mag': 0.03,  'constellation': 'Lyra',
    },

    # --- Draco ---
    'Thuban': {                                      # North star ~2700 BC (pyramid era)
        'ra': 211.097,  'dec': 64.376,  'dist':  120.0,
        'pm_ra':  -56.52,  'pm_dec':   17.19,
        'mag': 3.65,  'constellation': 'Draco',
    },

    # --- Ursa Minor ---
    'Polaris': {                                     # Current north star
        'ra': 37.955,   'dec': 89.264,  'dist':  133.0,
        'pm_ra':   44.22,  'pm_dec':  -11.74,
        'mag': 1.98,  'constellation': 'Ursa Minor',
    },

    # --- Cygnus ---
    'Deneb': {
        'ra': 310.358,  'dec': 45.280,  'dist':  802.0,
        'pm_ra':    1.56,  'pm_dec':    1.55,
        'mag': 1.25,  'constellation': 'Cygnus',
    },

    # --- Boötes ---
    'Arcturus': {
        'ra': 213.915,  'dec': 19.182,  'dist':   11.26,
        'pm_ra': -1093.45, 'pm_dec': -1999.40,    # One of the fastest-moving bright stars
        'mag': -0.05, 'constellation': 'Boötes',
    },

    # --- Leo ---
    'Regulus': {
        'ra': 152.093,  'dec': 11.967,  'dist':   24.31,
        'pm_ra': -248.73,  'pm_dec':    5.59,
        'mag': 1.35,  'constellation': 'Leo',
    },

    # --- Virgo ---
    'Spica': {                                       # Hipparchus used this to discover precession
        'ra': 201.298,  'dec': -11.161, 'dist':   77.0,
        'pm_ra':  -42.50,  'pm_dec':  -31.73,
        'mag': 1.04,  'constellation': 'Virgo',
    },

    # --- Scorpius ---
    'Antares': {
        'ra': 247.352,  'dec': -26.432, 'dist':  185.0,
        'pm_ra':  -12.11,  'pm_dec':  -23.30,
        'mag': 1.06,  'constellation': 'Scorpius',
    },

    # --- Auriga ---
    'Capella': {                                         # 6th brightest star; circumpolar from mid-latitudes
        'ra':  79.172,  'dec':  45.998, 'dist':   12.94,
        'pm_ra':   75.52,  'pm_dec': -427.13,
        'mag': 0.08,  'constellation': 'Auriga',
    },

    # --- Gemini ---
    'Castor': {
        'ra': 113.650,  'dec':  31.888, 'dist':   15.27,
        'pm_ra': -206.33,  'pm_dec': -148.18,
        'mag': 1.58,  'constellation': 'Gemini',
    },
    'Pollux': {                                          # Brightest star in Gemini; has a planet
        'ra': 116.329,  'dec':  28.026, 'dist':   10.34,
        'pm_ra': -626.55,  'pm_dec':  -45.80,
        'mag': 1.14,  'constellation': 'Gemini',
    },

    # --- Canis Minor ---
    'Procyon': {                                         # 8th brightest; rises just before Sirius
        'ra': 114.826,  'dec':   5.225, 'dist':    3.51,
        'pm_ra': -716.57,  'pm_dec': -1034.58,
        'mag': 0.34,  'constellation': 'Canis Minor',
    },

    # --- Canis Major (second brightest after Sirius) ---
    'Adhara': {
        'ra': 104.656,  'dec': -28.972, 'dist':  132.0,
        'pm_ra':    2.63,  'pm_dec':    3.33,
        'mag': 1.50,  'constellation': 'Canis Major',
    },

    # --- Eridanus ---
    'Achernar': {                                        # Marks the southern end of the river; key southern nav star
        'ra':  24.429,  'dec': -57.237, 'dist':   44.1,
        'pm_ra':   88.02,  'pm_dec':  -40.08,
        'mag': 0.46,  'constellation': 'Eridanus',
    },

    # --- Piscis Austrinus ---
    'Fomalhaut': {                                       # "Autumn Star"; one of the four Royal Stars of Persia
        'ra': 344.413,  'dec': -29.622, 'dist':    7.69,
        'pm_ra':  328.95,  'pm_dec': -164.67,
        'mag': 1.16,  'constellation': 'Piscis Austrinus',
    },

    # --- Centaurus ---
    'Alpha Centauri': {                                  # Closest star system to the Sun (4.37 ly)
        'ra': 219.902,  'dec': -60.835, 'dist':    1.34,
        'pm_ra': -3678.19,  'pm_dec':  481.84,
        'mag': -0.27, 'constellation': 'Centaurus',
    },
    'Hadar': {                                           # Beta Centauri; forms the Southern Pointer with Alpha Cen
        'ra': 210.956,  'dec': -60.373, 'dist':  120.0,
        'pm_ra':  -33.27,  'pm_dec':  -22.78,
        'mag': 0.61,  'constellation': 'Centaurus',
    },

    # --- Crux (Southern Cross) ---
    'Acrux': {                                           # Brightest star in the Southern Cross
        'ra': 186.650,  'dec': -63.099, 'dist':   98.3,
        'pm_ra':  -35.83,  'pm_dec':  -14.73,
        'mag': 0.77,  'constellation': 'Crux',
    },
    'Mimosa': {                                          # Second brightest in the Southern Cross
        'ra': 191.930,  'dec': -59.689, 'dist':  108.0,
        'pm_ra':  -48.22,  'pm_dec':  -12.14,
        'mag': 1.25,  'constellation': 'Crux',
    },

    # --- Ursa Major ---
    'Dubhe': {                                           # The outer pointer star to Polaris; part of the Big Dipper
        'ra': 165.932,  'dec':  61.751, 'dist':   37.9,
        'pm_ra': -134.11,  'pm_dec':  -34.72,
        'mag': 1.79,  'constellation': 'Ursa Major',
    },
}


# Monument orientation azimuths are the compass direction (degrees, 0=N clockwise)
# of the monument's primary astronomical sightline.  These are the directions you
# would look from the monument to observe the aligned celestial event.
MONUMENTS = {
    'Great Pyramid of Giza': {
        'lat': 29.9792, 'lon': 31.1342,
        'orientation_az': 0.0,
        'note': (
            'Sides aligned to true north within ~0.05°. '
            "King's Chamber south shaft targets Orion's Belt at transit (~2450 BC)."
        ),
    },
    'Stonehenge': {
        'lat': 51.1789, 'lon': -1.8262,
        'orientation_az': 51.2,
        'note': (
            'Heel stone marks midsummer sunrise (~51.2° in present era). '
            'The opposing axis (~231°) aligns with midwinter sunset.'
        ),
    },
    'Angkor Wat': {
        'lat': 13.4125, 'lon': 103.8670,
        'orientation_az': 90.0,
        'note': (
            'Main east-west axis unique for a Hindu temple (west-facing). '
            'Equinox sunrise aligns over the central tower as seen from the western entrance.'
        ),
    },
    'Pyramid of the Sun (Teotihuacan)': {
        'lat': 19.6925, 'lon': -98.8438,
        'orientation_az': 285.5,
        'note': (
            'West face oriented toward the Pleiades setting point (~285.5°). '
            'The Avenue of the Dead runs 15.5° east of north.'
        ),
    },
    'El Castillo (Chichen Itza)': {
        'lat': 20.6829, 'lon': -88.5686,
        'orientation_az': 25.7,
        'note': (
            'NNE staircase axis aligns with Venus at maximum northern elongation. '
            'Equinox serpent-shadow effect on the north staircase.'
        ),
    },
    'Göbekli Tepe': {
        'lat': 37.2232, 'lon': 38.9225,
        'orientation_az': 165.0,
        'note': (
            'Oldest known monumental complex (~9600 BC), predating agriculture. '
            'Enclosure D pillars may have targeted Sirius rising (~165°) at foundation date.'
        ),
    },
    'Newgrange': {
        'lat': 53.6947, 'lon': -6.4753,
        'orientation_az': 136.4,
        'note': (
            'Neolithic passage tomb (~3200 BC). The roof-box admits a narrow beam of '
            'winter solstice sunrise light that illuminates the inner chamber for 17 minutes.'
        ),
    },
    'Avebury': {
        'lat': 51.4292, 'lon': -1.8536,
        'orientation_az': 50.0,
        'note': (
            'Largest megalithic stone circle in the world (~2600 BC). '
            'The Kennet Avenue leads toward midsummer sunrise (~50°) from the southern entrance.'
        ),
    },
    'Carnac Stones': {
        'lat': 47.5989, 'lon': -2.9539,
        'orientation_az': 83.0,
        'note': (
            'Over 3,000 menhirs arranged in parallel rows (~4500-3300 BC). '
            'The Le Menec alignment runs roughly east (~83°), broadly toward equinox sunrise.'
        ),
    },
    'Mnajdra Temple (Malta)': {
        'lat': 35.8269, 'lon': 14.4367,
        'orientation_az': 94.0,
        'note': (
            'Neolithic temple (~3600 BC), among the oldest free-standing structures on Earth. '
            'Lower temple doorway precisely frames equinox sunrise; solstice sunlight hits side walls.'
        ),
    },
    'Machu Picchu': {
        'lat': -13.1631, 'lon': -72.5449,
        'orientation_az': 65.0,
        'note': (
            'Inca citadel (~1450 AD). The Torreon (Temple of the Sun) aligns with the '
            'June solstice sunrise (~65°). The Intihuatana stone served as a solar calendar.'
        ),
    },
    'Tiwanaku (Kalasasaya)': {
        'lat': -16.5544, 'lon': -68.6742,
        'orientation_az': 89.0,
        'note': (
            'Pre-Inca ceremonial center (~500-900 AD, possibly older). '
            'The eastern gate of the Kalasasaya temple is precisely aligned with the equinox sunrise.'
        ),
    },
    'Baalbek (Temple of Jupiter)': {
        'lat': 34.2036, 'lon': 36.2100,
        'orientation_az': 62.0,
        'note': (
            'Massive Roman temple complex built on a far older Phoenician foundation. '
            'The temple axis (~62°) aligns northeast; the monumental podium stones remain unexplained.'
        ),
    },
}
