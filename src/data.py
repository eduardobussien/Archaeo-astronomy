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

    # --- Aquila ---
    'Altair': {                                          # 12th brightest star; in Chinese myth the Cowherd star, separated from Vega by the Milky Way
        'ra': 297.696,  'dec':   8.868, 'dist':   5.13,
        'pm_ra':  536.82,  'pm_dec':  385.54,
        'mag': 0.77,  'constellation': 'Aquila',
    },

    # --- Gemini (addition) ---
    'Alhena': {                                          # Marks Gemini's foot; "the brand on a camel's neck" in Arabic
        'ra':  99.428,  'dec':  16.399, 'dist':  34.0,
        'pm_ra':    3.06,  'pm_dec':  -67.26,
        'mag': 1.93,  'constellation': 'Gemini',
    },

    # --- Canis Major (additions) ---
    'Wezen': {                                           # "The Weight" — a yellow-white supergiant 10,000× more luminous than the Sun
        'ra': 107.098,  'dec': -26.393, 'dist': 542.0,
        'pm_ra':   -2.73,  'pm_dec':    3.02,
        'mag': 1.84,  'constellation': 'Canis Major',
    },
    'Mirzam': {                                          # "The Announcer" — rises just before Sirius, heralding the brightest star
        'ra':  95.675,  'dec': -17.956, 'dist': 153.0,
        'pm_ra':   -3.55,  'pm_dec':   -0.61,
        'mag': 1.98,  'constellation': 'Canis Major',
    },

    # --- Carina (additions) ---
    'Avior': {                                           # Used in celestial navigation; marks the keel of the mythical ship Argo
        'ra': 125.628,  'dec': -59.510, 'dist': 200.0,
        'pm_ra':  -25.63,  'pm_dec':   22.72,
        'mag': 1.86,  'constellation': 'Carina',
    },
    'Miaplacidus': {                                     # Used by ancient Polynesian navigators; second brightest in Carina
        'ra': 138.300,  'dec': -69.717, 'dist':  34.1,
        'pm_ra': -157.66,  'pm_dec':  108.91,
        'mag': 1.67,  'constellation': 'Carina',
    },

    # --- Vela ---
    'Regor': {                                           # Brightest Wolf-Rayet star visible to the naked eye; once part of the ship Argo
        'ra': 122.383,  'dec': -47.337, 'dist': 260.0,
        'pm_ra':   -6.16,  'pm_dec':   10.14,
        'mag': 1.74,  'constellation': 'Vela',
    },
    'Suhail': {                                          # "The Bright One of the Ship" — used by Arab navigators of the Indian Ocean
        'ra': 136.999,  'dec': -43.433, 'dist': 168.0,
        'pm_ra':  -23.21,  'pm_dec':   14.29,
        'mag': 2.21,  'constellation': 'Vela',
    },

    # --- Puppis ---
    'Naos': {                                            # One of the hottest, most luminous naked-eye stars — a rare O-type supergiant
        'ra': 120.896,  'dec': -40.003, 'dist': 335.0,
        'pm_ra':  -30.82,  'pm_dec':   16.77,
        'mag': 2.25,  'constellation': 'Puppis',
    },

    # --- Hydra ---
    'Alphard': {                                         # "The Solitary One" — no bright neighbors; the Heart of the Sea Serpent
        'ra': 141.897,  'dec':  -8.658, 'dist':  54.0,
        'pm_ra':  -14.49,  'pm_dec':   33.25,
        'mag': 1.98,  'constellation': 'Hydra',
    },

    # --- Aries ---
    'Hamal': {                                           # The vernal equinox pointed here 2000 years ago, giving Aries its place as first sign
        'ra':  31.793,  'dec':  23.462, 'dist':  20.2,
        'pm_ra':  188.55,  'pm_dec': -148.08,
        'mag': 2.00,  'constellation': 'Aries',
    },

    # --- Perseus ---
    'Mirfak': {                                          # Heart of Perseus — the hero who slew Medusa and rescued Andromeda
        'ra':  51.081,  'dec':  49.861, 'dist': 155.0,
        'pm_ra':   24.11,  'pm_dec':  -26.01,
        'mag': 1.79,  'constellation': 'Perseus',
    },
    'Algol': {                                           # "The Demon Star" — its eerie dimming every 2.87 days was noticed and feared across cultures
        'ra':  47.042,  'dec':  40.956, 'dist':  28.2,
        'pm_ra':    2.39,  'pm_dec':   -1.44,
        'mag': 2.09,  'constellation': 'Perseus',
    },

    # --- Sagittarius ---
    'Kaus Australis': {                                  # Brightest in Sagittarius; the Milky Way's galactic core rises directly behind this star
        'ra': 276.043,  'dec': -34.385, 'dist':  44.5,
        'pm_ra':  121.17,  'pm_dec': -204.12,
        'mag': 1.79,  'constellation': 'Sagittarius',
    },
    'Nunki': {                                           # Babylonian "Yoke of the Sea" — they used it to predict flooding of the Euphrates
        'ra': 283.816,  'dec': -26.297, 'dist':  55.7,
        'pm_ra':   13.87,  'pm_dec':  -52.65,
        'mag': 2.05,  'constellation': 'Sagittarius',
    },

    # --- Pavo ---
    'Peacock': {                                         # One of the few stars with an official English proper name, assigned by British navigators in 1937
        'ra': 306.412,  'dec': -56.735, 'dist':  56.2,
        'pm_ra':    7.71,  'pm_dec':  -86.15,
        'mag': 1.94,  'constellation': 'Pavo',
    },

    # --- Triangulum Australe ---
    'Atria': {                                           # Brightest in the Southern Triangle; guides navigators toward the south celestial pole
        'ra': 252.166,  'dec': -69.028, 'dist': 127.0,
        'pm_ra':   17.99,  'pm_dec':  -32.92,
        'mag': 1.91,  'constellation': 'Triangulum Australe',
    },

    # --- Scorpius (additions) ---
    'Shaula': {                                          # "The Raised Tail" of the Scorpion; Polynesian navigators steered by this star across the Pacific
        'ra': 263.402,  'dec': -37.104, 'dist': 112.0,
        'pm_ra':   -3.09,  'pm_dec':  -29.95,
        'mag': 1.62,  'constellation': 'Scorpius',
    },
    'Sargas': {                                          # Ancient Sumerian name; part of the Scorpion's tail tracked by Mesopotamian astronomers
        'ra': 264.330,  'dec': -42.998, 'dist':  68.7,
        'pm_ra':   -5.58,  'pm_dec':  -25.71,
        'mag': 1.87,  'constellation': 'Scorpius',
    },
    'Dschubba': {                                        # "The Forehead of the Scorpion"; dramatically brightened in 2000 AD and remains variable
        'ra': 240.083,  'dec': -22.622, 'dist': 133.0,
        'pm_ra':   -8.36,  'pm_dec':  -26.58,
        'mag': 2.29,  'constellation': 'Scorpius',
    },

    # --- Ophiuchus ---
    'Rasalhague': {                                      # "Head of the Serpent Bearer" — Ophiuchus is the 13th zodiac constellation classical astrology omitted
        'ra': 263.734,  'dec':  12.560, 'dist':  14.6,
        'pm_ra':  108.07,  'pm_dec': -221.57,
        'mag': 2.08,  'constellation': 'Ophiuchus',
    },

    # --- Corona Borealis ---
    'Alphecca': {                                        # "The Jewel of the Northern Crown" — in Norse myth, the crown of Ariadne
        'ra': 233.672,  'dec':  26.715, 'dist':  22.9,
        'pm_ra':  120.35,  'pm_dec':  -89.17,
        'mag': 2.22,  'constellation': 'Corona Borealis',
    },

    # --- Andromeda ---
    'Alpheratz': {                                       # Corner star shared between Andromeda and Pegasus; anchors the Great Square navigation asterism
        'ra':   2.097,  'dec':  29.091, 'dist':  29.8,
        'pm_ra':  135.68,  'pm_dec': -162.95,
        'mag': 2.06,  'constellation': 'Andromeda',
    },

    # --- Pegasus ---
    'Enif': {                                            # "The Horse's Nose" — a bright orange supergiant marking the head of the winged Pegasus
        'ra': 326.046,  'dec':   9.875, 'dist': 211.0,
        'pm_ra':   26.92,  'pm_dec':    0.40,
        'mag': 2.38,  'constellation': 'Pegasus',
    },
    'Markab': {                                          # "The Saddle" of Pegasus — one corner of the Great Square used for navigation across cultures
        'ra': 346.190,  'dec':  15.205, 'dist':  42.5,
        'pm_ra':   61.10,  'pm_dec':  -42.56,
        'mag': 2.49,  'constellation': 'Pegasus',
    },

    # --- Corvus ---
    'Gienah': {                                          # Brightest in Corvus the Crow; Hipparchus used nearby stars to first measure precession of the equinoxes
        'ra': 183.786,  'dec': -17.542, 'dist':  50.0,
        'pm_ra': -158.61,  'pm_dec':   21.06,
        'mag': 2.59,  'constellation': 'Corvus',
    },

    # --- Cetus ---
    'Menkar': {                                          # "The Nostril of the Sea Monster"; in Persian star lore one of the four Royal Stars of Heaven
        'ra':  45.570,  'dec':   4.090, 'dist':  67.0,
        'pm_ra':  -11.81,  'pm_dec':  -78.32,
        'mag': 2.54,  'constellation': 'Cetus',
    },

    # --- Phoenix ---
    'Ankaa': {                                           # Brightest in Phoenix the Firebird; constellation introduced by Dutch navigators in the 1590s
        'ra':   6.571,  'dec': -42.306, 'dist':  23.8,
        'pm_ra':  232.76,  'pm_dec': -353.64,
        'mag': 2.40,  'constellation': 'Phoenix',
    },

    # --- Libra ---
    'Zubenelgenubi': {                                   # "The Southern Claw" — once part of Scorpius; detached into Libra by Roman astronomers
        'ra': 222.719,  'dec': -16.042, 'dist':  22.9,
        'pm_ra':   72.21,  'pm_dec': -103.82,
        'mag': 2.75,  'constellation': 'Libra',
    },

    # --- Crux (addition) ---
    'Gacrux': {                                          # Red giant at the top of the Southern Cross; key navigation star for indigenous Australians
        'ra': 187.792,  'dec': -57.113, 'dist':  27.2,
        'pm_ra':   27.94,  'pm_dec': -264.04,
        'mag': 1.59,  'constellation': 'Crux',
    },

    # --- Grus ---
    'Alnair': {                                          # "The Bright One" of Grus the Crane; used by Polynesian navigators in the south Pacific
        'ra': 332.058,  'dec': -46.961, 'dist':  31.0,
        'pm_ra':  127.60,  'pm_dec': -147.91,
        'mag': 1.73,  'constellation': 'Grus',
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
