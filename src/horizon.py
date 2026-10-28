"""
A horizon profile: the altitude of the visible skyline in every direction.

A profile is a list of (azimuth, altitude) points in degrees: azimuth from true
north through east, altitude as seen from the site (refraction included, as a
clinometer or theodolite reads it). Between points the skyline is interpolated
linearly in azimuth, wrapping through north; a single point is a level horizon.

Profiles arrive as text, one point per line, as a survey spreadsheet or a
Stellarium horizon file holds them ('120.5, 2.3', '120.5<tab>2.3' or
'120.5 2.3'; '#' starts a comment, an optional first line names the columns
and further columns are ignored), or compactly in a URL as 'az:alt,az:alt,...'.

Everywhere the app takes a horizon it takes either a profile or a number, the
altitude of a level horizon (as_profile turns one into the other).
"""
import re

import numpy as np

MAX_POINTS = 720               # one every half degree
ALTITUDE_RANGE = (-2.0, 45.0)  # below the sea horizon from a cliff, up to a steep mountainside

_SEPARATOR = re.compile(r'[,;\s]+')


class HorizonProfile:
    def __init__(self, azimuths, altitudes):
        az = np.asarray(azimuths, dtype=float)
        alt = np.asarray(altitudes, dtype=float)
        if az.ndim != 1 or az.shape != alt.shape or not 1 <= len(az) <= MAX_POINTS:
            raise ValueError(f'A horizon profile needs 1 to {MAX_POINTS} points')
        if not (np.isfinite(az).all() and np.isfinite(alt).all()):
            raise ValueError('Horizon azimuths and altitudes must be numbers')
        low, high = ALTITUDE_RANGE
        if not ((low <= alt) & (alt <= high)).all():
            raise ValueError(f'Horizon altitudes must be between {low:g} and {high:g} degrees')
        if not ((0.0 <= az) & (az <= 360.0)).all():
            raise ValueError('Horizon azimuths must be between 0 and 360 degrees')
        az = az % 360.0
        order = np.argsort(az, kind='stable')
        self.azimuths, self.altitudes = az[order], alt[order]
        if (np.diff(self.azimuths) == 0.0).any():
            raise ValueError('Each azimuth may appear only once in a horizon profile '
                             '(0 and 360 count as the same direction)')

    def __call__(self, azimuth):
        """Altitude of the skyline (degrees) at the given azimuth(s)."""
        if len(self.azimuths) == 1:
            return np.full_like(np.asarray(azimuth, dtype=float), self.altitudes[0])
        return np.interp(np.asarray(azimuth, dtype=float) % 360.0,
                         self.azimuths, self.altitudes, period=360.0)

    @property
    def is_level(self):
        """True when the skyline is at the same altitude in every direction."""
        return bool(np.ptp(self.altitudes) == 0.0)

    def outline(self, step=1.0):
        """(azimuths, altitudes) tracing the skyline once round, through every
        surveyed point and every step degrees between, for drawing."""
        az = np.append(np.union1d(np.arange(0.0, 360.0, step), self.azimuths), 360.0)
        return az, self(az)

    def compact(self):
        """The profile as 'az:alt,az:alt,...', for a URL."""
        return ','.join(f'{a:g}:{h:g}' for a, h in zip(self.azimuths, self.altitudes))

    @classmethod
    def from_compact(cls, text):
        try:
            pairs = [p.split(':') for p in text.strip().split(',')]
            az, alt = zip(*((float(a), float(h)) for a, h in pairs))
        except ValueError:
            raise ValueError("A compact horizon profile is 'azimuth:altitude' pairs separated "
                             "by commas, such as '0:1.5,90:3,180:0.5,270:2'") from None
        return cls(az, alt)

    @classmethod
    def from_text(cls, text):
        """One 'azimuth altitude' point per line; see the module docstring."""
        az, alt = [], []
        first_line = True
        for number, line in enumerate(text.splitlines(), start=1):
            fields = [f for f in _SEPARATOR.split(line.split('#', 1)[0].strip()) if f]
            if not fields:
                continue
            try:
                values = [float(f) for f in fields]
            except ValueError:
                if first_line:
                    first_line = False
                    continue                      # a header such as 'azimuth,altitude'
                raise ValueError(f'Line {number} of the horizon profile is not numbers: '
                                 f'{line.strip()!r}') from None
            first_line = False
            if len(values) < 2:
                raise ValueError(f'Line {number} of the horizon profile needs an azimuth and '
                                 f'an altitude: {line.strip()!r}')
            az.append(values[0])                  # further columns, such as distance, are ignored
            alt.append(values[1])
        if not az:
            raise ValueError('The horizon profile has no points')
        return cls(az, alt)


def as_profile(horizon):
    """A HorizonProfile as given, or a level one at the altitude of a number."""
    return horizon if isinstance(horizon, HorizonProfile) else HorizonProfile([0.0], [horizon])
