"""
Draw docs/orion_belt_giza.png: how high Orion's Belt climbed over Giza each
night, from 12,000 BC to today.

A star is highest when it crosses the meridian (culminates); for a star south
of the zenith that altitude is 90 - (latitude - declination). Precession
changes the declination over the centuries, so the culmination altitude traces
a slow wave with a 26,000-year period.

    python scripts/orion_belt_figure.py
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'src'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from alignment import _date_to_jd, format_year
from data import MONUMENTS, STARS
from precession import star_declination
from refraction import apparent_altitude

BELT = ('Alnitak', 'Alnilam', 'Mintaka')
COLORS = {'Alnitak': '#4fc3f7', 'Alnilam': '#f0c040', 'Mintaka': '#ef9a9a'}
PYRAMID_BUILT = -2559          # 2560 BC, the reign of Khufu
OUTPUT = os.path.join(ROOT, 'docs', 'orion_belt_giza.png')


def culmination_altitudes(star, years, lat, elevation_m):
    """Apparent altitude (degrees) of the star at culmination, mid-year of each year."""
    s = STARS[star]
    jd = np.array([_date_to_jd(int(y), 7, 2, 0.0) for y in years])
    dec = star_declination(s['ra'], s['dec'], s['pm_ra'], s['pm_dec'], s['dist'], s['rv'], jd)
    return apparent_altitude(90.0 - (lat - dec), elevation_m)


def main():
    giza = MONUMENTS['Great Pyramid of Giza']
    years = np.arange(-11999, 2027, 10)
    altitudes = {star: culmination_altitudes(star, years, giza['lat'], giza['elevation_m'])
                 for star in BELT}
    lowest = years[np.argmin(altitudes['Alnilam'])]
    at_build = {star: float(np.interp(PYRAMID_BUILT, years, alt)) for star, alt in altitudes.items()}

    fig, ax = plt.subplots(figsize=(11, 5.6), facecolor='#0d1117')
    ax.set_facecolor('#0d1117')
    for star in BELT:
        ax.plot(years, altitudes[star], color=COLORS[star], linewidth=2, label=star)

    ax.axvline(lowest, color='#8b949e', linestyle='--', linewidth=1)
    ax.annotate(f'Lowest culmination: {1 - lowest:,} BC\n'
                f'Alnilam {altitudes["Alnilam"].min():.1f} deg above the horizon',
                xy=(lowest, altitudes['Alnilam'].min()), xytext=(lowest + 900, 4),
                color='#c9d1d9', fontsize=10,
                arrowprops={'arrowstyle': '->', 'color': '#8b949e'})
    ax.axvline(PYRAMID_BUILT, color='#a8834f', linestyle=':', linewidth=1.5)
    ax.annotate(f'Great Pyramid built, c. 2560 BC: Alnilam {at_build["Alnilam"]:.1f} deg,\n'
                "about the slope of the King's Chamber southern shaft\n(Badawy and Trimble 1964)",
                xy=(PYRAMID_BUILT, at_build['Alnilam']), xytext=(PYRAMID_BUILT - 6400, 50),
                color='#c9d1d9', fontsize=10,
                arrowprops={'arrowstyle': '->', 'color': '#a8834f'})

    ax.set_xlim(years[0], years[-1])
    ax.set_ylim(0, 62)
    ticks = [-11999, -9999, -7999, -5999, -3999, -1999, 1, years[-1]]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{1 - t:,} BC' for t in ticks[:-2]] + ['AD 1', 'today'], color='#8b949e')
    ax.tick_params(colors='#8b949e')
    ax.set_ylabel('Altitude at culmination (deg)', color='#c9d1d9')
    ax.grid(color='#21262d', linewidth=0.6)
    for spine in ax.spines.values():
        spine.set_color('#30363d')
    ax.legend(loc='upper left', facecolor='#161b22', edgecolor='#30363d', labelcolor='#c9d1d9')
    ax.set_title("How high Orion's Belt climbed over Giza each night, 12,000 BC to today",
                 color='#e6edf3', fontsize=13, pad=12)
    fig.text(0.5, 0.005,
             'Precession (Vondrak 2011), proper motion and refraction included. The belt reached the '
             'bottom of its 26,000-year cycle around the date the "Orion correlation" theory gives;\n'
             'that is a fact about precession, and on its own says nothing about when, or why, the '
             'pyramids were planned.',
             ha='center', va='bottom', color='#8b949e', fontsize=8.5)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(OUTPUT, dpi=150, facecolor=fig.get_facecolor())
    print(f'Saved {os.path.relpath(OUTPUT, ROOT)}; lowest culmination {format_year(int(lowest))}, '
          f'Alnilam at {at_build["Alnilam"]:.2f} deg in 2560 BC')


if __name__ == '__main__':
    main()
