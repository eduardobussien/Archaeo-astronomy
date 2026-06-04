import argparse
import warnings
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from astropy.utils.exceptions import AstropyWarning
import erfa

from data import STARS, MONUMENTS
from alignment import calculate_alignments, check_alignments, get_observation_time

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=erfa.ErfaWarning)

# One color per constellation group
CONSTELLATION_COLORS = {
    'Orion':       '#4fc3f7',   # sky blue
    'Canis Major': '#fff176',   # yellow
    'Carina':      '#ce93d8',   # lavender
    'Taurus':      '#ffb74d',   # orange
    'Lyra':        '#80deea',   # cyan
    'Draco':       '#a5d6a7',   # green
    'Ursa Minor':  '#ef9a9a',   # pink
    'Cygnus':      '#b39ddb',   # purple
    'Boötes':      '#ff8a65',   # deep orange
    'Leo':         '#f48fb1',   # rose
    'Virgo':       '#c5e1a5',   # light green
    'Scorpius':    '#ff5252',   # red
}

DEFAULT_COLOR = '#e0e0e0'


def _magnitude_to_size(mag):
    """Convert apparent magnitude to a plot marker size (brighter = bigger)."""
    return max(20, (4.0 - mag) ** 2 * 8)


def plot_sky_map(results, monument_name=None, threshold_deg=2.0, save_path=None):
    """
    Draw a polar horizon chart showing star positions.

    Polar convention:
      - 0° at top  = North
      - Clockwise  = East
      - r = 0      = Zenith  (90° altitude)
      - r = 90     = Horizon (0° altitude)
    Only stars above the horizon (altitude > 0) are plotted.

    Parameters
    ----------
    results       : output of calculate_alignments()
    monument_name : if given, draw the monument's orientation line and flag matches
    threshold_deg : azimuth tolerance used for alignment matching
    save_path     : if given, save figure to this path instead of showing it
    """
    fig = plt.figure(figsize=(10, 10), facecolor='#0d1117')
    ax = fig.add_subplot(111, projection='polar')

    # --- polar axis cosmetics ---
    ax.set_facecolor('#0d1117')
    ax.set_theta_zero_location('N')   # 0° at top
    ax.set_theta_direction(-1)        # clockwise

    # radial axis: r=0 is zenith, r=90 is horizon
    ax.set_ylim(0, 90)
    ax.set_yticks([0, 30, 60, 90])
    ax.set_yticklabels(['Zenith', '60°', '30°', 'Horizon'],
                       color='#8b949e', fontsize=8)
    ax.yaxis.set_tick_params(labelcolor='#8b949e')

    # cardinal direction labels
    ax.set_xticks(np.radians([0, 45, 90, 135, 180, 225, 270, 315]))
    ax.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'],
                       color='#c9d1d9', fontsize=10, fontweight='bold')
    ax.grid(color='#21262d', linewidth=0.5, linestyle='--', alpha=0.6)
    for spine in ax.spines.values():
        spine.set_color('#21262d')

    # --- alignment-matched stars (precompute for highlighting) ---
    aligned_stars = set()
    if monument_name:
        orientation_az = MONUMENTS[monument_name]['orientation_az']
        matches = check_alignments(results, orientation_az, threshold_deg)
        aligned_stars = {m['star'] for m in matches}

        # Draw the monument orientation sightline
        theta = np.radians(orientation_az)
        ax.plot([theta, theta], [0, 90], color='#ffd700', linewidth=2,
                linestyle='--', alpha=0.85, zorder=2)
        ax.plot([theta, theta], [0, 90], color='#ffd700', linewidth=0.5,
                linestyle='--', alpha=0.3, zorder=1)

        # opposing direction (180° away)
        theta_opp = np.radians((orientation_az + 180) % 360)
        ax.plot([theta_opp, theta_opp], [0, 90], color='#ffd700', linewidth=1,
                linestyle=':', alpha=0.4, zorder=1)

    # --- plot each visible star ---
    legend_consts = {}
    for name, data in results['stars'].items():
        if not data['visible']:
            continue

        theta = np.radians(data['azimuth'])
        r     = 90 - data['altitude']          # zenith at center
        mag   = STARS[name]['mag']
        const = STARS[name]['constellation']
        color = CONSTELLATION_COLORS.get(const, DEFAULT_COLOR)
        size  = _magnitude_to_size(mag)

        if name in aligned_stars:
            # Highlight aligned stars with a gold ring
            ax.scatter(theta, r, s=size * 2.5, color='#ffd700',
                       alpha=0.35, zorder=3, linewidths=0)
            ax.scatter(theta, r, s=size, color=color,
                       edgecolors='#ffd700', linewidths=1.5, zorder=4)
        else:
            ax.scatter(theta, r, s=size, color=color,
                       alpha=0.90, zorder=3, edgecolors='none')

        # Label — offset slightly so it doesn't overlap the dot
        label_r = r - 3 if r > 5 else r + 4
        ax.annotate(
            name,
            xy=(theta, r),
            xytext=(theta, label_r),
            color='#c9d1d9',
            fontsize=7,
            ha='center', va='center',
            zorder=5,
        )

        if const not in legend_consts:
            legend_consts[const] = color

    # --- legend ---
    patches = [
        mpatches.Patch(color=c, label=cn)
        for cn, c in sorted(legend_consts.items())
    ]
    if monument_name:
        patches.append(
            mpatches.Patch(color='#ffd700',
                           label=f'Orientation: {MONUMENTS[monument_name]["orientation_az"]}°')
        )
    ax.legend(
        handles=patches,
        loc='lower left',
        bbox_to_anchor=(-0.12, -0.05),
        framealpha=0.15,
        labelcolor='#c9d1d9',
        fontsize=8,
    )

    # --- title ---
    era  = "BC" if results['_year'] < 0 else "AD"
    date = f"{abs(results['_year'])} {era}-{results['_month']:02d}-{results['_day']:02d}  {results['_hour']:.1f}h"
    loc  = f"({results['_lat']:.2f}°, {results['_lon']:.2f}°)"
    title_parts = [date, loc]
    if monument_name:
        title_parts.insert(0, monument_name)
    ax.set_title(
        '\n'.join(title_parts),
        color='#c9d1d9', fontsize=11, pad=18,
    )

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
        print(f"Saved: {save_path}")
    else:
        plt.show()
    plt.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Archaeo-Astronomy Sky Map Visualizer",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python src/visualize.py --monument Stonehenge --year -2500\n"
            "  python src/visualize.py --monument Giza --year -10499 --hour 22\n"
            "  python src/visualize.py --lat 51.17 --lon -1.82 --year -2500 --save outputs/map.png\n"
        ),
    )
    parser.add_argument("--lat",    type=float, default=29.9792)
    parser.add_argument("--lon",    type=float, default=31.1342)
    parser.add_argument("--year",   type=int,   default=-2500)
    parser.add_argument("--month",  type=int,   default=3)
    parser.add_argument("--day",    type=int,   default=20)
    parser.add_argument("--hour",   type=float, default=22.0,
                        help="Hour in 24h format (default: 22.0 = 10 PM)")
    parser.add_argument("--monument", type=str, default=None,
                        help="Monument name (overrides --lat/--lon, adds orientation line)")
    parser.add_argument("--threshold", type=float, default=2.0,
                        help="Azimuth tolerance in degrees for alignment highlights (default: 2.0)")
    parser.add_argument("--save",   type=str, default=None,
                        help="Save chart to this file path (PNG/PDF/SVG) instead of displaying it")

    args = parser.parse_args()

    lat, lon = args.lat, args.lon
    monument_name = None

    if args.monument:
        key = args.monument.lower()
        found = [k for k in MONUMENTS if key in k.lower()]
        if not found:
            print(f"Monument '{args.monument}' not found.")
            raise SystemExit(1)
        monument_name = found[0]
        lat = MONUMENTS[monument_name]['lat']
        lon = MONUMENTS[monument_name]['lon']
        print(f"Monument: {monument_name}  ({lat}, {lon})")

    results = calculate_alignments(lat, lon, args.year, args.month, args.day, args.hour)
    results.update({
        '_lat': lat,   '_lon': lon,
        '_year': args.year, '_month': args.month,
        '_day': args.day,   '_hour': args.hour,
    })

    plot_sky_map(results, monument_name=monument_name,
                 threshold_deg=args.threshold, save_path=args.save)
