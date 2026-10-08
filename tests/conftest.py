import sys
from pathlib import Path

from astropy.utils import iers

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT))

# astropy is only a reference near J2000, which the Earth orientation data bundled
# with astropy covers; never reach for the network during the tests.
iers.conf.auto_download = False
