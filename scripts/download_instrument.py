import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cosmic_anomaly_lab.provenance import DownloadCache

cache = DownloadCache(Path(__file__).resolve().parents[1])
cache.fetch(
    "https://raw.githubusercontent.com/gbrammer/msaexp/main/msaexp/data/jwst_nirspec_prism_disp.fits",
    "instrument/prism_disp.fits",
    survey="instrument-reference",
    version="msaexp-main-snapshot",
    kind="resolving_power",
)
