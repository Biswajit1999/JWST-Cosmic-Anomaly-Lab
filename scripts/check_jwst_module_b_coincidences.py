"""Measure incidental cross-channel matches after fixed short-wave sky shifts."""

import json

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from download_jwst_cross_detector import OUT
from search_jwst_cross_detector import pairs

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json

RADII = (3., 6., 9.)
ANGLES = np.arange(0, 360, 45)


def counts(peaks):
    result = pairs(peaks)
    if result.empty:
        return dict(pairs=0, two_filter=0, three_filter=0)
    confirmation = result.get("f356w_aperture_gate", pd.Series(False, index=result.index))
    return dict(pairs=len(result), two_filter=int(result.aperture_gate.sum()),
                three_filter=int((result.aperture_gate & confirmation.eq(True)).sum()))


def shifted_trials(peaks):
    choose = peaks.band == "F200W"
    source = SkyCoord(peaks.loc[choose, "ra"], peaks.loc[choose, "dec"], unit="deg")
    trial = peaks.copy()
    results = []
    for radius in RADII:
        for angle in ANGLES:
            moved = source.directional_offset_by(angle*u.deg, radius*u.arcsec)
            trial.loc[choose, "ra"] = moved.ra.deg
            trial.loc[choose, "dec"] = moved.dec.deg
            results.append(dict(radius_arcsec=radius, position_angle_deg=int(angle), **counts(trial)))
    return results


def main():
    peaks_path = OUT/"module_b_peaks.csv"
    search_path = OUT/"module_b_search.json"
    search = json.loads(search_path.read_text())
    peaks = pd.read_csv(peaks_path, dtype={"pointing": str})
    results = {}
    for pointing, group in peaks.groupby("pointing"):
        observed = counts(group)
        expected = search["pointings"][pointing]
        if tuple(observed.values()) != tuple(expected[key] for key in (
                "cross_channel_pairs", "passes_pair_aperture_gate", "passes_three_filter_gate")):
            raise ValueError(f"Unshifted counts differ from search ledger for {pointing}")
        results[pointing] = dict(observed=observed, shifted=shifted_trials(group))
        print(pointing, observed, "shifted max", {key: max(item[key] for item in results[pointing]["shifted"])
                                                  for key in observed}, flush=True)
    write_json(OUT/"module_b_coincidences.json", {
        "executed_utc": utcnow(), "pointings": results, "peaks_sha256": text_digest(peaks_path),
        "search_sha256": text_digest(search_path), "script_sha256": text_digest(__file__),
        "scope": "Incidental positional matches; shifts break co-spatial artifacts and can move edge peaks off the common footprint",
    })


if __name__ == "__main__":
    main()
