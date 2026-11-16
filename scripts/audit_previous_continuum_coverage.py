"""Audit all previous 288 queue measurements without altering the dated outputs."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u

from cosmic_anomaly_lab.continuum_audit import bilateral_core
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum

root = Path(__file__).resolve().parents[1]
out = root / "results/gn89464_followup"
source = root / "results/jades_balmer_census/signed_core_measurements.csv"
old = pd.read_csv(source)
products = json.loads((root / "results/jades_balmer_census/spectral_products.json").read_text())
accepted = {(p["dr4_id"], p["grating_filter"]): p for p in products if p["status"] == "DOWNLOADED"}
hashes = {}
for p in accepted.values():
    assert digest(root / p["path"]) == p["sha256"]
    hashes[p["path"]] = p["sha256"]
spectra, rows = {}, []
unit = u.erg / u.s / u.cm**2 / u.AA
for r in old.itertuples():
    result = dict(dr4_id=r.dr4_id, line=r.line, aperture_pixels=r.aperture_pixels, nods=r.nods,
                  grating_filter=r.grating_filter, half_width_kms=r.half_width_kms,
                  legacy_status=r.status, legacy_flux=r.flux, legacy_error=r.error)
    p = accepted.get((r.dr4_id, r.grating_filter))
    if p is None:
        result.update(status="NO_ACCEPTED_PRODUCT", flux=None, error=None)
    else:
        key = (p["path"], r.aperture_pixels, r.nods)
        if key not in spectra:
            s = read_dr4_spectrum(root / p["path"], r.aperture_pixels, r.nods, wavelength_range=(.6, 5.5))
            f = s.flux.to_value(unit, equivalencies=u.spectral_density(s.wavelength)) / 1e-20
            e = s.error.to_value(unit, equivalencies=u.spectral_density(s.wavelength)) / 1e-20
            spectra[key] = s, f, e
        s, f, e = spectra[key]
        v = (s.wavelength.value/r.center_um-1)*299792.458
        result.update(bilateral_core(v, f, e, s.good, r.center_um*1e4, r.half_width_kms))
        if result["status"] == "MEASURED_CORE":
            assert r.status == "MEASURED_CORE"
            np.testing.assert_allclose([result["flux"], result["error"]], [r.flux, r.error], rtol=1e-8, atol=1e-8)
    result["newly_rejected"] = r.status == "MEASURED_CORE" and result["status"] != "MEASURED_CORE"
    rows.append(result)
data = pd.DataFrame(rows)
data.to_csv(out / "previous_queue_coverage_audit.csv", index=False)
summary = dict(executed_utc=utcnow(), attempts=len(data), legacy_measured=int((old.status == "MEASURED_CORE").sum()),
               bilateral_measured=int((data.status == "MEASURED_CORE").sum()),
               newly_rejected=int(data.newly_rejected.sum()),
               affected_targets=sorted(data.loc[data.newly_rejected, "dr4_id"].unique().tolist()),
               retained_measurements_match_legacy=True, input_sha256=hashes,
               old_measurements_sha256=text_digest(source), script_sha256=text_digest(Path(__file__)),
               module_sha256=text_digest(root / "src/cosmic_anomaly_lab/continuum_audit.py"))
write_json(out / "previous_queue_summary.json", summary)
print({k: v for k, v in summary.items() if not k.endswith("sha256")})
print(data[data.newly_rejected].to_string(index=False))
