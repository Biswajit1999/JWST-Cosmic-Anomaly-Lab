"""Targeted post-selection check of published JADES-GS-159717, not a blind detection."""

from pathlib import Path

import pandas as pd
from matplotlib.figure import Figure

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.resolution_check import halpha_core_contrast
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
source = "CAL-82610407f3945338"
catalogue = pd.read_csv(root / "data/catalogues/observations.csv.gz")
rows = catalogue[(catalogue.source_id == source) & catalogue.grating.isin(["G395M", "G395H"])]
cache = DownloadCache(root, budget_bytes=4096 * 1024**2)
fig = Figure(figsize=(11, 5), layout="constrained")
axes = fig.subplots(1, 2)
records = []
for row in rows.itertuples():
    path = cache.fetch(
        row.url,
        "spectra/" + row.file,
        survey=row.survey,
        version="DJA-v4.4",
        kind="observed_spectrum",
    )
    spec = read_spectrum(path)
    diagnostics = [
        dict(redshift_basis=basis, **halpha_core_contrast(spec, z))
        for basis, z in [("archive", row.z_best), ("published_systemic", 5.07781)]
    ]
    records.append(
        dict(file=row.file, grating=row.grating, raw_sha256=digest(path), diagnostics=diagnostics)
    )
    velocity = (spec.wavelength.value / (0.656461 * (1 + 5.07781)) - 1) * 299792.458
    color = "#a84f2e" if row.grating == "G395H" else "#287b97"
    for ax, limit in zip(axes, [2000, 650], strict=True):
        mask = spec.good & (abs(velocity) < limit)
        ax.errorbar(
            velocity[mask],
            spec.flux.value[mask],
            yerr=spec.error.value[mask],
            fmt=".-",
            ms=3,
            lw=0.8,
            color=color,
            alpha=0.85,
            label=row.grating,
        )
        ax.set(
            xlim=(-limit, limit),
            xlabel="Velocity from Hα at published z=5.07781 (km/s)",
            ylabel="Released flux density (µJy)",
        )
        ax.axhline(0, color="grey", lw=0.6)
        ax.grid(alpha=0.15)
        ax.legend(fontsize=9)
for low, high in [(-450, -250), (250, 450)]:
    axes[1].axvspan(low, high, color="#d8dfe5", alpha=0.4)
axes[1].axvspan(-120, 120, color="#ead3c0", alpha=0.3)
axes[0].set_title("Observed profiles · same archive reduction")
axes[1].set_title("Fixed core and shoulder windows")
fig.suptitle(
    "JADES-GS-159717 · recovered known LRD\n"
    "Resolution-dependent local contrast; no new physical classification",
    fontsize=13,
)
save_figure(fig, root / "figures/recovered_lrd_resolution")
write_json(
    root / "results/open_search/recovered_lrd_resolution.json",
    dict(
        executed_utc=utcnow(),
        source_id=source,
        spectra=records,
        literature="https://arxiv.org/abs/2506.14870v2",
        design="Targeted after literature lookup; signed mean shoulder minus core contrast, "
        "fixed windows and two specified redshifts, no search significance",
        limitations=[
            "Same DJA reduction; not independent extraction",
            "No full emission+absorption+LSF physical model",
            "Positive contrast is a local dip, not identification of its cause",
            "Quoted errors omit resampling/calibration systematics",
        ],
        code_sha256=text_digest(root / "src/cosmic_anomaly_lab/resolution_check.py"),
        script_sha256=text_digest(Path(__file__)),
    ),
)
