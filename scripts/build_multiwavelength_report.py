"""Build a compact comparison figure and report from executed results."""

import json
from pathlib import Path

import pandas as pd
from matplotlib.figure import Figure

from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/multiwavelength"
ratios = pd.read_csv(out / "dja_core_ratios.csv")
fig = Figure(figsize=(10, 4.6), layout="constrained")
for ax, (number, line, ref) in zip(fig.subplots(1, 3),
                                  [(66336, "ha", 2.86), (10009106, "hg", .468), (89464, "hg", .468)], strict=True):
    sub = ratios[ratios.dr4_id.str.endswith("_" + str(number)) & (ratios.line == line)]
    for i, (reduction, color) in enumerate([("DJA", "#126e82"), ("DR4 5px/3nod", "#b94d36"), ("DR4 5px/2nod", "#9b8b4e")]):
        for half, offset, marker in [(400, -.09, "o"), (600, .09, "s")]:
            r = sub[(sub.reduction == reduction) & (sub.half_width_kms == half)].iloc[0]
            if r.eligible:
                ax.errorbar(i+offset, r.ratio, yerr=r.error, fmt=marker, color=color, ms=5, capsize=3)
    ax.axhline(ref, color="#666666", ls="--", lw=1)
    ax.set(xticks=[0, 1, 2], xticklabels=["DJA", "DR4\n3 nods", "DR4\n2 nods"],
           ylabel="Hα/Hβ" if line == "ha" else "Hγ/Hβ", title=str(number), xlim=(-.5, 2.5))
fig.suptitle("Observed line-core ratios across reductions\nCircles: ±400 km/s; squares: ±600 km/s; conditional diagonal errors", fontsize=12)
save_figure(fig, root / "figures/multiwavelength/reduction_ratio_comparison")
context = json.loads((out / "context_summary.json").read_text())
text = """# Cross-reduction and multi-telescope validation

Research snapshot: **16 November 2026**. Author: Biswajit Jana.

## Result that changes the shortlist

Six additional released DJA medium-resolution spectra were checked for four DR4 targets.
Source number, programme, catalogue coordinates and FITS coordinates/grating were checked.
The same signed-core algorithm, fixed DR4 redshift and 0.6–5.5 µm mask were used.
There are 72 DJA/DR4 core attempts and 60 usable measurements. Reductions share
observations, and this is not replication with independent photons.

| Target | Diagnostic | DJA ±400 / ±600 km/s | DR4 5px, 3 nods ±400 / ±600 km/s |
|---|---|---|---|
"""
for number, line in [(66336, "ha"), (10009106, "hg"), (89464, "hg")]:
    sub = ratios[ratios.dr4_id.str.endswith("_" + str(number)) & (ratios.line == line)]
    values = []
    for reduction in ["DJA", "DR4 5px/3nod"]:
        values.append(" / ".join(f"{r.ratio:.3f} ± {r.error:.3f}" for r in sub[sub.reduction == reduction].itertuples()))
    text += f"| {number} | {line.upper()}/HB | {values[0]} | {values[1]} |\n"
text += """
The extreme ratios for GN-66336 and GS-10009106 are not reproduced by DJA. The
different absolute line fluxes must be investigated before physical interpretation.
GN-89464 retains an elevated Hγ/Hβ core ratio in the tested reductions; this remains
a follow-up diagnostic, not a black-hole-star identification. GN-25734 has no usable
DJA Hα core under the fixed coverage rule and is retained as an extraction-sensitive control.
The errors above omit cross-pixel and cross-reduction covariance. No independent
significance is assigned to reduction differences. Reference ratios 2.86 and 0.468
in the figure are illustrative nebular values, not universal physical limits.

![Reduction test](../figures/multiwavelength/reduction_ratio_comparison.png)

![Absolute spectral comparison](../figures/multiwavelength/goods-s-mediumjwst_10009106_reductions.png)

## Hubble, Spitzer and Chandra context

The GOODS-N/S 3D-HST v4.1 catalogues contain 38,279 and 50,507 rows. We searched
them for all 12 queued targets, retaining every candidate within 0.3 arcsec without
merging identities. Eight targets each have one positional candidate. Four do not.
These 88,786 catalogue rows were searched for neighbours, not individually classified.

The eight candidates supply 88 band rows: 74 have valid signed flux/error/weight
values, including 17 negative fluxes. Applying use_phot=1 leaves 65 displayable
measurements. Fluxes use the documented AB=25 zeropoint: one unit is 0.3631 µJy.
Nonpositive errors/weights, nonfinite values and the -99 sentinel are missing.
Signed low-S/N values are retained, not converted into positive detections or limits.
Both plotted low-redshift candidates carry contam_flag=1 and extraction flags=3;
their IRAC photometry requires deblending review. The GS-10009106 candidate is
0.294 arcsec away, near our search boundary. No aperture matching or SED fit is claimed.

![Published photometry](../figures/multiwavelength/hst_spitzer_context.png)

The CDS Chandra endpoints timed out; NASA HEASARC supplied the 683-record CDF-N
and 1,008-record CDF-S main catalogues through TAP. An initial unsupported CSV-format
response was rejected by schema/row-count checks; VOTable retrieval succeeded.
None of the 12 target positions has a main-catalogue neighbour within 2 arcsec.
This does **not** measure a source-specific X-ray upper limit or rule out an AGN.
No event files, local backgrounds, sensitivity maps or exposure maps were analysed.

## The paper behind the screenshot

The MIT headline dated 12 August 2026 corresponds to Naidu et al.,
[A gas-enshrouded and gas-reddened black hole at cosmic dawn](https://doi.org/10.1038/s41586-026-10846-4).
It concerns MoM-BH*-1, not one of our newly ranked targets. Its combination of
continuum and resolved Balmer features supports the authors' gas-enshrouded
accretion interpretation. We have not independently reproduced its physical fit.
The [GLIMPSE-17775 study](https://arxiv.org/abs/2511.07515) is a separate benchmark.

## Reproduction and limitations

Run the five scripts listed in `experiments/multiwavelength_validation.yaml`.
Raw URLs, SHA256 values, retrieval failures and measurements are retained in
`results/multiwavelength/` and `data/manifests/downloads.json`. The original full
DR4 census remains unchanged as a historical result; this experiment updates its
interpretation. No new source count, discovery, or all-archive coverage is asserted.
"""
(root / "reports/multiwavelength_validation.md").write_text(text, encoding="utf-8")
