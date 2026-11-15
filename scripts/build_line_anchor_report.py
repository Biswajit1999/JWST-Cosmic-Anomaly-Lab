"""Summarize line constraints and inventory broader existing follow-up priorities."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/line_anchor"
summary = json.loads((out / "summary.json").read_text())
widths = pd.read_csv(out / "width_constraints.csv")
assert summary["attempted_windows"] == 22 and summary["fitted_windows"] == 12
selected = widths[(widths.half_window_kms == 900) & (widths.status == "CONDITIONAL_GAUSSIAN_WIDTH")]
fig = Figure(figsize=(11, 5), layout="constrained")
axes = fig.subplots(1, 2)
for ax, label in zip(axes, ["DJA G395H", "DJA G235M"], strict=True):
    row = selected[selected.label == label].iloc[0]
    q = np.linspace(0.8, 3, 500)
    instrumental_sigma = 299792.458 / (2.354820045 * row.reference_R_at_oiii * q)

    def intrinsic(observed):
        variance = observed**2 - instrumental_sigma**2
        return np.sqrt(np.where(variance >= 0, variance, np.nan))

    ax.plot(q, intrinsic(row.observed_sigma_kms), color="#287a8c", label="Best observed width")
    ax.fill_between(
        q,
        np.sqrt(np.maximum(0, row.sigma_delta3p84_lower**2 - instrumental_sigma**2)),
        intrinsic(row.sigma_delta3p84_upper),
        color="#287a8c",
        alpha=0.2,
        label="Formal observed-width range",
    )
    ax.axvline(row.conditional_q_lower, color="#a25133", ls="--", label="Conditional lower q bound")
    ax.set(
        title=label,
        xlabel="Assumed R_eff / reference R at [OIII]",
        ylabel="Compatible intrinsic gas σ (km/s)",
        ylim=(0, 100),
    )
    ax.grid(alpha=0.15)
    ax.legend(fontsize=8)
fig.suptitle(
    "CONDITIONAL decomposition · many gas widths and instrument responses fit one observed width\nGaussian variance addition; curves are not measurements of the instrument or gas velocity",
    fontsize=11,
)
save_figure(fig, root / "figures/line_width_degeneracy")

# Archive listings for existing priorities: no automatic source matching or new detections.
old = pd.read_csv(root / "results/expanded/followup_priority.csv")
new = pd.read_csv(root / "results/open_search/shape_followup.csv")
priorities = pd.concat(
    [old.assign(queue="joint_screen"), new.assign(queue="shape_search")], ignore_index=True
)
assert priorities.source_id.is_unique
observations = pd.read_csv(root / "data/catalogues/observations.csv.gz")
inventory = []
products = []
for source in priorities.itertuples():
    obs = observations[observations.source_id == source.source_id]
    gratings = sorted(set(obs.grating.dropna()))
    high = [g for g in gratings if g in ["G140H", "G235H", "G395H"]]
    medium = [g for g in gratings if g in ["G140M", "G235M", "G395M"]]
    inventory.append(
        dict(
            source_id=source.source_id,
            queue=source.queue,
            survey=source.survey,
            z_best=source.z_best,
            halpha_observed_um=0.656461 * (1 + source.z_best),
            oiii_observed_um=0.500824 * (1 + source.z_best),
            listed_high_gratings=";".join(high),
            listed_medium_gratings=";".join(medium),
            listed_spectra=len(obs),
            has_listed_high_grating=bool(high),
            has_listed_medium_grating=bool(medium),
            usable_line_coverage="NOT_ESTABLISHED_BY_CATALOGUE",
        )
    )
    for product in obs.itertuples():
        products.append(
            dict(
                source_id=source.source_id,
                file=product.file,
                grating=product.grating,
                url=product.url,
                cached=(root / "data/raw/spectra" / product.file).is_file(),
            )
        )
inventory = pd.DataFrame(inventory)
inventory.to_csv(out / "priority_followup_inventory.csv", index=False)
pd.DataFrame(products).to_csv(out / "priority_archive_products.csv", index=False)
write_json(
    out / "inventory_provenance.json",
    dict(
        executed_utc=utcnow(),
        sources=len(inventory),
        with_listed_high_grating=int(inventory.has_listed_high_grating.sum()),
        with_listed_medium_grating=int(inventory.has_listed_medium_grating.sum()),
        interpretation="Archive product inventory only; not inspected coverage or new objects",
        input_sha256={
            p.relative_to(root).as_posix(): (digest(p) if p.suffix == ".gz" else text_digest(p))
            for p in [
                root / "results/expanded/followup_priority.csv",
                root / "results/open_search/shape_followup.csv",
                root / "data/catalogues/observations.csv.gz",
            ]
        },
        script_sha256=text_digest(Path(__file__)),
    ),
)
table = "\n".join(
    f"| {r.label} | {r.observed_sigma_kms:.2f} | {r.sigma_delta3p84_lower:.2f}–{r.sigma_delta3p84_upper:.2f} | {r.reduced_chi2:.2f} | {r.conditional_R_lower_at_oiii:.0f} |"
    for r in selected.itertuples()
)
h = selected[selected.label == "DJA G395H"].iloc[0]
medium = selected[selected.label == "DJA G235M"].iloc[0]
dr4 = selected[selected.label == "DR4 G235M 3px/3-nod"].iloc[0]
report = f"""# Can another line calibrate the instrumental blur?

Research snapshot: **15 November 2026**, following the [Balmer benchmark](balmer_audit.md).
Author: Biswajit Jana. Observed measurements and explicitly conditional models.

**The [OIII] doublet provides a conditional lower limit on resolving power, not a
unique instrument calibration.** The same observed width can be produced by
different combinations of intrinsic gas motion and instrument blur. G235M also
shows a reduction-dependent width difference. We therefore do not use this test
to claim a uniquely calibrated absorber or a new black-hole star.

## Data and coverage

We checked eleven extraction variants: DJA G395H, G235M and G395M, plus four
DR4 aperture/nod variants for each medium-resolution grating. Both [OIII] lines
are fitted in ±500 and ±900 km/s windows, giving **22 attempted measurements,
12 fitted and 10 explicitly unavailable**. The unavailable measurements are all
G395M: the released data contain zero usable pixels in these line windows.
A grating name alone would have incorrectly implied a usable calibration line.

Two new observed products (DJA and JADES DR4 G235M), an actual DR4 directory
listing, and the STScI G395H resolution table were retrieved through the finite
4 GiB cache. No archive neighbour was substituted. DR4 target, programme and
coordinates were checked against the existing GS-159717 association. All paths
and hashes are retained in the [input record](../results/line_anchor/inputs.json)
and [download record](../results/line_anchor/downloads.json).

## Observed narrow-line measurements

The model has two independently fitted nonnegative line areas, a common observed
Gaussian width and velocity, and a linear continuum. The area ratio is free.
We integrate the model over the spectral pixels, retain signed flux and missing
values, and use the released diagonal errors without rescaling or residual
clipping. No Hα pixels enter these fits. Wavelengths are vacuum 4960.295 and
5008.240 Å, with z=5.07781 as the velocity reference and a free fitted shift.

Results for ±900 km/s windows:

| Extraction | Observed σ (km/s) | Formal Δχ²≤3.84 range | χ²/(n−6) | Conditional lower R at [OIII] |
| --- | ---: | ---: | ---: | ---: |
{table}

The profile statistic refits continuum, amplitudes and velocity at fixed width.
Its threshold is a formal one-parameter likelihood convention, **not a calibrated
95% confidence statement** when correlated errors, non-Gaussian profiles or
calibration systematics are omitted. Grid limits that fail to bracket a crossing
produce null endpoints, not invented limits. All profiles and both window sizes
are stored in [fits](../results/line_anchor/fits.json) and the
[measurement table](../results/line_anchor/width_constraints.csv).

![Observed doublet profiles](../figures/oiii_line_anchor.png)

The high-resolution width is stable between the two windows (60.14 versus
60.22 km/s), but its reduced χ² is 2.59 in the narrower window and 1.77 in the
wider one. The Gaussian model and diagonal errors do not fully describe the
data; stable widths alone do not certify the model. We retain that limitation
when discussing instrument bounds.

DJA G235M gives {medium.observed_sigma_kms:.2f} km/s; DR4 3px/3-nod gives
{dr4.observed_sigma_kms:.2f} km/s, about {(dr4.observed_sigma_kms / medium.observed_sigma_kms - 1) * 100:.1f}% broader.
This is a comparison of reductions/extractions with shared observations, not an
independent difference-significance test. We do not attribute the discrepancy
uniquely to a pipeline defect, physical variability or gas kinematics.

![Conditional width profiles](../figures/oiii_width_profiles.png)

## What the width can and cannot establish

Under a local Gaussian convolution model,

`sigma_observed² = sigma_intrinsic² + sigma_instrument²`.

Nonnegative intrinsic variance implies
`R_effective ≥ c / (2.35482 * sigma_observed_upper)`.
It supplies **no finite upper limit on R**, because a narrower instrument profile
can be offset by greater intrinsic gas width. Neither term should be fixed to
zero merely to force a unique answer.

![Instrument/gas degeneracy](../figures/line_width_degeneracy.png)

For G395H the formal [OIII] bound is R≥{h.conditional_R_lower_at_oiii:.0f}.
The [STScI instrument documentation](https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-instrumentation/nirspec-dispersers-and-filters)
provides wavelength-dependent reference resolving powers. Its tabulated values
at this source's [OIII] and Hα wavelengths are {h.reference_R_at_oiii:.0f} and
{h.reference_R_at_halpha:.0f}. **Only if** the ratio of effective to reference
resolution is constant within G395H can the lower bound transfer to approximately
R(Hα)≥{h.conditional_same_grating_R_lower_at_halpha:.0f}. This assumption does not
determine the source's illumination profile or guarantee that [OIII] and Balmer
components have the same effective response. Model mismatch further limits the
formal bound. It is not a hard rejection threshold for the previous Hα fits.

G235M cannot be extrapolated to Hα in G395M. The gratings and wavelength coverage
differ, and the common scale factor would be an additional untested assumption.
We leave that transfer null. The medium-resolution absorption sensitivity
identified previously therefore remains unresolved.

The [published source study](https://arxiv.org/html/2506.14870v2) also discusses
line-spread-function and spatial-extent complications. This experiment measures
our released products; it does not reproduce its full multi-line Bayesian model.

## Returning to the wider candidate set

We also inventoried all **{len(inventory)} existing priorities** from the two
search queues. {int(inventory.has_listed_high_grating.sum())} have a listed
high-resolution grating product and {int(inventory.has_listed_medium_grating.sum())}
have a listed medium-resolution product; these counts can overlap. Product
listings are not evidence of usable line coverage. The GS-159717 G395M gap
demonstrates why coverage must be checked from pixels before further fitting.

The [priority inventory](../results/line_anchor/priority_followup_inventory.csv)
and [individual archive products](../results/line_anchor/priority_archive_products.csv)
record exact source IDs, URLs, cached-file status and expected observed line
wavelengths. No new source has been claimed or automatically merged. The total
search remains 4,731 analysed records.

The [subsequent high-resolution follow-up](high_resolution_followup.md) now
executes pixel-level checks on all twelve listed high-resolution products and
conditional Hα model comparisons where coverage permits. The known benchmark
exercises expose limitations that should be resolved before applying discovery
language to unfamiliar objects.

## Reproduction

```sh
python scripts/replay_manifest.py --budget-mib 4096
python scripts/download_line_anchor.py
python scripts/analyse_line_anchor.py
python scripts/build_line_anchor_report.py
pytest -q
ruff check src tests scripts
```

See the [run summary](../results/line_anchor/summary.json),
[inventory provenance](../results/line_anchor/inventory_provenance.json) and
[execution ledger](../experiments/line_anchor.yaml). Figures have PNG and PDF
versions. Conference posters remain the dated 12 November snapshot until the
final synthesis; this report adds reproducible measurements for that revision.
"""
(root / "reports/line_anchor.md").write_text(report, encoding="utf-8")
print(
    inventory[["source_id", "listed_high_gratings", "listed_medium_gratings"]].to_string(
        index=False
    )
)
