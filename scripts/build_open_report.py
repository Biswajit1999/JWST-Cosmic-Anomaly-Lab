"""Build an auditable report and overview from the completed open-search run."""

import json
import platform
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.cosmology import Planck18
from matplotlib.figure import Figure

from cosmic_anomaly_lab.provenance import digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/open_search"
write_json(
    out / "software_versions.json",
    dict(
        recorded_utc=utcnow(),
        python=platform.python_version(),
        platform=platform.system(),
        packages={
            name: version(name)
            for name in [
                "numpy",
                "scipy",
                "pandas",
                "astropy",
                "matplotlib",
                "pyyaml",
                "requests",
            ]
        },
        prism_resolution_reference_sha256=digest(root / "data/raw/instrument/prism_disp.fits"),
    ),
)
summary = json.loads((out / "summary.json").read_text())
inspection = json.loads((out / "inspection_audit.json").read_text())["sources"]
validation = json.loads((out / "synthetic_validation.json").read_text())
resolution_check = json.loads((out / "recovered_lrd_resolution.json").read_text())
scores = pd.read_csv(out / "shape_scores.csv")
queue = pd.read_csv(out / "shape_followup.csv")
new = pd.read_csv(out / "new_measurements.csv")
measurements = pd.read_csv(out / "combined_measurements.csv")
selection = pd.concat(
    [
        pd.read_csv(root / "data/catalogues/expanded_selection.csv.gz"),
        pd.read_csv(root / "data/catalogues/open_selection.csv.gz"),
    ]
)
notes_path = out / "review_notes.json"
notes = json.loads(notes_path.read_text()) if notes_path.exists() else {}
fig = Figure(figsize=(12, 9), layout="constrained")
axes = fig.subplots(2, 2)
colors = dict(CEERS="#286b9b", JADES="#4c8d75", RUBIES="#b46e3e", UNCOVER="#845f9c")
for ax, mode in zip(axes[0], ["optical", "ultraviolet"], strict=True):
    for survey, group in scores[scores["mode"] == mode].groupby("survey"):
        ax.scatter(
            group.z_best,
            group.perturbation_q10,
            s=8,
            alpha=0.45,
            color=colors[survey],
            label=survey,
        )
    selected = queue[queue["mode"] == mode]
    ax.scatter(
        selected.z_best,
        selected.perturbation_q10,
        s=70,
        facecolors="none",
        edgecolors="black",
        label="Inspection queue",
    )
    ax.set(
        xlabel="Archive redshift",
        ylabel="Perturbed-score 10th percentile",
        title=f"{mode.capitalize()} shapes · descriptive rarity",
    )
    ax.legend(fontsize=8)
    ax.grid(alpha=0.15)
for label, subset, color in [
    ("Previous: archive S/N ≥ 5", selection[selection.sn50 >= 5], "#286b9b"),
    ("Added: 2 ≤ archive S/N < 5", selection[selection.sn50 < 5], "#b46e3e"),
]:
    axes[1, 0].hist(
        Planck18.age(subset.z_best.to_numpy()).value,
        bins=np.arange(0, 9.1, 0.5),
        histtype="step",
        lw=1.7,
        color=color,
        label=label,
    )
axes[1, 0].set(
    xlabel="Cosmic age at observed epoch (Gyr; Planck18)",
    ylabel="Selected source records",
    title="Sample counts, not cosmic abundance",
)
axes[1, 0].legend(fontsize=8)
axes[1, 1].axis("off")
status_text = "COVERAGE AUDIT\n\n" + "\n".join(
    f"{mode}: "
    + ", ".join(
        f"{r['count']:,} {r['status'].lower().replace('_', ' ')}"
        for r in summary["shape_status_counts"]
        if r["mode"] == mode
    ).replace(", ", "\n  ")
    for mode in ["optical", "ultraviolet"]
)
axes[1, 1].text(
    0.02,
    0.97,
    status_text + "\n\nRarity is not a discovery probability.",
    va="top",
    fontsize=9,
    linespacing=1.4,
)
fig.suptitle(
    f"OBSERVED expanded search · {summary['combined_analysed']:,} analysed source records",
    fontsize=14,
)
save_figure(fig, root / "figures/open_search_overview")
lines = [
    "# Open search: fainter sources and unfamiliar spectral shapes",
    "",
    "Research owner: Biswajit Jana. Research snapshot: 13 November 2026.",
    "No new physical class, black-hole star, or cosmological anomaly is established.",
    "",
    "## What was executed",
    "",
    f"The added fixed sample contains **{summary['new_attempted']:,} sources**, all with "
    "PRISM grade 3, 2 ≤ archive median pixel S/N < 5, and 0.5 ≤ z < 12. "
    "The highest-S/N product per upstream source ID was selected, excluding the previous "
    "sample. No colour or physical-class filter was applied.",
    "",
    f"**{summary['new_analysed']:,} new sources analysed; {summary['new_failed']} failures.** "
    f"The combined experiment now covers **{summary['combined_analysed']:,} upstream source "
    f"records**, {100 * summary['fraction_of_parent']:.1f}% of the 20,144-record parent catalogue. "
    "This is the full stated S/N ≥ 2 reliable-redshift subset, not a complete JWST census. "
    "Upstream IDs are not independently deduplicated physical galaxies.",
    "",
    f"The unchanged earlier continuum-plus-blend cuts add **{summary['new_joint_priorities']} "
    f"joint priorities**; the combined count is {summary['combined_joint_priorities']}. "
    "A non-selection in this fainter sample is not proof of physical absence: continuum "
    "coverage and measurement uncertainty limit sensitivity.",
    "",
    "| Assigned survey | Combined analysed | Added sources |",
    "| --- | ---: | ---: |",
]
for survey, group in measurements.groupby("survey"):
    lines.append(f"| {survey} | {len(group):,} | {int((new.survey == survey).sum()):,} |")
lines += [
    "",
    "![Coverage and rarity](../figures/open_search_overview.png)",
    "",
    "## Search without an assumed object class",
    "",
    "Each spectrum is integrated into 12 logarithmically spaced rest-frame bands in "
    "0.32–0.68 µm (optical) and separately 0.14–0.30 µm (ultraviolet). Signed fν is "
    "preserved. Pixel-overlap weights propagate both variance and the covariance "
    "from pixels shared between adjacent bands. Every band needs ≥80% coverage and "
    "at least two contributing pixels; median signed band S/N must be ≥3. Missing "
    "bands are never imputed. These bands include emission lines and are not pure "
    "continuum measurements.",
    "",
    "Shapes are divided by RMS flux, retaining slope and signed structure. The "
    "target is excluded from its peers. Peers must be from the same assigned survey, "
    "within |Δz|/(1+z_target)<0.15 and a factor of two in archive S/N; at least 30 "
    "are required. The score is the RMS residual relative to the peer median, "
    "scaled by peer MAD (minimum 0.05) and approximate target band errors. It is "
    "not a likelihood or a physical-template fit. Shared normalisation creates "
    "bin dependence; the score is not chi-square.",
    "",
    "We perturb measured band fluxes 128 times using their covariance, normalise "
    "each draw, and rank by the 10th percentile of the descriptive score. This "
    "is a conditional noise-sensitivity check, **not** a posterior confidence bound "
    "or a false-discovery calibration. Rank percentiles are calculated separately "
    "for the two wavelength windows, sorted together, deduplicated by upstream ID, "
    "and the top 12 selected for inspection. A fixed queue length does not assert "
    "that 12 scientifically significant anomalies exist.",
    "",
    f"**{summary['unique_shape_scored_sources']:,} unique sources were shape-scored**, "
    f"with {summary['scored_source_mode_pairs']:,} usable source/window pairs. "
    f"Of the 12 queued sources, {summary['shape_queue_overlap_joint_screen']} also pass "
    "the earlier joint screen. The two methods ask different questions.",
    "",
    "| Window | Status | Source records |",
    "| --- | --- | ---: |",
]
for r in summary["shape_status_counts"]:
    lines.append(f"| {r['mode']} | {r['status']} | {r['count']:,} |")
lines += [
    "",
    "## Measured follow-up queue",
    "",
    "Scores express relative rarity within restricted empirical peer groups. "
    "They do not demonstrate novelty in the literature or failure of established physics.",
    "",
    "| Source | Survey ID | z | Window | Score | After removing strongest band | Other-survey score |",
    "| --- | --- | ---: | --- | ---: | ---: | ---: |",
]
for r in inspection:
    cross = "unavailable" if r["other_survey_score"] is None else f"{r['other_survey_score']:.2f}"
    lines.append(
        f"| {r['source_id']} | {r['survey']} {r['source_number']} | {r['z_best']:.3f} | "
        f"{r['mode']} | {r['score']:.2f} | {r['score_without_dominant_bin']:.2f} | {cross} |"
    )
lines += [
    "",
    "![Observed shape comparisons](../figures/open_shape_comparison.png)",
    "",
    "Grey regions show empirical peer scatter, not uncertainty on a physical model. "
    "Target bars show released statistical error divided by RMS; the full "
    "normalisation covariance is handled only by the perturbation test.",
    "",
    "## Attempts to break the anomalies",
    "",
    "For every queued source, we remove the largest-residual band, renormalise "
    "both target and peers using the remaining bands, and repeat the distance "
    "calculation. We also shift (1+z) by ±1% and reintegrate the "
    "spectrum against fixed peers. This is a deliberately broad stress test, "
    "not an inferred redshift error. Other-survey comparison uses the same "
    "redshift/S/N tolerances and at least 30 peers. These checks were added "
    "after fixing the search method and do not retroactively change the queue.",
    "",
    "| Source | Score at −1% / +1% in (1+z) | Alternate PRISM check |",
    "| --- | --- | --- |",
]
for r in inspection:
    shift = " / ".join(
        f"{s['score']:.2f}" if "score" in s else s["status"]
        for s in r["redshift_stress_fixed_peers"]
    )
    alt = r["alternate_product_check"]
    alt_text = alt["status"] + (f"; score {alt['score']:.2f}" if "score" in alt else "")
    lines.append(f"| {r['source_id']} | {shift} | {alt_text} |")
lines += [
    "",
    "Alternate products share the upstream source association and reduction; "
    "they are not demonstrated independent exposures or independent extractions. "
    "A missing counterpart is missing evidence, not confirmation.",
    "",
    "The UV shape of JADES-GS-13905 changes substantially between its two products "
    "(score 2.66 versus 0.73); the lower-S/N alternative does not support the same "
    "degree of rarity. This is a measurement/depth concern, not established "
    "variability. JADES-GN-38147 is more repeatable between products (3.85 versus "
    "3.71), but its score drops to 1.54 after the strongest band is removed and "
    "the remaining shape renormalized. Neither check proves a new continuum family.",
    "",
    "## Recovering a known LRD and testing resolution",
    "",
    "The queue recovers CEERS-3506, an already studied extreme emission-line "
    "source ([Zhu et al.](https://arxiv.org/abs/2410.12198v2)), and the known LRD "
    "JADES-GS-159717 ([D'Eugenio et al.](https://arxiv.org/abs/2506.14870v2)). "
    "The latter passed our old broad-blend cut but failed the continuum-ratio "
    "requirement: its measured ratio is 1.33 rather than >2. The class-agnostic "
    "search therefore recovers an interesting known object missed by that conjunction. "
    "JADES-GS-13905 is also present at matching coordinates in "
    "[Chen et al.](https://arxiv.org/html/2505.24080v2). The full queue has not "
    "undergone an exhaustive novelty audit.",
    "",
    "After the literature lookup, we downloaded the matched DJA G395M and G395H "
    "spectra of JADES-GS-159717. This is explicitly a targeted benchmark check. "
    "We fixed shoulder windows [-450,-250] and [250,450] km/s and a central "
    "window [-120,120] km/s, and measured (mean shoulders minus core) at the "
    "archive redshift and published systemic z=5.07781. Pixel-overlap covariance "
    "is propagated. A positive value would describe a local central dip, not "
    "identify its physical cause.",
    "",
    "| Grating | Redshift basis | Shoulder minus core (µJy) | Statistical error (µJy) |",
    "| --- | --- | ---: | ---: |",
]
for spectrum in resolution_check["spectra"]:
    for check in spectrum["diagnostics"]:
        lines.append(
            f"| {spectrum['grating']} | {check['redshift_basis']} | "
            f"{check['contrast_uJy']:.3f} | {check['statistical_error_uJy']:.3f} |"
        )
lines += [
    "",
    "![Observed grating comparison](../figures/recovered_lrd_resolution.png)",
    "",
    "**The simple central-dip test fails here:** the central narrow emission "
    "peak dominates, giving negative contrasts in both gratings. The high-resolution "
    "profile resolves additional troughs flanking that peak. These measurements "
    "neither reproduce nor refute the published emission-plus-absorption "
    "decomposition. They show why a fixed core statistic cannot replace a full "
    "multi-component fit with the instrumental LSF and extraction checks. "
    "The paper's physical conclusions must not be attributed to this simpler test.",
    "",
    "## Synthetic method check",
    "",
    "A separate, explicitly synthetic test uses 300 noisy log-quadratic "
    "continuum peers, 500 calibration trials, and 500 independent trials per "
    "scenario. A 99th-percentile synthetic-null threshold is used only for this "
    "test; it is not applied as a discovery threshold to the real catalogue.",
    "",
    "| Synthetic scenario | Selected / trials |",
    "| --- | ---: |",
]
for r in validation["outcomes"]:
    lines.append(f"| {r['scenario']} | {r['selected']} / {r['trials']} |")
lines += [
    "",
    "The deliberately strong artificial defects also trigger the detector. "
    "Thus score alone cannot distinguish unusual astrophysics from defects. "
    "This simplified experiment omits extraction/calibration systematics and "
    "real population diversity; it does not estimate the real-source false-positive rate. "
    "Exact binomial intervals and assumptions are in the "
    "[validation record](../results/open_search/synthetic_validation.json).",
    "",
    "## Actual source images and spectra",
    "",
    "![Observed image atlas](../figures/open_rgb_atlas.png)",
    "",
    "The cutouts contain observed NIRCam pixels, with display stretches and "
    "unmatched PSFs stated on the figures. Image brightness cannot be compared "
    "directly between independently stretched panels. Nearby components are not "
    "automatically merged or assigned to the spectrum.",
    "",
]
for r in inspection:
    lines += [
        f"### {r['survey']} {r['source_number']} — {r['source_id']}",
        "",
        notes.get(r["source_id"], "Physical explanation and novelty remain unresolved."),
        "",
    ]
    if "evidence_png" in r:
        lines += [f"![Observed evidence for {r['source_id']}](../{r['evidence_png']})", ""]
lines += [
    "## Evolution and cosmology",
    "",
    "The [evolution and unknown-object research note](evolution_and_unknown_objects.md) "
    "compares competing published interpretations, defines tests for unfamiliar "
    "phenomena, and explains the conditional growth calculation. We computed 180 "
    "seed/efficiency/duty/redshift scenarios. No mass or future fate of a source "
    "was inferred from its anomaly score.",
    "",
    "![Conditional cosmic growth scenarios](../figures/cosmic_growth_context.png)",
    "",
    "## Reproduction and remaining limits",
    "",
    "The [configuration](../configs/open_search.yaml), "
    "[run summary](../results/open_search/summary.json), "
    "[all signed measurements](../results/open_search/combined_measurements.csv), "
    "[shape scores](../results/open_search/shape_scores.csv), "
    "[full inspection audit](../results/open_search/inspection_audit.json) and "
    "[download ledger](../results/open_downloads.json) preserve denominators and failures. "
    "Compressed JSON Lines files retain band covariance, source hashes and peer IDs. "
    "Raw files remain outside Git under a 4 GiB cap.",
    "",
    "After restoring the earlier experiment and pinned raw inputs:",
    "",
    "```sh",
    "python scripts/replay_manifest.py --budget-mib 4096",
    "python scripts/download_open_sample.py",
    "python scripts/analyse_open_sample.py",
    "python scripts/validate_open_search.py",
    "python scripts/inspect_open_sources.py",
    "python scripts/check_recovered_lrd_resolution.py",
    "python scripts/cosmic_growth_context.py",
    "python scripts/build_open_report.py",
    "```",
    "",
    "The method remains restricted to reliable archive redshifts and sufficiently "
    "covered, positive median band S/N. It can miss narrow-line anomalies, faint "
    "sources, redshift mistakes and rare sources without enough peers. Same-survey "
    "matching does not eliminate within-survey depth or calibration differences. "
    "Selection effects prevent inferring abundance or cosmic evolution from the "
    "histogram. A complete physical-model comparison and independent reductions "
    "remain necessary before any new-object claim.",
    "",
]
(root / "reports/open_search.md").write_text("\n".join(lines), encoding="utf-8")
print("Wrote reports/open_search.md and figures/open_search_overview")
