"""Bounded pixel-level coverage audit of all listed high-resolution priorities."""

from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.figure import Figure

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.spectroscopy import pixel_edges, read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/high_resolution_coverage"
out.mkdir(parents=True, exist_ok=True)
products = pd.read_csv(root / "results/line_anchor/priority_archive_products.csv")
products = products[products.grating.isin(["G140H", "G235H", "G395H"])]
selection = pd.concat(
    [
        pd.read_csv(root / "results/expanded/followup_priority.csv"),
        pd.read_csv(root / "results/open_search/shape_followup.csv"),
    ]
).set_index("source_id")
cache = DownloadCache(root, budget_bytes=4096 * 1024**2)
records, coverage = [], []
for product in products.itertuples():
    target = selection.loc[product.source_id]
    record = dict(
        source_id=product.source_id,
        file=product.file,
        grating=product.grating,
        url=product.url,
        z_fixed=float(target.z_best),
    )
    try:
        path = cache.fetch(
            product.url,
            "spectra/" + product.file,
            survey=target.survey,
            version="DJA-v4.4",
            kind="observed_spectrum",
        )
        spec = read_spectrum(path)
        prism = read_spectrum(root / "data/raw/spectra" / target.file)

        def coords(s):
            return SkyCoord(s.metadata["SRCRA"] * u.deg, s.metadata["SRCDEC"] * u.deg)

        separation = coords(spec).separation(coords(prism)).arcsec
        if separation > 0.05 or spec.metadata["GRATING"] != product.grating:
            raise ValueError("Header association check failed")
        record.update(
            status="INSPECTED",
            path=path.relative_to(root).as_posix(),
            sha256=digest(path),
            header_separation_arcsec=float(separation),
        )
        fig = Figure(figsize=(12, 5), layout="constrained")
        axes = fig.subplots(1, 2)
        for ax, (name, rest, limit) in zip(
            axes, [("Hα", 0.656461, 4500), ("[OIII] 5008", 0.500824, 900)], strict=True
        ):
            center = rest * (1 + target.z_best)
            velocity = (spec.wavelength.value / center - 1) * 299792.458
            edges = pixel_edges(velocity)
            overlap = np.maximum(0, np.minimum(edges[1:], limit) - np.maximum(edges[:-1], -limit))
            good = spec.good & (overlap > 0)
            fraction = float(overlap[spec.good].sum() / (2 * limit))
            status = (
                "COVERAGE_QUALIFIED"
                if fraction >= 0.8 and good.sum() >= 12
                else "INSUFFICIENT_COVERAGE"
            )
            row = dict(
                source_id=product.source_id,
                file=product.file,
                grating=product.grating,
                line=name,
                fixed_redshift=float(target.z_best),
                center_observed_um=float(center),
                half_window_kms=limit,
                good_coverage_fraction=fraction,
                valid_pixels=int(good.sum()),
                status=status,
                median_pixel_flux_over_error=float(
                    np.median(spec.flux.value[good] / spec.error.value[good])
                )
                if good.any()
                else None,
            )
            coverage.append(row)
            if good.any():
                ax.errorbar(
                    velocity[good],
                    spec.flux.value[good],
                    yerr=spec.error.value[good],
                    fmt=".",
                    ms=3,
                    lw=0.6,
                    color="#286d85",
                )
                lo = np.min((spec.flux.value - spec.error.value)[good])
                hi = np.max((spec.flux.value + spec.error.value)[good])
                margin = max((hi - lo) * 0.08, 1e-5)
                ax.set_ylim(lo - margin, hi + margin)
            else:
                ax.text(
                    0.5, 0.5, "No usable pixels", ha="center", va="center", transform=ax.transAxes
                )
            ax.set(
                xlim=(-limit, limit),
                xlabel="Velocity at fixed catalogue z (km/s)",
                ylabel="Released flux density (µJy)",
                title=f"{name} · {fraction:.0%} good coverage",
            )
            ax.axhline(0, color="gray", lw=0.6)
            ax.grid(alpha=0.15)
        fig.suptitle(
            f"OBSERVED {product.source_id} · z={target.z_best:.4f}\n{product.file}\nCoverage inspection only; no absorption classification",
            fontsize=10,
        )
        save_figure(
            fig, root / "figures/high_resolution_sources" / product.file.removesuffix(".spec.fits")
        )
    except (OSError, ValueError, KeyError, TypeError) as exc:
        record.update(status="UNAVAILABLE", reason=str(exc))
    records.append(record)
    write_json(out / "products.json", records)
    print(product.file, record["status"], flush=True)
table = pd.DataFrame(coverage)
table.to_csv(out / "line_coverage.csv", index=False)
write_json(
    out / "summary.json",
    dict(
        executed_utc=utcnow(),
        research_snapshot="2026-11-15",
        source_records=int(products.source_id.nunique()),
        listed_products=len(products),
        inspected_products=sum(r["status"] == "INSPECTED" for r in records),
        line_windows=len(coverage),
        coverage_qualified_windows=int((table.status == "COVERAGE_QUALIFIED").sum()),
        interpretation="Pixel-level coverage only; median pixel flux/error includes continuum, not line significance",
        physical_classification="NOT_ESTABLISHED",
        script_sha256=text_digest(Path(__file__)),
        inventory_sha256=text_digest(root / "results/line_anchor/priority_archive_products.csv"),
    ),
)
print("Completed", len(records), "products and", len(coverage), "line windows", flush=True)
