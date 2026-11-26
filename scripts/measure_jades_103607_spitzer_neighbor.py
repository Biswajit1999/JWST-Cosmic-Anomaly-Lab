"""Test the February 2006 IRAC residual against neighbouring image structure."""

import json

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from download_jwst_cross_detector import OUT, ROOT
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from measure_jades_103607_spitzer_epochs import GRID, SCALE, SOURCE, inspect_year, released_mosaics
from scipy.ndimage import gaussian_filter
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

STEM = "module_b_001_spitzer_neighbor"
HISTORY = OUT / "module_b_001_spitzer_epochs.json"
WIDTHS = (1.55, 1.93, 2.35)
INNER = (4.2, 7.2)
N = len(GRID) // 2
YY, XX = np.meshgrid(GRID * SCALE, GRID * SCALE, indexing="ij")
R = np.hypot(XX, YY)
INJECTIONS_NJY = (25, 250, 500, 1000)


def spatial_basis(reference):
    valid = np.isfinite(reference)
    blurred = gaussian_filter(np.where(valid, reference, 0), 0.65)
    denominator = gaussian_filter(valid.astype(float), 0.65)
    blurred = np.divide(
        blurred,
        denominator,
        out=np.zeros_like(blurred),
        where=denominator > 0.9,
    )
    dy, dx = np.gradient(blurred)
    dyy, dyx = np.gradient(dy)
    _, dxx = np.gradient(dx)
    return np.stack(
        (
            np.ones_like(reference),
            XX,
            YY,
            blurred,
            dx,
            dy,
            dxx,
            dyy,
            dyx,
        ),
        axis=-1,
    ), denominator > 0.995


def match_background(reference, epoch, inner):
    """Predict the epoch image from a reference with gain, shift and blur terms."""
    basis, supported = spatial_basis(reference)
    finite = np.isfinite(epoch) & supported
    zone = finite & (R >= inner) & (R <= 27)
    feature, target = basis[zone], epoch[zone]
    scales = np.sqrt(np.mean(feature**2, axis=0))
    if len(target) < 500 or not np.all(scales > 0):
        raise ValueError("Insufficient common coverage for IRAC nuisance fit")
    normalized = feature / scales
    used = np.ones(len(target), dtype=bool)
    for _ in range(4):
        fitted = np.linalg.lstsq(normalized[used], target[used], rcond=None)[0]
        residual = target - normalized @ fitted
        center = float(np.median(residual[used]))
        scatter = float(1.4826 * np.median(abs(residual[used] - center)))
        used = abs(residual - center) <= 4.5 * scatter
    coefficients = fitted / scales
    model = basis @ coefficients
    return np.where(finite, epoch - model, np.nan), dict(
        training_pixels=len(target),
        clipped_background_scatter_njy_per_pixel=scatter,
        retained_pixels=int(used.sum()),
        coefficients=coefficients.tolist(),
    )


def psf_template(x, y, width):
    sigma = width / 2.354820045
    s2 = (XX - x) ** 2 + (YY - y) ** 2
    return np.exp(-0.5 * s2 / sigma**2) * SCALE**2 / (2 * np.pi * sigma**2)


def fit_source(data, x, y, width):
    dx, dy = XX - x, YY - y
    keep = (np.hypot(dx, dy) <= 2.7) & np.isfinite(data)
    if keep.sum() < 52:
        return np.nan
    template = psf_template(x, y, width)
    features = np.column_stack((template[keep], np.ones(keep.sum()), dx[keep], dy[keep]))
    return float(np.linalg.lstsq(features, data[keep], rcond=None)[0][0])


def evaluate(reference, epoch, inner):
    delta = epoch - reference
    corrected, diagnostics = match_background(reference, epoch, inner)
    rows = []
    for width in WIDTHS:
        controls = np.array(
            [
                fit_source(corrected, i * SCALE, j * SCALE, width)
                for i in range(-40, 41, 8)
                for j in range(-40, 41, 8)
                if i * i + j * j >= 17**2
            ]
        )
        controls = controls[np.isfinite(controls)]
        if len(controls) < 35:
            raise ValueError("Insufficient colocated control apertures")
        center = float(np.median(controls))
        noise = float(1.4826 * np.median(abs(controls - center)))
        flux = fit_source(corrected, 0, 0, width) - center
        row = dict(
            fwhm_arcsec=width,
            point_amplitude_njy=flux,
            control_median_njy=center,
            control_scatter_njy=noise,
            control_count=len(controls),
            control_fluxes_njy=controls.tolist(),
        )
        if width == WIDTHS[1]:
            injected_results = []
            source = psf_template(0, 0, width)
            for amplitude in INJECTIONS_NJY:
                trial, _ = match_background(reference, epoch + amplitude * source, inner)
                recovery = fit_source(trial, 0, 0, width) - center - flux
                detection = np.mean(controls - center + recovery > 3 * noise)
                injected_results.append(
                    dict(
                        input_same_band_flux_njy=amplitude,
                        recovered_increment_njy=recovery,
                        control_fraction_above_three_scatter=float(detection),
                    )
                )
            row["synthetic_injection"] = injected_results
        rows.append(row)
    return (
        dict(protected_radius_arcsec=inner, fit=diagnostics, template_checks=rows),
        delta,
        corrected,
    )


def figure(records, panels):
    fig = Figure(figsize=(9.4, 6.3), layout="constrained")
    axes = fig.subplots(len(records), 3, squeeze=False)
    k = 15
    cut = np.s_[N - k : N + k + 1, N - k : N + k + 1]
    extent = np.array([-1, 1, -1, 1]) * (k + 0.5) * SCALE
    for axset, (rec, (epoch, raw, corrected)) in zip(
        axes, zip(records, panels, strict=True), strict=True
    ):
        intensity = float(np.nanpercentile(abs(epoch[cut]), 99))
        limit = float(np.nanpercentile(abs(raw[cut]), 97))
        labels = (
            "February 2006",
            f"2006 minus {rec['reference_year']}",
            "After local PSF matching",
        )
        for j, (ax, data, title) in enumerate(
            zip(axset, (epoch, raw, corrected), labels, strict=True)
        ):
            ax.imshow(
                data[cut],
                origin="lower",
                extent=extent,
                cmap="magma" if j == 0 else "RdBu_r",
                vmin=-0.1 * intensity if j == 0 else -limit,
                vmax=intensity if j == 0 else limit,
            )
            ax.add_patch(
                Circle((0, 0), 1.8, fill=False, color="#62d3d3" if j == 0 else "#244957", lw=1)
            )
            ax.set(title=title, xticks=[], yticks=[])
        m = rec["variants"][-1]["template_checks"][1]
        axset[-1].set_xlabel(
            f"{m['point_amplitude_njy']:+.0f} ± {m['control_scatter_njy']:.0f} nJy"
        )
    fig.suptitle("JADES 103607 · Spitzer 4.5 µm residual and neighbouring light", fontsize=13)
    save_figure(fig, FIGURES / STEM)


def main():
    position = json.loads(SOURCE.read_text())["icrs_deg"]
    sky = SkyCoord(*position, unit="deg")
    grid = sky.spherical_offsets_by(-XX * u.arcsec, YY * u.arcsec)
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    historical, receipts = released_mosaics(cache, grid)
    warm_sources, warm = inspect_year(cache, sky, grid, "2011")
    epoch = historical["2"]["2006"][0]
    references = {"2004": historical["2"]["2004"][0], "2011": warm["2"][0]}
    records, panels = [], []
    for year, ref in references.items():
        variants = []
        for inner in INNER:
            row, raw, residual = evaluate(ref, epoch, inner)
            variants.append(row)
            measure = row["template_checks"][1]
            print(
                year,
                inner,
                measure["point_amplitude_njy"],
                measure["control_scatter_njy"],
                flush=True,
            )
        records.append(dict(reference_year=year, variants=variants))
        panels.append((epoch, raw, residual))
    figure(records, panels)
    write_json(
        OUT / f"{STEM}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=position,
            source_position_sha256=text_digest(SOURCE),
            prior_spitzer_comparison_sha256=text_digest(HISTORY),
            cold_image_ranges=receipts["2"],
            warm_2011_products=warm_sources,
            template_widths_arcsec=list(WIDTHS),
            comparisons=records,
            scope="First-order spatial matching with target excluded. Fits remain conditional on local PSF and "
            "background models. Injections are synthetic unresolved sources in the IRAC channel and "
            "are not a conversion of the October 2026 JWST photometry.",
        ),
    )


if __name__ == "__main__":
    main()
