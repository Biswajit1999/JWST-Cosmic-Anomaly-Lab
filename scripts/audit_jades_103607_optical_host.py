"""Check the public 3D-HST catalogue and MUSE-Wide spectrum near JADES 103607."""

import json
from io import StringIO
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from scipy.special import ndtr
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.spectroscopy import pixel_edges
from cosmic_anomaly_lab.visual_evidence import save_figure

STEM = "module_b_001_optical_host"
VIZIER = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"
MUSE = "https://musewide.aip.de"
CATALOGUES = {
    "three_d_hst": "J/ApJS/225/27/3dhst_sp",
    "muse_photometric": "J/A+A/624/A141/photcat",
    "muse_emission": "J/A+A/624/A141/emobj",
}
# Vacuum rest wavelengths in angstroms, with equal amplitudes for listed doublet components.
LINES = {
    "C IV": (1548.20, 1550.77),
    "He II": (1640.42,),
    "O III]": (1660.81, 1666.15),
    "C III]": (1906.68, 1908.73),
    "C II]": (2326.0,),
    "Mg II": (2796.35, 2803.53),
}
SIGMAS = (1.25, 2.1)
FIT_RADIUS_A = 15
COVERAGE = 0.98
REDSHIFT_GRID = np.round(np.arange(1.85, 2.601, 0.001), 3)


def catalogue(cache, target, tag, table):
    params = {
        "-source": table,
        "-c": f"{target.ra.deg:.8f} {target.dec.deg:.8f}",
        "-c.rs": 8,
        "-out.max": 100,
    }
    url = VIZIER + "?" + urlencode(params)
    path = cache.fetch(
        url,
        f"metadata/jades_103607_{tag}.tsv",
        survey="CDS VizieR",
        version=table,
        kind="published_counterpart_catalogue",
    )
    content = path.read_text(encoding="utf-8")
    if f"#Name: {table}" not in content or "#Truncated" in content:
        raise ValueError(f"Incomplete or unidentified VizieR response for {table}")
    try:
        data = pd.read_csv(StringIO(content), sep="\t", comment="#")
    except pd.errors.EmptyDataError:
        data = pd.DataFrame()
    if not data.empty:
        for c in ("RAJ2000", "DEJ2000"):
            data[c] = pd.to_numeric(data[c], errors="coerce")
        data = data.dropna(subset=["RAJ2000", "DEJ2000"]).copy()
        data["separation_arcsec"] = target.separation(
            SkyCoord(data.RAJ2000, data.DEJ2000, unit="deg")
        ).arcsec
        data = data.sort_values("separation_arcsec")
    return dict(
        product=receipt(path, url),
        radius_arcsec=8,
        parsed_rows=len(data),
        nearby_rows=json.loads(data.head(6).to_json(orient="records", double_precision=10)),
    )


def air(vacuum):
    wavenumber_squared = (1e4 / vacuum) ** 2
    n = (
        1
        + 0.05792105 / (238.0185 - wavenumber_squared)
        + 0.00167917 / (57.362 - wavenumber_squared)
    )
    return vacuum / n


def line_fit(wave, edges, flux, error, rest, z, sigma):
    components = air(np.array(rest) * (1 + z))
    centre = float(components.mean())
    lo = int(np.searchsorted(wave, components[0] - FIT_RADIUS_A))
    hi = int(np.searchsorted(wave, components[-1] + FIT_RADIUS_A, side="right"))
    if lo == 0 or hi == len(wave):
        return None
    local_wave, y, unc = wave[lo:hi], flux[lo:hi], error[lo:hi]
    bounds = edges[lo : hi + 1]
    good = np.isfinite(y) & np.isfinite(unc) & (unc > 0)
    if good.sum() < 15:
        return None
    area = np.diff(ndtr((bounds[:, None] - components) / sigma), axis=0).mean(axis=1)
    if area[good].sum() < COVERAGE:
        return None
    design = np.column_stack(
        (area / np.diff(bounds), np.ones(len(y)), (local_wave - centre) / FIT_RADIUS_A)
    )[good]
    weighted = design / unc[good, None]
    covariance = np.linalg.inv(weighted.T @ weighted)
    values = covariance @ weighted.T @ (y[good] / unc[good])
    e = float(np.sqrt(covariance[0, 0]))
    return dict(
        centre_air_a=centre,
        sigma_a=sigma,
        flux_native_integral=float(values[0]),
        error_native_integral=e,
        signed_snr=float(values[0] / e),
        used_bins=int(good.sum()),
    )


def plot(wave, flux, error, scans, fixed):
    fig = Figure(figsize=(10, 6.3), layout="constrained")
    top, bottom = fig.subplots(2, 1)
    valid = np.isfinite(flux) & np.isfinite(error) & (error > 0)
    top.plot(wave, np.where(valid, flux / error, np.nan), lw=0.55, color="#275f77")
    top.axhline(0, color=".55", lw=0.7)
    for ion, m in fixed.items():
        if m:
            top.axvline(m["centre_air_a"], color="#ad5b32", lw=0.5, alpha=0.6)
            top.annotate(
                ion,
                (m["centre_air_a"], 0.98),
                xycoords=("data", "axes fraction"),
                rotation=90,
                ha="right",
                va="top",
                fontsize=8,
                color="#8b482a",
            )
    top.set(
        ylabel="Released flux / error",
        xlabel="Observed air wavelength (Å)",
        title="MUSE-Wide spectrum at CANDELS 10658 · marked wavelengths use z = 2.47",
    )
    x = [row["z"] for row in scans]
    for key, c, style in (("max_snr", "#216782", "-"), ("second_snr", "#a05b37", "-")):
        bottom.plot(
            x,
            [row[key] for row in scans],
            color=c,
            ls=style,
            lw=0.85,
            label="Strongest complex" if key == "max_snr" else "Second strongest",
        )
    bottom.axhline(5, color=".45", ls=":", lw=0.7)
    bottom.axhline(3.5, color=".45", ls=":", lw=0.7)
    bottom.legend(frameon=False, loc="upper right")
    bottom.set(
        xlabel="Trial redshift",
        ylabel="Formal line S/N",
        title="Search for two emission complexes at one redshift",
    )
    save_figure(fig, FIGURES / STEM)


def main():
    source = OUT / "module_b_001_native_epoch_flux.json"
    target = SkyCoord(*json.loads(source.read_text())["icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, 30_000 * 1024**2)
    context = {tag: catalogue(cache, target, tag, table) for tag, table in CATALOGUES.items()}
    hst = context["three_d_hst"]["nearby_rows"][0]
    muse = context["muse_photometric"]["nearby_rows"][0]
    if (
        int(hst["ID"]) != 21067
        or int(muse["GuoId"]) != 10658
        or max(hst["separation_arcsec"], muse["separation_arcsec"]) > 0.3
    ):
        raise ValueError("Inconsistent catalogue counterpart positions")
    page = MUSE + "/dr1/photometric_catalog/10658/"
    product = (
        MUSE
        + "/files/dr1/photometric_catalog_spectra/candels-cdfs-26/tdose_spectrum_candels-cdfs-26_10658.fits"
    )
    metadata = cache.fetch(
        page,
        "metadata/musewide_10658.html",
        survey="MUSE-Wide",
        version="DR1",
        kind="public_spectrum_listing",
    )
    if product not in metadata.read_text(encoding="utf-8"):
        raise ValueError("Spectrum link is absent from the MUSE-Wide source page")
    path = cache.fetch(
        product,
        "spectra/musewide_10658.fits",
        survey="MUSE-Wide",
        version="DR1 TDOSE",
        kind="observed_optical_spectrum",
    )
    with fits.open(path) as hdus:
        wave, flux, error = (
            np.array(hdus["SPEC1D"].data[k], float) for k in ("wave", "flux", "fluxerror")
        )
        units = {key: hdus["SPEC1D"].columns[key].unit for key in ("wave", "flux", "fluxerror")}
    edges = pixel_edges(wave)
    rows = []
    for z in REDSHIFT_GRID:
        selected = {
            ion: fit
            for ion, lam in LINES.items()
            if (fit := line_fit(wave, edges, flux, error, lam, z, SIGMAS[0]))
        }
        rank = sorted(selected, key=lambda k: selected[k]["signed_snr"], reverse=True)
        rows.append(
            dict(
                z=float(z),
                usable_complexes=len(selected),
                max_snr=selected[rank[0]]["signed_snr"] if rank else None,
                max_ion=rank[0] if rank else None,
                second_snr=selected[rank[1]]["signed_snr"] if len(rank) > 1 else None,
                second_ion=rank[1] if len(rank) > 1 else None,
            )
        )
    hits = [
        item
        for item in rows
        if item["max_snr"] is not None
        and item["max_snr"] >= 5
        and item["second_snr"] is not None
        and item["second_snr"] >= 3.5
    ]
    fixed = {
        str(z): {
            ion: [line_fit(wave, edges, flux, error, lam, z, sigma) for sigma in SIGMAS]
            for ion, lam in LINES.items()
        }
        for z in (2.2459, 2.263, 2.47, 2.4858)
    }
    plot(wave, flux, error, rows, {ion: f[0] for ion, f in fixed["2.47"].items()})
    write_json(
        OUT / f"{STEM}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            target_input_sha256=text_digest(source),
            position_icrs_deg=[target.ra.deg, target.dec.deg],
            catalogues=context,
            spectrum=receipt(path, product),
            spectrum_listing=receipt(metadata, page),
            fits_units=units,
            wavelength_range_air_a=[float(wave.min()), float(wave.max())],
            rest_vacuum_lines_a=LINES,
            sigma_a=list(SIGMAS),
            local_continuum_window_padding_air_a=FIT_RADIUS_A,
            template_minimum_coverage=COVERAGE,
            fixed_photometric_redshift_fits=fixed,
            screen_hits=hits,
            screen_rows=rows,
            scope="Nearby host continuum extraction from 2015. Coordinates differ from the near-infrared transient. "
            "The spectrum has native flux units and formal diagonal errors. Line selection and bounds are limited "
            "to the listed complexes, widths and trial redshifts. This search is not a global false-alarm calibration.",
        ),
    )
    print(
        "fixed",
        [
            (z, [(ion, round(m[0]["signed_snr"], 2)) for ion, m in lines.items() if m[0]])
            for z, lines in fixed.items()
        ],
    )
    print("grid", len(rows), "hits", hits, "best", max(rows, key=lambda r: r["max_snr"] or -np.inf))


if __name__ == "__main__":
    main()
