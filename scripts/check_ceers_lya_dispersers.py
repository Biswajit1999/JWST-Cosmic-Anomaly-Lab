"""Apply the Lyα screen to optical-line-anchored CEERS G140M and PRISM spectra."""

from itertools import product

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
from astropy.io import fits
from check_jades_reionization_prism import (
    BREAK_OFFSETS_A,
    DISP,
    NOMINAL,
    RESOLUTION_FACTORS,
    SLOPES,
    local_fit,
)
from compare_jades_reionization_dja import INVENTORY
from download_jades_reionization_sample import ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from measure_jades_reionization import C_KMS, LYA, WIDTHS, best_model, line_models, spectral_arrays

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

OUT = ROOT / "results/jades_reionization_full"
STEM = "jwst_ceers_lya_crosscheck"
LITERATURE = {
    1027: {
        "result": "published Lyα emitter",
        "reference": "Tang et al. 2023",
        "url": "https://doi.org/10.1093/mnras/stad2763",
    },
    1023: {
        "result": "published Lyα non-detection",
        "reference": "Napolitano et al. 2024",
        "url": "https://doi.org/10.1051/0004-6361/202449644",
    },
}


def matched_sources():
    table = pd.read_csv(INVENTORY, low_memory=False)
    table = table[
        (table.survey == "CEERS")
        & table.grade.eq(3)
        & table.z_best.between(5.5, 9, inclusive="left")
    ].copy()
    table["lya_um"] = LYA * (1 + table.z_best) / 10000
    cov = (table.wmin < table.lya_um - 0.01) & (table.wmax > table.lya_um + 0.03)
    g140m = table[
        cov
        & table.grating.eq("G140M")
        & ((table["filter"] == "F070LP") | ((table["filter"] == "F100LP") & table.lya_um.ge(0.98)))
    ]
    prism = table[cov & table.grating.eq("PRISM") & table["filter"].eq("CLEAR")]
    oxy = table[
        (table.line_oiii_5007_err > 0)
        & (table.line_oiii_5007 / table.line_oiii_5007_err >= 5)
        & (table.wmax > 0.500824 * (1 + table.z_best))
    ].copy()
    oxy["oxygen_snr"] = oxy.line_oiii_5007 / oxy.line_oiii_5007_err
    for sid in sorted(set(g140m.source_id) & set(prism.source_id) & set(oxy.source_id)):
        g = (
            g140m[g140m.source_id == sid]
            .sort_values(["exptime", "file"], ascending=[False, True])
            .iloc[0]
        )
        p = (
            prism[prism.source_id == sid]
            .sort_values(["exptime", "file"], ascending=[False, True])
            .iloc[0]
        )
        o = (
            oxy[oxy.source_id == sid]
            .sort_values(["oxygen_snr", "file"], ascending=[False, True])
            .iloc[0]
        )
        sep = SkyCoord(g.ra, g.dec, unit="deg").separation(SkyCoord(p.ra, p.dec, unit="deg")).arcsec
        sep_o = (
            SkyCoord(g.ra, g.dec, unit="deg").separation(SkyCoord(o.ra, o.dec, unit="deg")).arcsec
        )
        redshift_spread = max(abs(g.z_best - p.z_best), abs(g.z_best - o.z_best)) / (1 + o.z_best)
        selection = dict(
            coordinate_spread_arcsec=float(max(sep, sep_o)),
            relative_redshift_spread=float(redshift_spread),
            inventory_products=[g.file, p.file, o.file],
            passed=bool(max(sep, sep_o) <= 0.15 and redshift_spread <= 0.002),
        )
        yield sid, g, p, o, selection


def download(cache, row, position):
    path = cache.fetch(
        row.url,
        f"spectra/{row.file}",
        survey="DAWN JWST Archive",
        version="v4",
        kind="observed_spectrum",
    )
    spec = read_spectrum(path, wavelength_range=(0.65, 5.4))
    sep = position.separation(
        SkyCoord(spec.metadata["SRCRA"], spec.metadata["SRCDEC"], unit="deg")
    ).arcsec
    if (
        sep > 0.15
        or str(spec.metadata["GRATING"]).upper() != row.grating
        or str(spec.metadata["FILTER"]).upper() != row["filter"]
    ):
        raise ValueError(f"Spectrum header inconsistent with CEERS archive index: {row.file}")
    return spectral_arrays(spec), dict(
        product=receipt(path, row.url),
        root=row.root,
        archive_exposure_identifier=row.dataset,
        header_separation_arcsec=float(sep),
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cache = DownloadCache(ROOT, 28_000 * 1024**2)
    path = cache.fetch(
        DISP,
        "metadata/jwst_nirspec_prism_disp.fits",
        survey="STScI",
        version="linked JDox curve",
        kind="instrument_reference",
    )
    with fits.open(path) as hdus:
        wl, r = (np.array(hdus[1].data[name], float) for name in ("WAVELENGTH", "R"))
    results, displays, excluded = [], [], []
    for sid, g, p, o, selection in sorted(matched_sources(), key=lambda entry: entry[1].srcid):
        if not selection["passed"]:
            excluded.append(dict(source_id=sid, ceers_id=int(g.srcid), selection=selection))
            continue
        position = SkyCoord(g.ra, g.dec, unit="deg")
        gdata, gs = download(cache, g, position)
        pdata, ps = download(cache, p, position)
        z = float(o.z_best)
        zero = LYA * (1 + z)
        models = line_models(gdata, zero, WIDTHS[0])
        line = best_model(models)
        offset = line["offset_kms"] if line is not None else 350
        resolving_power = float(np.interp(zero / 1e4, wl, r))
        fit, ppm = local_fit(pdata, z, offset, resolving_power, *NOMINAL)
        prism_variants = []
        for scenario in product(RESOLUTION_FACTORS, SLOPES, BREAK_OFFSETS_A):
            variant, _ = local_fit(pdata, z, offset, resolving_power, *scenario)
            if variant is not None:
                prism_variants.append(variant)
        matched_flux = (
            dict(
                offset_kms=offset,
                flux=line["flux"],
                error=line["error"],
                snr=line["flux"] / line["error"],
                coverage=line["coverage"],
            )
            if line
            else None
        )
        item = dict(
            source_id=sid,
            ceers_id=int(g.srcid),
            ra=float(g.ra),
            dec=float(g.dec),
            z=z,
            cosmic_age_gyr=float(Planck18.age(z).value),
            oxygen_anchor=dict(
                file=o.file, url=o.url, archive_oxygen_snr=float(o.oxygen_snr), grade=int(o.grade)
            ),
            grating=gs | {"measurement": matched_flux},
            prism=ps | {"measurement": fit, "sensitivity_fits": prism_variants},
            searched_offset_kms=offset,
            grating_offset_count=len(models),
            selection=selection,
            literature_context=LITERATURE.get(int(g.srcid)),
        )
        results.append(item)
        displays.append((item, gdata, ppm))
        print(
            item["ceers_id"],
            item["z"],
            None if matched_flux is None else matched_flux["snr"],
            None if fit is None else fit["snr"],
            flush=True,
        )
    fig = Figure(figsize=(9.6, 2.15 * len(displays) + 0.6), layout="constrained")
    axes = fig.subplots(len(displays), 2, squeeze=False)
    for index, ((row, data, prism), (left, right)) in enumerate(zip(displays, axes, strict=True)):
        w, _, f, e, good = data
        zero = LYA * (1 + row["z"])
        v = (w / zero - 1) * C_KMS
        near = (v > -1600) & (v < 3800)
        ff = np.where(good, f, np.nan)
        left.step(v[near], ff[near], where="mid", color="#267184")
        left.fill_between(
            v[near], (ff - e)[near], (ff + e)[near], step="mid", color="#267184", alpha=0.18
        )
        left.axvline(row["searched_offset_kms"], color=".4", lw=0.7, ls=":")
        left.set(
            xlim=(-1600, 3800),
            ylabel="Fλ (10⁻²⁰ cgs/Å)",
            title=f"CEERS {row['ceers_id']} · z={row['z']:.3f} · G140M",
        )
        if prism:
            wave, f, e, ok, pred, null = prism
            rest = wave / (1 + row["z"])
            right.errorbar(
                rest[ok],
                f[ok],
                yerr=e[ok],
                fmt="o",
                ms=3,
                capsize=2,
                color="#267184",
                label="Observed",
            )
            right.step(rest, pred, where="mid", color="#aa513c", label="Break + line")
            right.step(rest, null, where="mid", color="#775da2", ls="--", label="Break only")
            if index == 0:
                right.legend(frameon=False, loc="best", fontsize=7)
        else:
            right.text(
                0.5, 0.5, "Insufficient PRISM coverage", transform=right.transAxes, ha="center"
            )
        right.set(title="Observed PRISM")
    axes[-1, 0].set_xlabel("Velocity relative to systemic Lyα (km/s)")
    axes[-1, 1].set_xlabel("Rest wavelength (Å)")
    fig.suptitle("CEERS · NIRSpec spectra around Lyα", fontsize=13)
    save_figure(fig, ROOT / f"figures/{STEM}")
    write_json(
        OUT / "ceers_crosscheck.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            inventory_sha256=digest(INVENTORY),
            instrument_curve=receipt(path, DISP),
            targets=results,
            inventory_matches_excluded=excluded,
            selection="Grade 3, optical oxygen S/N >=5, distinct disperser spectra, bounded wavelength and coordinate overlap.",
            scope="Released reductions. Formal diagonal-noise significance; neither source priority nor catalogue completeness inferred.",
        ),
    )


if __name__ == "__main__":
    main()
