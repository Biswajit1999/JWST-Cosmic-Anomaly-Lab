"""Search independent JWST visits for signed, position-matched flux changes."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib.figure import Figure
from scipy.spatial import cKDTree

from cosmic_anomaly_lab.epoch_search import aperture_score, fit_registration_residual, local_peaks
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/jwst_epoch_search"
config_path = root / "configs/jwst_epoch_search.yaml"
config = yaml.safe_load(config_path.read_text())
receipt = json.loads((out / "products.json").read_text())
assert receipt["config_sha256"] == text_digest(config_path)
images = receipt["images"]
assert {(r["band"], r["epoch"]) for r in images} == {
    (b, e) for b in config["filters"] for e in ["older", "newer"]
}
reference = images[0]["sci"]
with fits.open(root / reference["path"], memmap=False) as h:
    grid = WCS(h[0].header).celestial
scale_values = proj_plane_pixel_scales(grid)*3600
assert np.allclose(scale_values, .04, atol=1e-5)
scale = float(scale_values.mean())
cx, cy = grid.world_to_pixel_values(config["query_ra_deg"], config["query_dec_deg"])
rad = round(config["crop_half_width_arcsec"]/scale)
x0, y0 = int(round(float(cx)))-rad, int(round(float(cy)))-rad
side = 2*rad+1
corners_x, corners_y = np.array([0, side-1, 0, side-1]), np.array([0, 0, side-1, side-1])
sky = grid.pixel_to_world_values(corners_x+x0, corners_y+y0)


def load(record, kind):
    product = record[kind]
    path = root / product["path"]
    assert digest(path) == product["sha256"]
    with fits.open(path, memmap=False) as h:
        header = h[0].header
        assert str(header["FILTER"]).upper() == record["band"].split("-")[0]
        jwcs = WCS(header).celestial
        px, py = jwcs.world_to_pixel_values(*sky)
        sx, sy = int(round(px[0])), int(round(py[0]))
        residual = float(np.max(np.hypot(px-sx-corners_x, py-sy-corners_y)))
        if residual > .025:
            raise ValueError(f"{product['path']} is not on the common integer grid: {residual:.4f} px")
        assert 0 <= sx and 0 <= sy and sx+side <= header["NAXIS1"] and sy+side <= header["NAXIS2"]
        # A whole read decompresses the gzip member once instead of seeking row by row.
        result = np.array(h[0].data[sy:sy+side, sx:sx+side], dtype=np.float32)
        calibration = float(header["PHOTFNU"])*1e9
    assert 0 < calibration < 1e4
    return (result*calibration if kind == "sci" else result/calibration**2), residual


def stamp_measurements(band, y, x, sign):
    data = arrays[band]
    radius = config["aperture_radius_arcsec"]/scale
    half = int(np.ceil(.36/scale))
    sl = np.s_[y-half:y+half+1, x-half:x+half+1]
    yy, xx = np.mgrid[-half:half+1, -half:half+1]
    r = np.hypot(xx, yy)
    core = r <= radius
    annulus = (r >= .2/scale) & (r <= .36/scale)
    mask = data["valid"][sl]
    residual = data["difference"][sl]
    sign_sum = float(sign*np.sum(residual[core])) if mask[core].all() else None
    opposite = float(np.sum(np.maximum(-sign*residual[annulus & mask], 0)))
    annular_coverage = float(mask[annulus].mean())
    photometry = {}
    for epoch in ["older", "newer"]:
        image = data[epoch][sl]
        if mask[core].all() and annular_coverage >= .8:
            photometry[epoch+"_aperture_njy"] = float(image[core].sum()-core.sum()*np.median(image[mask & annulus]))
        else:
            photometry[epoch+"_aperture_njy"] = None
    return dict(score=float(data["score"][y, x]), delta_njy=float(data["flux"][y, x]),
                local_sigma_njy=float(data["noise"][y, x]), annular_coverage=annular_coverage,
                counter_flux_ratio=(opposite / sign_sum if sign_sum and sign_sum > 0 else None),
                **photometry)


arrays, calibrations, candidates = {}, [], []
for band in config["filters"]:
    records = {p["epoch"]: p for p in images if p["band"] == band}
    loaded = [load(records[e], k) for k in ["sci", "wht"] for e in ["older", "newer"]]
    old, new, w_old, w_new = [v[0] for v in loaded]
    good = np.isfinite(old+new+w_old+w_new) & (w_old > 0) & (w_new > 0)
    # Fill only for computing local derivatives. Erosion excludes the filled edges.
    old[~good], new[~good] = 0, 0
    difference, valid, calibration = fit_registration_residual(old, new, good)
    calibration["band"] = band
    calibration["max_header_registration_px"] = max(v[1] for v in loaded)
    calibrations.append(calibration)
    flux, score, noise = aperture_score(difference, w_old, w_new, valid, config["aperture_radius_arcsec"]/scale)
    arrays[band] = dict(older=old, newer=new, difference=difference, score=score, flux=flux,
                       noise=noise, valid=valid)
    for sign in [1, -1]:
        peaks = local_peaks(score, sign, config["peak_to_local_noise_min"])
        calibration[f"peaks_sign_{sign}"] = len(peaks)
        for rank, (y, x) in enumerate(peaks[:config["max_ranked_per_sign_and_band"]], 1):
            if min(y, x, side-y-1, side-x-1) < 16:
                continue
            ra, dec = grid.pixel_to_world_values(x+x0, y+y0)
            values = stamp_measurements(band, y, x, sign)
            candidates.append(dict(band=band, sign=sign, rank=rank, ra=float(ra), dec=float(dec),
                                   x=int(x), y=int(y), **values))
    print(band, calibration, flush=True)

table = pd.DataFrame(candidates)
table.to_csv(out / "band_peaks.csv", index=False)
blue, red = config["filters"][1], config["filters"][0]
matches = []
for sign in [1, -1]:
    first = table[(table.band == red) & (table.sign == sign)].reset_index(drop=True)
    second = table[(table.band == blue) & (table.sign == sign)].reset_index(drop=True)
    if first.empty or second.empty:
        continue
    distance, near = cKDTree(second[["x", "y"]]).query(first[["x", "y"]])
    for i, (d, j) in enumerate(zip(distance, near, strict=True)):
        if d*scale > config["match_radius_arcsec"]:
            continue
        a, b = first.iloc[i], second.iloc[j]
        if min(sign*a.score, sign*b.score) < config["aperture_to_blank_noise_min"]:
            continue
        x, y = int(a.x), int(a.y)
        morphology = [a.counter_flux_ratio, b.counter_flux_ratio]
        both_clean = bool(np.isfinite(morphology).all() and max(morphology) < 2
                          and min(a.annular_coverage, b.annular_coverage) >= .8)
        matches.append(dict(candidate=f"GS-JWST-{float(a.ra):.6f}{float(a.dec):+.6f}", sign=sign,
            ra=float(a.ra), dec=float(a.dec), x=x, y=y, separation_arcsec=float(d*scale),
            f444w_score=float(a.score), f356w_score=float(b.score),
            minimum_absolute_score=float(min(abs(a.score), abs(b.score))),
            compactness_gate=both_clean, f444w_counter_flux_ratio=float(a.counter_flux_ratio),
            f356w_counter_flux_ratio=float(b.counter_flux_ratio),
            f444w_delta_njy=float(a.delta_njy), f356w_delta_njy=float(b.delta_njy),
            f444w_old_njy=float(a.older_aperture_njy), f444w_new_njy=float(a.newer_aperture_njy),
            f356w_old_njy=float(b.older_aperture_njy), f356w_new_njy=float(b.newer_aperture_njy),
            disposition="REQUIRES_INDEPENDENT_EVIDENCE"))
matched = pd.DataFrame(matches).sort_values(["compactness_gate", "minimum_absolute_score"], ascending=False) if matches else pd.DataFrame()
matched.to_csv(out / "matched_peaks.csv", index=False)
common = np.isfinite(arrays[blue]["score"]+arrays[red]["score"])
write_json(out / "summary.json", dict(executed_utc=utcnow(), band_peaks=len(table), two_band_peaks=len(matches),
    positive_matches=sum(r["sign"] == 1 for r in matches), negative_matches=sum(r["sign"] == -1 for r in matches),
    compact_positive=sum(r["sign"] == 1 and r["compactness_gate"] for r in matches),
    compact_negative=sum(r["sign"] == -1 and r["compactness_gate"] for r in matches),
    searched_area_sq_arcmin=float(common.sum()*scale**2/3600), pixel_scale_arcsec=scale,
    calibration=calibrations, paired_by="Nearest same-sign peak, no position merging",
    classification="UNVALIDATED", products_sha256=text_digest(out / "products.json"),
    script_sha256=text_digest(__file__), diagnostic_module_sha256=text_digest(root / "src/cosmic_anomaly_lab/epoch_search.py")))

for sign in [1, -1]:
    if matched.empty:
        break
    chosen = matched[matched.sign == sign].head(5)
    if chosen.empty:
        continue
    figure = Figure(figsize=(12.5, 2.35*len(chosen)+.5), layout="constrained")
    axes = figure.subplots(len(chosen), 6, squeeze=False)
    for ax_group, row in zip(axes, chosen.itertuples(), strict=True):
        sl = np.s_[row.y-16:row.y+17, row.x-16:row.x+17]
        for k, band in enumerate(config["filters"]):
            a = arrays[band]
            old, new, diff = [np.where(a["valid"][sl], a[key][sl], np.nan)
                              for key in ["older", "newer", "difference"]]
            low, high = np.nanpercentile(np.r_[old.ravel(), new.ravel()], [5, 99.5])
            lim = float(np.nanpercentile(abs(diff), 99))
            for ax, image, label in zip(ax_group[3*k:3*k+3], [old, new, diff],
                    ["2025 observed", "2026 observed", "2026 minus model"], strict=True):
                ax.imshow(image, origin="lower", cmap="RdBu_r" if label.endswith("model") else "gray",
                          vmin=-lim if label.endswith("model") else low, vmax=lim if label.endswith("model") else high)
                ax.scatter([16], [16], facecolors="none", edgecolors="#3de2b9", s=150, linewidths=.8)
                ax.set(xticks=[], yticks=[], title=f"{band.split('-')[0]} · {label}")
            ax_group[0].set_ylabel(f"{row.ra:.6f}\n{row.dec:+.6f}", fontsize=8)
    figure.suptitle(("Brightening" if sign == 1 else "Fading")+
        " trial locations · actual JWST NIRCam visit cutouts\n1.32″ per panel · unvalidated difference peaks", fontsize=12)
    save_figure(figure, root / f"figures/jwst_epoch_search/{'brightening' if sign == 1 else 'fading'}_atlas")
print(matched.head(12).to_string(index=False))
