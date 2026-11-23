"""Test the location of the October 2026 excess relative to the quiescent light."""

import json

import numpy as np
from astropy.coordinates import SkyCoord
from download_jwst_cross_detector import OUT
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import aperture, pixel_scale, show_stamp
from measure_jwst_module_b_brightening import ANNULUS, NEW, OLD, local_flux_reference
from scipy.ndimage import gaussian_filter, shift
from scipy.optimize import least_squares
from search_jwst_cross_detector import check_same_grid
from vet_jwst_field_b_candidate import FIGURES, vet_mosaic

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

BANDS = {"F150W": .18, "F200W": .20, "F277W": .23}
BASELINES = ("jw01180010001", "jw01180012001", OLD)
COMMON_SMOOTHING = (0., .03)


def patch(image, x, y, scale, radius):
    half = int(np.ceil(radius/scale)) + 1
    ix, iy = round(x), round(y)
    yy, xx = np.mgrid[iy-half:iy+half+1, ix-half:ix+half+1]
    return (xx-x)*scale, (yy-y)*scale, image[yy, xx]


def gaussian(x, y, amplitude, dx, dy, sigma):
    return amplitude*np.exp(-((x-dx)**2+(y-dy)**2)/(2*sigma**2))


def fit_core(x, y, image, radius):
    valid = (np.hypot(x, y) <= radius) & np.isfinite(image)
    if valid.sum() < .95*(np.hypot(x, y) <= radius).sum():
        return None
    xx, yy, values = x[valid], y[valid], image[valid]
    background = np.median(values[np.hypot(xx, yy) >= .75*radius])
    peak = max(values.max()-background, 1e-6)

    def residual(p):
        amplitude, dx, dy, sigma, sky = p
        return gaussian(xx, yy, amplitude, dx, dy, sigma) + sky - values

    centre_limit = .65*radius
    result = least_squares(
        residual, (peak, 0., 0., .30*radius, background),
        bounds=((0., -centre_limit, -centre_limit, .012, -np.inf),
                (np.inf, centre_limit, centre_limit, .78*radius, np.inf)),
    )
    if not result.success or np.any(result.active_mask):
        return None
    amplitude, dx, dy, sigma, sky = map(float, result.x)
    return dict(amplitude_njy_per_pixel=amplitude, offset_xy_arcsec=[dx, dy],
                core_sigma_arcsec=sigma, background_njy_per_pixel=sky)


def local_alignment(old, new, x, y, scale, radius):
    valid = np.isfinite(old) & np.isfinite(new)
    photometric = local_flux_reference(old, new, valid, x, y, scale)
    references = []
    for ref in photometric["references"]:
        if not .85 <= ref["new_to_old"] <= 1.15:
            continue
        rx, ry = np.array([x, y]) + np.array(ref["offset_arcsec"])/scale
        old_core, new_core = (fit_core(*patch(image, rx, ry, scale, radius), radius) for image in (old, new))
        if old_core is None or new_core is None:
            continue
        displacement = np.subtract(new_core["offset_xy_arcsec"], old_core["offset_xy_arcsec"])
        if np.linalg.norm(displacement) <= .08:
            references.append(dict(
                **ref, new_minus_old_centroid_xy_arcsec=displacement.tolist(),
                new_minus_old_variance_arcsec2=new_core["core_sigma_arcsec"]**2-old_core["core_sigma_arcsec"]**2,
            ))
    if len(references) < 6:
        raise ValueError(f"Insufficient local core matches ({len(references)})")
    displacements = np.array([r["new_minus_old_centroid_xy_arcsec"] for r in references])
    old_to_new = np.median(displacements, axis=0)
    offsets = np.array([v-np.median(np.delete(displacements, j, axis=0), axis=0)
                        for j, v in enumerate(displacements)])
    deltas = [r["new_minus_old_variance_arcsec2"] for r in references]
    return dict(reference_count=len(references), reference_sources=references,
                old_to_new_xy_arcsec=old_to_new.tolist(),
                held_out_reference_radius68_arcsec=float(np.quantile(np.linalg.norm(offsets, axis=1), .68)),
                matched_variance_arcsec2=float(np.median(deltas)),
                new_to_old_scale=photometric["median_gain"])


def convolve(image, displacement_yx, sigma):
    good = np.isfinite(image)
    clean = np.where(good, image, 0.)
    clean = shift(clean, displacement_yx, order=3, mode="constant", prefilter=True)
    mask = shift(good.astype(float), displacement_yx, order=1, mode="constant")
    if sigma > .10:
        clean = gaussian_filter(clean, sigma, mode="constant")
        mask = gaussian_filter(mask, sigma, mode="constant")
    clean[mask < .999] = np.nan
    return clean


def match_pair(old, new, alignment, scale, extra_blur):
    sigma2 = alignment["matched_variance_arcsec2"]
    o = convolve(old, (0., 0.), np.sqrt(max(sigma2, 0.) + extra_blur**2)/scale)
    shift_yx = -np.array(alignment["old_to_new_xy_arcsec"][::-1])/scale
    n = convolve(new, shift_yx, np.sqrt(max(-sigma2, 0.) + extra_blur**2)/scale)
    return o, n, n-alignment["new_to_old_scale"]*o


def injected_recoveries(difference, x, y, scale, radius, source):
    """Recover a synthetic positive core on independent observed residual patches."""
    offsets = [(dx, dy) for dy in range(-7, 8) for dx in range(-7, 8) if 2.0 <= np.hypot(dx, dy) <= 7.]
    injected = ((0., 0.), (.04, 0.), (-.04, 0.), (0., .04), (0., -.04))
    measured = []
    for dx, dy in offsets:
        xx, yy, observed = patch(difference, x+dx/scale, y+dy/scale, scale, radius)
        if not np.isfinite(observed).all():
            continue
        for bx, by in injected:
            signal = gaussian(xx, yy, source["amplitude_njy_per_pixel"], bx, by, source["core_sigma_arcsec"])
            fit = fit_core(xx, yy, observed+signal, radius)
            if fit is not None:
                measured.append(dict(control_position_offset_xy_arcsec=[dx, dy], injected_xy_arcsec=[bx, by],
                                     recovered_xy_arcsec=fit["offset_xy_arcsec"]))
    if len(measured) < 100:
        raise ValueError(f"Too few recovered mock core placements ({len(measured)})")
    errors = np.array([np.subtract(v["recovered_xy_arcsec"], v["injected_xy_arcsec"]) for v in measured])
    norms = np.linalg.norm(errors, axis=1)
    return dict(placements=len(offsets)*len(injected), usable_recoveries=len(measured),
                median_xy_error_arcsec=np.median(errors, axis=0).tolist(),
                radius68_arcsec=float(np.quantile(norms, .68)), radius95_arcsec=float(np.quantile(norms, .95)),
                measurements=measured)


def sky_offsets(wcs, origin, px, py, scale, host, event):
    (hx, hy), (ex, ey) = [np.array([px, py])+np.array(v["offset_xy_arcsec"])/scale for v in (host, event)]
    coords = wcs.pixel_to_world(np.array([hx, ex])+origin[0], np.array([hy, ey])+origin[1])
    east, north = coords[0].spherical_offsets_to(coords[1])
    return [float(east.arcsec), float(north.arcsec)]


def prepare_pair(old, new, radius):
    o, ow, ox, oy, ox0, oy0 = old
    n, nw, _, _, nx0, ny0 = new
    scale = pixel_scale(ow)
    check_same_grid((o, ow, (ox0, oy0), scale), (n, nw, (nx0, ny0), pixel_scale(nw)))
    alignment = local_alignment(o, n, ox, oy, scale, radius)
    return o, n, ow, ox, oy, (ox0, oy0), scale, alignment


def measure_pair(old, new, radius):
    o, n, ow, ox, oy, origin, scale, alignment = prepare_pair(old, new, radius)
    results, views = [], []
    for blur in COMMON_SMOOTHING:
        template, _, difference = match_pair(o, n, alignment, scale, blur)
        host, event = (fit_core(*patch(im, ox, oy, scale, radius), radius) for im in (template, difference))
        if host is None or event is None:
            raise ValueError("Unconstrained Gaussian core fit")
        sky = sky_offsets(ow, origin, ox, oy, scale, host, event)
        controls = injected_recoveries(difference, ox, oy, scale, radius, event) if blur == 0. else None
        results.append(dict(common_smoothing_sigma_arcsec=blur, host_core=host, event_core=event,
                            excess_minus_host_east_north_arcsec=sky, projected_offset_arcsec=float(np.linalg.norm(sky)),
                            excess_aperture_njy=aperture(difference, ox, oy, scale, .24, ANNULUS),
                            injected_position_trials=controls))
        views.append((template, difference, ox, oy, scale))
    return dict(registration=alignment, cases=results), views[0]


def inject_at_host(old, new, radius, event):
    """Place a centered mock brightening over measured pre-2026 host residuals."""
    o, n, wcs, ox, oy, origin, scale, alignment = prepare_pair(old, new, radius)
    template, _, difference = match_pair(o, n, alignment, scale, 0.)
    host = fit_core(*patch(template, ox, oy, scale, radius), radius)
    if host is None:
        raise ValueError("Unconstrained host fit in the quiescent comparison")
    x, y, sky = patch(difference, ox, oy, scale, radius)
    hx, hy = host["offset_xy_arcsec"]
    fits = []
    for strength in (.8, 1., 1.2):
        for width in (.8, 1., 1.2):
            synthetic = gaussian(x, y, strength*event["amplitude_njy_per_pixel"], hx, hy,
                                 width*event["core_sigma_arcsec"])
            recovered = fit_core(x, y, sky+synthetic, radius)
            if recovered is None:
                raise ValueError("Unconstrained injected fit in the quiescent host")
            offset = sky_offsets(wcs, origin, ox, oy, scale, host, recovered)
            fits.append(dict(relative_peak=strength, relative_sigma=width, recovered_east_north_arcsec=offset,
                             projected_offset_arcsec=float(np.linalg.norm(offset))))
    return dict(registration=alignment, difference_aperture_njy=aperture(difference, ox, oy, scale, .24, ANNULUS),
                injections=fits)


def plot(results, views, date_labels):
    figure = Figure(figsize=(10.7, 8.0), layout="constrained")
    axes = figure.subplots(3, 3, width_ratios=(1, 1, 1.4))
    for i, band in enumerate(BANDS):
        chosen = (BASELINES[1], band)
        template, diff, x, y, scale = views[chosen]
        result = results[chosen]["cases"][0]
        for ax, image, title, cmap in zip(axes[i, :2], (template, diff),
                                         ("2023 baseline", "2026 minus 2023"), ("gray", "RdBu_r"), strict=True):
            local = patch(image, x, y, scale, .38)[2]
            peak = np.nanpercentile(abs(local), 98)
            bounds = (np.nanpercentile(local, 5), peak) if cmap == "gray" else (-peak, peak)
            show_stamp(ax, image, x, y, scale, bounds, cmap, .16, half_arcsec=.38)
            ax.set_title(f"{band} · {title}", fontsize=10)
            for field, colour, mark in (("host_core", "#36b9df", "+"), ("event_core", "#d81387", "x")):
                ax.scatter(*result[field]["offset_xy_arcsec"], s=56, marker=mark, color=colour, linewidths=1.4)
        ax = axes[i, 2]
        for visit, colour in zip(BASELINES, ("#355a96", "#268069", "#ac4827"), strict=True):
            r = results[visit, band]["cases"]
            xy = np.array([v["excess_minus_host_east_north_arcsec"] for v in r])*1000
            ax.plot(*xy.T, color=colour, alpha=.7, lw=1)
            ax.scatter(*xy[0], color=colour, label=date_labels[visit])
            ax.scatter(*xy[1], facecolors="none", edgecolors=colour)
        radius95 = max(results[v, band]["cases"][0]["injected_position_trials"]["radius95_arcsec"] for v in BASELINES)
        lim = max(.05, radius95)*1000
        ax.axhline(0, color=".7", lw=.7)
        ax.axvline(0, color=".7", lw=.7)
        ax.set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="East offset (mas)", ylabel="North offset (mas)")
        ax.set_aspect("equal")
        if i == 0:
            ax.legend(loc="upper right", fontsize=7, frameon=False)
    figure.suptitle("JADES 103607 · Position of the October 2026 excess\n"
                    "Cyan + quiescent core · magenta × excess · open circles add common smoothing", fontsize=12)
    save_figure(figure, FIGURES/"module_b_001_excess_position")


def main():
    photo_path, vet_path, search_path = (OUT/name for name in (
        "module_b_001_six_filter_change.json", "module_b_vetting.json", "module_b_products.json"))
    photo, vet, search = [json.loads(path.read_text()) for path in (photo_path, vet_path, search_path)]
    (lead,) = (row for row in vet["candidates"] if row["sign"] > 0)
    sky = photo["icrs_deg"]
    assert SkyCoord(*sky, unit="deg").separation(SkyCoord(*lead["icrs_deg"], unit="deg")).arcsec < 1e-6
    products = [
        *photo["additional_products"], *lead["extra_products"],
        *[r for r in search["images"] if r["pointing"] == "001" and r["detector"] == "NRCB4"],
    ]
    chosen = {(r["visit"], r["band"]): r for r in products if r["visit"] in (*BASELINES, NEW) and r["band"] in BANDS}
    if len(chosen) != 12:
        raise ValueError(f"Expected three reference visits and one bright visit per filter, got {len(chosen)}")
    mosaics = {key: vet_mosaic(r, sky)[0] for key, r in chosen.items()}
    results, views, host_controls = {}, {}, []
    for band, radius in BANDS.items():
        for visit in BASELINES:
            print(band, visit, flush=True)
            results[visit, band], views[visit, band] = measure_pair(mosaics[visit, band], mosaics[NEW, band], radius)
        event_core = results[BASELINES[1], band]["cases"][0]["event_core"]
        for first, second in zip(BASELINES[:-1], BASELINES[1:], strict=True):
            host_controls.append(dict(band=band, older_visit=first, newer_visit=second,
                                      **inject_at_host(mosaics[first, band], mosaics[second, band], radius, event_core)))
    plot(results, views, {v: chosen[v, "F200W"]["dateobs"][:4] for v in BASELINES})
    write_json(OUT/"module_b_001_excess_position.json", dict(
        executed_utc=utcnow(), candidate=photo["candidate"], candidate_icrs_deg=sky, bright_visit=NEW,
        measurements=[dict(band=band, baseline_visit=visit, core_fit_radius_arcsec=BANDS[band], **data)
                      for (visit, band), data in results.items()],
        mosaic_products=list(chosen.values()), quiescent_host_controls=host_controls,
        six_filter_sha256=text_digest(photo_path), vetting_sha256=text_digest(vet_path),
        search_products_sha256=text_digest(search_path), script_sha256=text_digest(__file__),
        scope="Gaussian cores describe the light distribution. The three comparisons share the bright images."
              " The synthetic tests sample noise in field residuals and a centered Gaussian on quiescent host"
              " residuals. The location of the dynamical nucleus and chromatic PSF effects remain unmeasured.",
    ))


if __name__ == "__main__":
    main()
