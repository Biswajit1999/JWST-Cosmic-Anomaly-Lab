"""Measure released eROSITA and ROSAT X-ray images and events at JADES 103607."""

import json
import warnings
from datetime import datetime
from io import BytesIO

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
from download_jades_103607_xray_surveys import BANDS
from download_jades_103607_xray_surveys import TAG as INPUT_TAG
from download_jwst_cross_detector import OUT
from matplotlib.colors import PowerNorm
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from measure_jwst_epoch_candidate import verified_path
from scipy.stats import binom, chi2
from unlzw3 import unlzw
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_xray_surveys"
PRIOR = OUT / "goods_s_prior_records.json"
ERS_RADIUS = 40.0
ROS_RADIUS = 120.0
ROS_ANNULUS = (180.0, 300.0)
ROS_PI = (11, 240)
CONFIDENCE = 0.9987


def image_wcs(header):
    # Astropy repairs the legacy DATE-OBS and RADECSYS cards in archived images.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        return WCS(header).celestial


def cutout(hdu, point, half_width):
    wcs = image_wcs(hdu.header)
    x, y = map(float, wcs.world_to_pixel(point))
    x0, x1 = max(0, int(x)-half_width), min(hdu.data.shape[1], int(x)+half_width+1)
    y0, y1 = max(0, int(y)-half_width), min(hdu.data.shape[0], int(y)+half_width+1)
    if x1 <= x0 or y1 <= y0:
        raise ValueError("Requested position is outside the published image")
    xx, yy = np.meshgrid(np.arange(x0, x1), np.arange(y0, y1))
    offset = wcs.pixel_to_world(xx, yy).transform_to(point.skyoffset_frame())
    east = offset.lon.wrap_at(180*u.deg).arcsec
    north = offset.lat.arcsec
    extent = [x0-.5, x1-.5, y0-.5, y1-.5]
    return dict(wcs=wcs, extent=extent, x=xx, y=yy, radius=np.hypot(east, north),
                data=np.array(hdu.data[y0:y1, x0:x1], float))


def count_limit(count):
    # Assign every observed source-aperture photon to the target for a conservative bound.
    return float(0.5 * chi2.ppf(CONFIDENCE, 2*(count+1)))


def verified_fits(product):
    path = verified_path(product)
    return fits.open(BytesIO(unlzw(path)) if path.suffix == ".Z" else path)


def e_time(header, seconds):
    if header["TIMEUNIT"] != "s" or header["TIMESYS"] != "TT":
        raise ValueError("Unexpected eROSITA clock definition")
    reference = Time(header["MJDREF"], format="mjd", scale="tt")
    return (reference + (np.asarray(seconds) + header["TIMEZERO"])*u.s).utc.isot.tolist()


def erosita_band(products, point):
    with verified_fits(products["counts"]) as image:
        counts = cutout(image[0], point, 45)
        h = image[0].header
        events = image["EVENTS"].data
        # The released band's post-GTI EVENTS ext maps to the count pixels via physical WCS.
        x = np.floor((events["X"] - h["CRVAL1P"])/h["CDELT1P"] + h["CRPIX1P"]-.5).astype(int)
        y = np.floor((events["Y"] - h["CRVAL2P"])/h["CDELT2P"] + h["CRPIX2P"]-.5).astype(int)
        left, right, bottom, top = counts["extent"]
        inside = (x > left) & (x < right) & (y > bottom) & (y < top)
        xx, yy = x[inside]-int(left+.5), y[inside]-int(bottom+.5)
        radius = counts["radius"][yy, xx]
        local_times = np.sort(events["TIME"][inside][radius <= 180])
        dates = np.unique(np.array(e_time(h, local_times), dtype="U23").astype("datetime64[D]"))
        within = radius <= ERS_RADIUS
        aperture = counts["radius"] <= ERS_RADIUS
        total = int(counts["data"][aperture].sum())
        if total != within.sum() or len(local_times) == 0:
            raise ValueError("eROSITA event and image aperture cross-check failed")
        timing = dict(
            tile_first_last_utc=e_time(h, [h["TSTART"], h["TSTOP"]]),
            local_check_radius_arcsec=180, local_event_count=len(local_times),
            local_first_last_event_utc=e_time(h, local_times[[0, -1]]),
            local_photon_utc_dates=dates.astype(str).tolist(),
            in_aperture_event_utc=e_time(h, np.sort(events["TIME"][inside][within])),
        )
    maps = {}
    with verified_fits(products["exposure"]) as exposure, verified_fits(products["background"]) as bg:
        for key, hdu in [("exposure", exposure[0]), ("background", bg[0])]:
            grid = cutout(hdu, point, 45)
            if grid["extent"] != counts["extent"] or not np.allclose(
                grid["radius"], counts["radius"], rtol=0, atol=0.01
            ):
                raise ValueError(f"eROSITA count and {key} maps have different grids")
            maps[key] = grid["data"]
    exp = float(maps["exposure"][aperture].mean())
    bkg = float(maps["background"][aperture].sum())
    if exp <= 0 or not np.isfinite(exp+bkg):
        raise ValueError("eROSITA aperture has invalid exposure or background")
    return dict(
        radius_arcsec=ERS_RADIUS, counts=total, expected_model_background_counts=bkg,
        signed_net_counts=total-bkg, mean_corrected_exposure_s=exp,
        conservative_aperture_rate_upper_ct_s=count_limit(total)/exp, timing=timing,
    ), counts


def erosita(products, point):
    limits = {}
    for name, meta in products["limits"].items():
        path = verified_path(meta)
        record = json.loads(path.read_text(encoding="utf-8"))
        if not record["de_sky"] or str(record["field"]).zfill(6) != products["tile_id"]:
            raise ValueError(f"Inconsistent eRODat limit cell for {name}")
        limits[name] = dict(
            band_keV=BANDS[record["band"]], aperture_counts=record.get("Count"),
            background_counts=record.get("Bkg_counts"),
            source_map_background_counts=record.get("Bkg_SourceMap"),
            exposure_s=record["Exposure"], smooth_background_upper_erg_cm2_s=record["UL_B"],
            source_map_upper_erg_cm2_s=record["UL_S"], catalog_position_flag=record["Flag_pos"],
        )
    catalogues = {}
    for name, meta in products["catalogues"].items():
        table = Table.read(verified_path(meta))
        if table.colnames[:1] != ["iauname"] or not {"ra", "dec", "pos_err"} <= set(table.colnames):
            raise ValueError(f"Unexpected eRODat cone result: {name}")
        catalogues[name] = [
            dict(name=str(row["iauname"]), icrs_deg=[float(row["ra"]), float(row["dec"])],
                 position_error_arcsec=float(row["pos_err"]),
                 separation_arcsec=float(point.separation(
                     SkyCoord(row["ra"], row["dec"], unit="deg")
                 ).arcsec))
            for row in table
        ]
    bands, images = {}, {}
    for band, meta in products["bands"].items():
        bands[band], images[band] = erosita_band(meta, point)
        bands[band]["energy_keV"] = BANDS[band]
    # The full calibrated list has looser time screening than the released band images.
    with verified_fits(products["events"]) as events:
        if events["EVENTS"].header["NAXIS2"] < max(
            row["timing"]["local_event_count"] for row in bands.values()
        ):
            raise ValueError("Unexpected eROSITA full-tile event inventory")
        event_rows = events["EVENTS"].header["NAXIS2"]
    return dict(tile_id=products["tile_id"], archive_limits=limits, cone_radius_arcmin=3,
                catalogue_matches=catalogues, dr1_maps=bands, dr1_full_tile_event_rows=event_rows), images


def rosat_events(hdu):
    header = fits.Header()
    for axis in (1, 2):
        for key, column_key in [
            ("CTYPE", "TCTYP"), ("CRPIX", "TCRPX"), ("CRVAL", "TCRVL"), ("CDELT", "TCDLT"),
        ]:
            header[f"{key}{axis}"] = hdu.header[f"{column_key}{axis}"]
        header[f"CUNIT{axis}"] = "deg"
    header["RADESYS"] = hdu.header["RADECSYS"]
    header["EQUINOX"] = hdu.header["EQUINOX"]
    data = hdu.data
    valid = (data["PI"] >= ROS_PI[0]) & (data["PI"] <= ROS_PI[1]) & (data["STATUS"] == 0)
    # STDEVT's integer X/Y values are FITS pixel coordinates with origin 1.
    coords = WCS(header).pixel_to_world(data["X"][valid]-1, data["Y"][valid]-1)
    return coords, data["TIME"][valid]


def rosat_observation(meta, point):
    with verified_fits(meta["merged_exposure"]) as exposure:
        h = exposure[0].header
        grid = cutout(exposure[0], point, 24)
        valid = np.isfinite(grid["data"]) & (grid["data"] > 0)
        r = grid["radius"]
        source = r <= ROS_RADIUS
        annulus = (r >= ROS_ANNULUS[0]) & (r <= ROS_ANNULUS[1])
        area_fraction = float(valid[source].mean())
        source_exp = float(grid["data"][source & valid].sum())
        annulus_exp = float(grid["data"][annulus & valid].sum())
        mean_exp = source_exp/source.sum()
        first = datetime.strptime(h["DATE-OBS"], "%d/%m/%y").date().isoformat()
        last = datetime.strptime(h["DATE_END"], "%d/%m/%y").date().isoformat()
    with verified_fits(meta["standard_events"]) as event_file:
        events = event_file["STDEVT"]
        coords, scc = rosat_events(events)
        offset = coords.transform_to(point.skyoffset_frame())
        radius = np.hypot(offset.lon.wrap_at(180*u.deg).arcsec, offset.lat.arcsec)
        local = radius <= ROS_ANNULUS[1]
        if local.any():
            ref = events.header["MJDREFI"] + events.header["MJDREFF"]
            nominal = Time(ref + np.array([scc[local].min(), scc[local].max()])/86400,
                           format="mjd", scale="utc")
            nominal_dates = [d[:10] for d in nominal.isot]
        else:
            nominal_dates = None
        x, y = grid["wcs"].world_to_pixel(coords)
        left, right, bottom, top = grid["extent"]
        pixels, _, _ = np.histogram2d(
            y, x, bins=[np.arange(bottom, top+1), np.arange(left, right+1)],
        )
        on_map = (x >= left) & (x < right) & (y >= bottom) & (y < top)
        ix = np.floor(x[on_map]-left).astype(int)
        iy = np.floor(y[on_map]-bottom).astype(int)
        usable = valid[iy, ix]
        radii = radius[on_map]
    measurement = dict(
        sequence=meta["sequence"], center_offset_deg=meta["center_distance_deg"],
        archive_field_dates_utc=[first, last],
        nominal_scc_dates_within_300_arcsec=nominal_dates, events_within_300_arcsec=int(local.sum()),
        positive_exposure_source_fraction=area_fraction,
        positive_exposure_annulus_fraction=float(valid[annulus].mean()),
        mean_aperture_exposure_s=float(mean_exp), counts=None, annulus_counts=None,
        expected_annulus_background_counts=None, signed_net_counts=None, conditional_p=None,
        conservative_aperture_rate_upper_ct_s=None,
    )
    grid["data"] = pixels
    if area_fraction < 1 or annulus_exp <= 0:
        return dict(status="no usable full source aperture", **measurement), grid
    n = int(np.count_nonzero(usable & (radii <= ROS_RADIUS)))
    m = int(np.count_nonzero(usable & (radii >= ROS_ANNULUS[0]) & (radii <= ROS_ANNULUS[1])))
    bkg = m*source_exp/annulus_exp
    p_value = float(binom.sf(n-1, n+m, source_exp/(source_exp+annulus_exp)))
    measurement.update(
        counts=n, annulus_counts=m, expected_annulus_background_counts=bkg, signed_net_counts=n-bkg,
        conditional_p=p_value, conservative_aperture_rate_upper_ct_s=count_limit(n)/mean_exp,
    )
    return dict(status="measured", **measurement), grid


def render(erosita_images, rosat_images, rows, point, neighbor):
    # The selected ROSAT panels are the all-sky exposure and the deepest qualifying pointing.
    by_seq = {row["sequence"]: row for row in rows if row["status"] == "measured"}
    survey, = [seq for seq in by_seq if seq.startswith("RS")]
    deepest = max((seq for seq in by_seq if seq.startswith("RP")),
                  key=lambda seq: by_seq[seq]["mean_aperture_exposure_s"])
    panels = [
        (erosita_images["024"], "eROSITA DR1 · 0.2–2.3 keV", ERS_RADIUS),
        (erosita_images["023"], "eROSITA DR1 · 2.3–5.0 keV", ERS_RADIUS),
        *((rosat_images[seq], f"ROSAT PSPC · {seq}", ROS_RADIUS) for seq in (survey, deepest)),
    ]
    fig = Figure(figsize=(10.9, 9.4))
    fig.subplots_adjust(left=.115, right=.906, bottom=.13, top=.89, wspace=.47, hspace=.29)
    for i, (image, title, radius) in enumerate(panels):
        wcs = image["wcs"]
        ax = fig.add_subplot(2, 2, i+1, projection=wcs)
        data = image["data"]
        plotted = ax.imshow(data, origin="lower", extent=image["extent"], interpolation="nearest",
                            cmap="Greys", norm=PowerNorm(.55, vmin=0, vmax=max(1, data.max())))
        x, y = wcs.world_to_pixel(point)
        nx, ny = wcs.world_to_pixel(neighbor)
        scale = float(wcs.proj_plane_pixel_scales()[0].to_value(u.arcsec))
        ax.add_patch(Circle((x, y), radius/scale, lw=1.35, edgecolor="#bb3262", facecolor="none"))
        ax.plot(x, y, "+", color="#bb3262", markersize=6)
        ax.plot(nx, ny, "x", color="#087e86", markersize=6)
        for coord in ax.coords:
            coord.set_ticks(number=4)
            coord.set_ticklabel(size=9)
        ax.coords[0].set_axislabel("Right ascension (J2000)" if i >= 2 else "")
        ax.coords[1].set_axislabel("Declination (J2000)" if i % 2 == 0 else "")
        ax.set_title(title, fontsize=10)
        bar = fig.colorbar(plotted, ax=ax, shrink=.79, ticks=range(int(plotted.norm.vmax)+1))
        bar.ax.tick_params(labelsize=9)
        if i % 2:
            bar.set_label("Observed counts per pixel", fontsize=9)
    fig.suptitle("Historical X-ray counts at the JADES 103607 position", fontsize=14, y=.965)
    fig.legend(handles=[
        Line2D([], [], color="#bb3262", marker="+", ls="", label="JWST position and aperture"),
        Line2D([], [], color="#087e86", marker="x", ls="", label="Nearest catalogued Chandra source"),
    ], loc="lower center", bbox_to_anchor=(.5, .045), ncol=2, frameon=False, fontsize=9)
    save_figure(fig, FIGURES / TAG)


if __name__ == "__main__":
    input_path = OUT / f"{INPUT_TAG}.json"
    inputs = json.loads(input_path.read_text(encoding="utf-8"))
    prior = json.loads(PRIOR.read_text(encoding="utf-8"))
    match, = [row["nearest"][0] for tag, row in
              prior["literature_catalogues"]["luo_2017"]["positions"].items()
              if tag.startswith("GS-001-B")]
    neighbor = SkyCoord(match["RAJ2000"], match["DEJ2000"], frame="fk5", unit="deg").icrs
    point = SkyCoord(*inputs["icrs_deg"], unit="deg")
    e, e_images = erosita(inputs["erosita"], point)
    r, r_images = [], {}
    for meta in inputs["rosat"]["sequences"]:
        row, image = rosat_observation(meta, point)
        r.append(row)
        r_images[row["sequence"]] = image
    render(e_images, r_images, r, point, neighbor)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), icrs_deg=inputs["icrs_deg"],
        script_sha256=text_digest(__file__), inputs_sha256=text_digest(input_path),
        chandra_prior_sha256=text_digest(PRIOR), upper_confidence=CONFIDENCE,
        nearest_chandra_icrs_deg=[neighbor.ra.deg, neighbor.dec.deg],
        nearest_chandra_separation_arcsec=point.separation(neighbor).arcsec,
        erosita=e, rosat=dict(source_radius_arcsec=ROS_RADIUS, background_annulus_arcsec=ROS_ANNULUS,
                            pi_channels_inclusive=ROS_PI, observations=r),
        limitations=[
            "eROSITA DR2 includes DR1 observations and yields an aggregate three-scan flux limit.",
            "Local photon dates describe event sampling. Per-visit sensitivity was not estimated.",
            "Measured count-rate limits conservatively assign every aperture count to the position. "
            "They omit encircled-energy and spectral-response corrections.",
            "ROSAT annular diagnostics assume background proportional to the exposure map and "
            "have no multiple-epoch correction.",
            "ROSAT local dates use the nominal spacecraft-clock reference without drift "
            "calibration and are reported at date precision.",
            "Both instruments' apertures encompass the neighbouring Chandra source.",
        ],
    ))
