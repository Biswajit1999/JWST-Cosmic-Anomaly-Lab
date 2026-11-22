"""Cross-match the transient with public bright-epoch MAST slitless extractions."""

import json
from io import StringIO
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
from astropy.io import fits
from astropy.table import Table
from download_jwst_cross_detector import API, OUT, ROOT, download
from download_jwst_field_b_candidate import STORAGE_MIB
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import pixel_scale, show_stamp
from vet_jwst_field_b_candidate import FIGURES, vet_mosaic

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

MAST = "https://mast.stsci.edu/api/v0.1/Download/file?"
VISITS = ("jw04549005001", "jw04549006001")
RADIUS = .25


def mast_product(cache, name, kind):
    url = MAST + urlencode(dict(uri="mast:JWST/product/"+name))
    path = cache.fetch(url, "cross_detector/"+name, survey="MAST JWST",
                       version="public GO 4549 slitless reduction", kind=kind)
    return path, dict(filename=name, url=url, path=path.relative_to(ROOT).as_posix(), sha256=digest(path))


def spectrum_position(path, sky):
    with fits.open(path) as h:
        header = h[0].header
        if header["PROGRAM"] != "04549" or header["PUPIL"] != "GRISMR" or header["DETECTOR"] != "NRCALONG":
            raise ValueError(f"Unexpected WFSS instrument settings in {path.name}")
        rows = h["EXTRACT1D"].data
        finite = np.isfinite(rows["SOURCE_RA"]+rows["SOURCE_DEC"])
        coords = SkyCoord(rows["SOURCE_RA"][finite], rows["SOURCE_DEC"][finite], unit="deg")
        separation = sky.separation(coords).arcsec
        i = int(np.argmin(separation))
        source = rows[np.flatnonzero(finite)[i]]
        good = (np.isfinite(source["WAVELENGTH"]+source["FLUX"]+source["FLUX_ERROR"])
                & (source["WAVELENGTH"] > 0) & (source["FLUX_ERROR"] > 0) & (source["DQ"] == 0))
        return dict(source_rows=len(rows), positioned_rows=int(finite.sum()), date_obs=header["DATE-OBS"],
                    first_exposure_mjd=header["EXPSTART"], source_catalogue=Path(header["SCATFILE"]).name,
                    matches_within_radius=int((separation <= RADIUS).sum()),
                    nearest=dict(id=int(source["SOURCE_ID"]), icrs_deg=[source["SOURCE_RA"], source["SOURCE_DEC"]],
                                 separation_arcsec=float(separation[i]), usable_spectral_pixels=int(good.sum())))


def catalogue_position(path, sky):
    table = Table.read(path)
    coords = table["sky_centroid"]
    finite = np.isfinite(coords.ra.deg+coords.dec.deg)
    values = sky.separation(coords[finite]).arcsec
    i = int(np.argmin(values))
    row = table[np.flatnonzero(finite)[i]]
    return dict(source_rows=len(table), positioned_rows=int(finite.sum()),
                matches_within_radius=int((values <= RADIUS).sum()),
                nearest=dict(id=int(row["label"]), icrs_deg=[row["sky_centroid"].ra.deg, row["sky_centroid"].dec.deg],
                             separation_arcsec=float(values[i])))


def optical_neighbours(cache, sky):
    url = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?" + urlencode({
        "-source": "J/MNRAS/425/2116/table2", "-c": f"{sky.ra.deg:.9f} {sky.dec.deg:.9f}",
        "-c.rs": 5, "-out.max": 25, "-out.all": 1,
    })
    receipt = download(cache, url, "field_b_transient_aces.tsv", "published_optical_redshift_catalogue")
    table = pd.read_csv(StringIO((ROOT / receipt["path"]).read_text(encoding="utf-8")), sep="\t", comment="#")
    for col in ("Seq", "z", "zhel", "Q", "RAJ2000", "DEJ2000", "MJD", "_r"):
        table[col] = pd.to_numeric(table[col], errors="coerce")
    secure = table[(table.Q >= 3) & table.zhel.notna()].copy()
    estimates = []
    for identifier, rows in secure.groupby("Seq"):
        point = rows.iloc[0]
        estimates.append(dict(
            aces_id=int(identifier), ra_j2000_deg=float(point.RAJ2000), dec_j2000_deg=float(point.DEJ2000),
            separation_arcsec=float(point["_r"]), median_heliocentric_redshift=float(rows.zhel.median()),
            observations=[dict(mask=r.Mask.strip(), mjd=r.MJD, z=r.z, z_heliocentric=r.zhel, quality=int(r.Q))
                          for r in rows.itertuples()],
        ))
    return receipt, sorted(estimates, key=lambda v: v["separation_arcsec"])


def conditional_luminosity(source, closest):
    """Report the rest-frame shifted-band absolute magnitude conditional on association."""
    base, bright = (source["measurements"][i]["apertures"][-1] for i in (0, 1))
    excess = bright["flux_njy"]-base["flux_njy"]
    sigma = np.hypot(bright["local_mad_sigma_njy"], base["local_mad_sigma_njy"])
    z = closest["median_heliocentric_redshift"]
    m_ab = float(31.4-2.5*np.log10(excess))
    return dict(conditional_host_aces_id=closest["aces_id"], redshift=z, bright_minus_baseline_njy=excess,
                local_aperture_scatter_njy=float(sigma), observed_aperture_ab_mag=m_ab,
                luminosity_distance_mpc=float(Planck18.luminosity_distance(z).value),
                approximate_rest_wavelength_um=1.50/(1+z),
                shifted_passband_absolute_ab_mag=float(m_ab-Planck18.distmod(z).value+2.5*np.log10(1+z)),
                assumptions="Planck18, redshifted same passband, no extinction or aperture correction; host unconfirmed")


def plot_neighbours(source, sky, neighbours):
    figure = Figure(figsize=(8.0, 4.6), layout="constrained")
    for ax, index in zip(figure.subplots(1, 2), (1, 5), strict=True):
        product = source["products"][index]
        (im, wcs, x, y, x0, y0), _ = vet_mosaic(product, (sky.ra.deg, sky.dec.deg))
        scale = pixel_scale(wcs)
        show_stamp(ax, np.arcsinh(im/.32), x, y, scale, (-.20, 4), "gray", .24, half_arcsec=2.7)
        for n in neighbours:
            coord = SkyCoord(n["ra_j2000_deg"], n["dec_j2000_deg"], unit="deg", frame="fk5").icrs
            gx, gy = wcs.world_to_pixel_values(coord.ra.deg, coord.dec.deg)
            px, py = (gx-x0-x)*scale, (gy-y0-y)*scale
            ax.scatter(px, py, marker="x", c="#f7bd6a", s=58, linewidths=1.5)
            ax.annotate(f"ACES {n['aces_id']}", (px, py), xytext=(5, 6 if py >= 0 else -16),
                        textcoords="offset points", color="#ffd38a", fontsize=8)
        ax.set_title(product["dateobs"][:10]+" · F150W", fontsize=10)
    figure.suptitle("GOODS-S transient and nearby galaxies\n"
                    "Observed JWST images · cyan = transient aperture · gold = ACES positions", fontsize=12)
    save_figure(figure, FIGURES / "field_b_transient_environment")


def main():
    source_path = OUT / "field_b_transient_f150w.json"
    source = json.loads(source_path.read_text())
    ra, dec = source["position_degrees"]
    sky = SkyCoord(ra, dec, unit="deg")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    url = API + "/exposures?" + urlencode(dict(ra=f"{ra:.9f}", dec=f"{dec:.9f}", arcmin=".03",
                                                 filters="F356W-GRISMR", output="csv"))
    listing = download(cache, url, "field_b_transient_wfss_exposures.csv", "archival_slitless_exposure_listing")
    frame = pd.read_csv(ROOT / listing["path"])
    selected = frame[(frame.extension == "rate") & frame.file.str.startswith(VISITS)
                     & (frame.detector == "NRCALONG") & (frame["filter"] == "F356W-GRISMR")]
    if len(selected) != 11 or selected.file.nunique() != 11:
        raise ValueError("The pinned WFSS search expects 11 exposure footprints")
    spectra = []
    for row in selected.sort_values("expstart").itertuples():
        path, product = mast_product(cache, row.file+"_x1d.fits", "observed_native_wfss_spectrum")
        evidence = spectrum_position(path, sky)
        spectra.append({**product, **evidence})
        print(product["filename"], evidence["matches_within_radius"], evidence["nearest"]["separation_arcsec"], flush=True)
    catalogues = []
    for name in sorted({r["source_catalogue"] for r in spectra}):
        path, product = mast_product(cache, name, "observed_wfss_source_catalogue")
        catalogues.append({**product, **catalogue_position(path, sky)})
    aces, neighbours = optical_neighbours(cache, sky)
    plot_neighbours(source, sky, neighbours)
    write_json(OUT / "field_b_transient_wfss.json", {
        "executed_utc": utcnow(), "candidate": source["candidate"], "position_degrees": [ra, dec],
        "crossmatch_radius_arcsec": RADIUS, "exposure_footprint_listing": listing,
        "x1d_exposures": spectra, "direct_image_catalogues": catalogues,
        "aces_query": aces, "aces_nearby_galaxies": neighbours,
        "conditional_luminosity": conditional_luminosity(source, neighbours[0]) if neighbours else None,
        "scope": "Pipeline source positions in archived X1D and direct catalogues; no forced extraction or spectral limit",
        "f150w_input_sha256": text_digest(source_path), "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()
