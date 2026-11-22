"""Compare the GOODS-S transient with archival Chandra, Hubble, VLA and MUSE data."""

import gzip
import json
import tarfile
from io import BytesIO
from urllib.parse import urlencode

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from download_jwst_cross_detector import API, OUT, ROOT, download
from download_jwst_transient_grism_geometry import STORAGE_MIB, fetch
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from measure_jwst_epoch_candidate import photometry, pixel_scale, show_stamp, verified_path
from measure_jwst_transient_grism_geometry import direct_cutout
from scipy.ndimage import map_coordinates
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

BUNDLE = "luo2017-cdfs-data-products03.tar.gz"
BUNDLE_MD5 = "md5:1d40972eded449f6f921482ffe82f4c9"
SOFT_IMAGE = "CDFS-7Ms-0p5to2-asca-im-bin1-01.fits.gz"
SOFT_BACKGROUND = "CDFS-7Ms-0p5to2-bin1-bkg-01.fits.gz"
SOFT_EXPOSURE = "CDFS-7Ms-0p5to2-bin1-02.emap.fits.gz"
MUSE_NAMES = ("musewide_field_mapping.fits", "expmap_dr1_seconds.fits")
RADIO_RADIUS = 20
BRIGHT = "2025-09-30"


def records(cache, sky):
    ra, dec = f"{sky.ra.deg:.10f}", f"{sky.dec.deg:.10f}"
    url = f"https://zenodo.org/records/8021940/files/{BUNDLE}?download=1"
    path = cache.fetch(url, "cross_detector/"+BUNDLE, survey="Chandra CDF-S",
                       version="7 Ms Luo et al. 2017; Zenodo 8021940", kind="published_mosaic_and_calibration",
                       checksum=BUNDLE_MD5)
    bundle = dict(url=url, path=path.relative_to(ROOT).as_posix(), sha256=digest(path), checksum=BUNDLE_MD5)
    cda = "https://cda.cfa.harvard.edu/fps-dbdal-rest/jax-rs/main/footprints?"
    # The public service also includes stacked/duplicate rows, handled by instrument and ObsId below.
    footprints_url = cda+urlencode(dict(POS=f"{ra},{dec}", SIZE="0.0", inst="ACIS-I,ACIS-S"))
    footprints = fetch(cache, footprints_url, "transient_chandra_footprints.xml", "observation_footprints",
                       survey="Chandra Data Archive", version="public ACIS")
    query = f"SELECT * FROM vlaecdfs1p4 WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{RADIO_RADIUS/3600:.8f}))"
    radio_url = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync?"+urlencode(dict(
        REQUEST="doQuery", LANG="ADQL", FORMAT="votable", QUERY=query, MAXREC=200))
    radio = fetch(cache, radio_url, "transient_vla_1p4.xml", "published_radio_catalogue_query",
                  survey="VLA E-CDF-S", version="Miller et al. 2013 through HEASARC")
    hst_url = API+"/thumb?"+urlencode(dict(ra=ra, dec=dec, size=14, all_filters=False,
                                           filters="f606w,f814w", output="fits_weight"))
    hst = download(cache, hst_url, "transient_hst_optical.fits", "observed_optical_cutout")
    maps = [fetch(cache, "https://musewide.aip.de/files/dr1/"+name, name, "coverage_map",
                  survey="VLT MUSE-Wide", version="DR1") for name in MUSE_NAMES]
    return dict(chandra_bundle=bundle, chandra_footprints=footprints, radio_query=radio, optical_cutout=hst,
                muse_maps=maps)


def pixels(hdu, sky, half=24):
    full = WCS(hdu.header).celestial
    x, y = map(float, full.world_to_pixel(sky))
    x0, y0 = round(x)-half, round(y)-half
    if min(x0, y0) < 0 or x0+2*half >= hdu.shape[1] or y0+2*half >= hdu.shape[0]:
        raise ValueError("Required local cutout falls outside the image")
    region = np.s_[y0:y0+2*half+1, x0:x0+2*half+1]
    wcs = full.slice(region)
    return hdu.data[region].copy(), wcs, (x-x0, y-y0)


def unpack_chandra(receipt, sky):
    wanted = {SOFT_IMAGE, SOFT_BACKGROUND, SOFT_EXPOSURE, "tbl04.fits", "tbl05.fits"}
    images, tables = {}, {}
    with tarfile.open(verified_path(receipt), "r|gz") as archive:
        for member in archive:
            name = member.name.rsplit("/", 1)[-1]
            if name not in wanted or not member.isfile():
                continue
            with archive.extractfile(member) as stream:
                data = stream.read()
            with fits.open(BytesIO(gzip.decompress(data) if name.endswith(".gz") else data)) as h:
                if name.endswith(".gz"):
                    images[name] = pixels(h[0], sky)
                else:
                    tables[name] = Table(h[1].data)
    if set(images)|set(tables) != wanted or [len(tables[k]) for k in ("tbl04.fits", "tbl05.fits")] != [1008, 47]:
        raise ValueError("Missing or incomplete published CDF-S products")
    source, wcs, xy = images[SOFT_IMAGE]
    for name, (data, other_wcs, other_xy) in images.items():
        if data.shape != source.shape or np.linalg.norm(np.subtract(other_xy, xy)) > .01:
            raise ValueError(f"Chandra mosaic grid differs in {name}")
        if np.max(abs(wcs.pixel_scale_matrix-other_wcs.pixel_scale_matrix)) > 1e-8:
            raise ValueError(f"Chandra pixel scales differ in {name}")
    return images, tables


def xray_sources(tables, sky, galaxy):
    result = []
    for table, rows in tables.items():
        coords = SkyCoord(rows["RA"], rows["DEC"], unit="deg")
        distance = sky.separation(coords).arcsec
        for i in np.argsort(distance)[:3]:
            row = rows[i]
            counterpart = str(row["CP_CAT"]).strip()
            cp_xy = ([float(row[f"{counterpart}_{key}"]) for key in ("RA", "DEC")]
                     if counterpart in ("WFI", "GOODSS", "GEMS", "CANDELS", "TENIS", "SEDS", "VLA") else None)
            cp_dist = galaxy.separation(SkyCoord(*cp_xy, unit="deg")).arcsec if cp_xy else None
            result.append(dict(
                table=table, id=int(row["ID"]), icrs=[coords[i].ra.deg, coords[i].dec.deg],
                separation_from_transient_arcsec=float(distance[i]),
                separation_from_grism_galaxy_arcsec=float(galaxy.separation(coords[i]).arcsec),
                positional_error_radius_68_arcsec=float(row["SIGMAX"]),
                counterpart_catalogue=counterpart, counterpart_icrs=cp_xy,
                counterpart_separation_from_grism_galaxy_arcsec=float(cp_dist) if cp_dist is not None else None,
                redshift=float(row["ZFINAL"]) if row["ZFINAL"] > 0 else None, published_type=str(row["TYPE"]).strip(),
                soft_flux_signed_erg_per_cm2_s=float(row["SB_FLUX"]),
            ))
    return sorted(result, key=lambda row: row["separation_from_transient_arcsec"])


def check_xray_pixels(images, positions):
    image, wcs, _ = images[SOFT_IMAGE]
    bkg = images[SOFT_BACKGROUND][0]
    exp = images[SOFT_EXPOSURE][0]
    y, x = np.indices(image.shape)
    target = {}
    for label, sky in positions.items():
        sx, sy = wcs.world_to_pixel(sky)
        radius = np.hypot(x-sx, y-sy)*pixel_scale(wcs)
        target[label] = []
        for r in (1., 1.5):
            mask = radius <= r
            target[label].append(dict(radius_arcsec=r, pixel_count=int(mask.sum()), observed_counts=int(image[mask].sum()),
                                      background_map_counts=int(bkg[mask].sum()),
                                      median_effective_exposure_cm2_s=float(np.median(exp[mask]))))
    return target


def optical_stamps(receipt, sky):
    results, figures = {}, {}
    with fits.open(verified_path(receipt)) as h:
        for band in ("F606W", "F814W"):
            sci, weight = h[band, "SCI"], h[band, "WHT"]
            cal = float(sci.header["PHOTFNU"])*1e9
            if not np.isclose(cal, 10) or sci.header["INSTRUME"] != "ACS":
                raise ValueError(f"Unexpected calibration or instrument for {band}")
            wcs = WCS(sci.header).celestial
            xy = tuple(float(v) for v in wcs.world_to_pixel(sky))
            image = np.where(weight.data > 0, sci.data*cal, np.nan)
            results[band] = dict(apertures=photometry(image, *xy, pixel_scale(wcs), (.12, .18, .24), (.30, .48)),
                                 exposure_members=[sci.header[k] for k in sci.header if k.startswith("FLT")])
            figures[band] = (image, wcs, xy)
    return results, figures


def plot(optical, direct, chandra, positions, sigma, date):
    fig = Figure(figsize=(10.2, 3.9), layout="constrained")
    axes = fig.subplots(1, 3)
    specs = ((optical, "Hubble ACS F814W", "gray", (-.18, 4), .9),
             (direct, f"JWST F356W · {date}", "gray", (-.18, 4.6), .4),
             (chandra, "Chandra · 0.5–2 keV", "magma", (0, 4), None))
    for ax, (stamp, title, cmap, limits, stretch) in zip(axes, specs, strict=True):
        values, wcs, (x, y) = stamp
        scale = pixel_scale(wcs)
        show_stamp(ax, np.arcsinh(values/stretch) if stretch else values, x, y, scale, limits, cmap, .24, 5)
        for label, sky, color, marker in positions:
            u, v = (np.array(wcs.world_to_pixel(sky))-[x, y])*scale
            ax.scatter(u, v, marker=marker, color=color, s=38, linewidths=1.2, label=label)
            if label == "Chandra 726":
                ax.add_patch(Circle((u, v), sigma, fill=False, linewidth=1, color=color))
        ax.set(xlim=(-5, 5), ylim=(-5, 5), title=title)
    handles, labels = axes[0].get_legend_handles_labels()
    handles.insert(0, axes[0].patches[0])
    fig.legend(handles, ["Transient aperture", *labels], loc="outside lower center", ncols=4, fontsize=8)
    fig.suptitle("GOODS-S transient and archival counterparts\n"
                 "Observed images · 10 arcseconds across · north up and east left", fontsize=12)
    save_figure(fig, FIGURES/"field_b_transient_ancillary")


def main():
    overlap_path = OUT/"field_b_transient_grism_overlap.json"
    overlap = json.loads(overlap_path.read_text())
    source = overlap["position_degrees"]
    sky = SkyCoord(*source, unit="deg")
    (g,) = [r for r in overlap["overlap_results"][0]["overlapping_sources"] if r["catalogue_label"] == 1178]
    galaxy = SkyCoord(*g["icrs"], unit="deg")
    products = records(DownloadCache(ROOT, STORAGE_MIB*1024**2), sky)
    images, tables = unpack_chandra(products["chandra_bundle"], sky)
    xray = xray_sources(tables, sky, galaxy)
    closest = xray[0]
    if (closest["table"], closest["id"]) != ("tbl04.fits", 726):
        raise ValueError("Unexpected closest source in the fixed CDF-S release")
    xray_pos = SkyCoord(*closest["icrs"], unit="deg")
    apertures = check_xray_pixels(images, {"transient": sky, "chandra_726": xray_pos})
    footprints = Table.read(verified_path(products["chandra_footprints"]), format="votable")
    visits = sorted({(str(r["obs_date"]), int(r["ObsId"]), str(r["Instrument"])) for r in footprints
                     if r["Instrument"] in ("ACIS-I", "ACIS-S")})
    if not visits or len(visits) != len({r[1] for r in visits}):
        raise ValueError("Expected unambiguous ACIS observation dates in the footprint response")
    radio = Table.read(verified_path(products["radio_query"]), format="votable")
    (rrow,) = radio
    radio_pos = SkyCoord(rrow["ra"], rrow["dec"], unit="deg")
    radio_match = dict(name=str(rrow["name"]), icrs=[radio_pos.ra.deg, radio_pos.dec.deg],
                       separation_from_transient_arcsec=float(sky.separation(radio_pos).arcsec),
                       separation_from_grism_galaxy_arcsec=float(galaxy.separation(radio_pos).arcsec),
                       catalogue_snr=float(rrow["snr_1p4_ghz"]),
                       integrated_flux=float(rrow["int_flux_1p4_ghz"]),
                       integrated_flux_error=float(rrow["int_flux_1p4_ghz_error"]),
                       flux_unit=str(radio["int_flux_1p4_ghz"].unit))
    optical, stamps = optical_stamps(products["optical_cutout"], sky)
    muse = {}
    for receipt in products["muse_maps"]:
        with fits.open(verified_path(receipt)) as h:
            x, y = np.rint(WCS(h[0].header).celestial.world_to_pixel(sky)).astype(int)
            muse[receipt["path"].rsplit("/", 1)[-1]] = float(h[0].data[y, x])
    jwst_path = OUT/"field_b_transient_grism_products.json"
    jwst = json.loads(jwst_path.read_text())
    obs = next(row for row in jwst["direct_images"] if "o006" in row["path"])
    (image, x, y, _), direct_metadata = direct_cutout(obs, source)
    # Use the north-up HST grid for the displayed JWST panel. Photometry uses native data.
    hst_image, hst_wcs, hst_xy = stamps["F814W"]
    hy, hx = np.indices(hst_image.shape, dtype=float)
    with fits.open(verified_path(obs), memmap=True) as h:
        wcs = WCS(h["SCI"].header).celestial
        origin = np.rint(np.array(wcs.world_to_pixel(sky))-[x, y]).astype(int)
        jx, jy = wcs.world_to_pixel_values(*hst_wcs.pixel_to_world_values(hx, hy))
    reprojected = map_coordinates(image, [jy-origin[1], jx-origin[0]], order=1, mode="constant", cval=np.nan)
    direct = reprojected, hst_wcs, hst_xy
    plot(stamps["F814W"], direct, images[SOFT_IMAGE], [
        ("Galaxy 1178", galaxy, "#ffca76", "+"),
        ("Chandra 726", xray_pos, "#ff8f9a", "x"),
        ("VLA 1.4 GHz", radio_pos, "#a8cff5", "1"),
    ], closest["positional_error_radius_68_arcsec"], direct_metadata["date_utc"][:10])
    write_json(OUT/"field_b_transient_ancillary.json", dict(
        executed_utc=utcnow(), candidate=overlap["candidate"], position_degrees=source,
        products=products, chandra_catalogue_rows={"main": 1008, "supplementary": 47},
        nearest_xray_sources=xray, soft_image_apertures=apertures,
        acis_footprints=dict(response_rows=len(footprints), distinct_observations=len(visits),
                             first_obs=visits[0], last_obs=visits[-1],
                             since_transient_bright_epoch=sum(row[0][:10] >= BRIGHT for row in visits)),
        radio_radius_arcsec=RADIO_RADIUS, radio_matches=[radio_match], hst_apertures=optical,
        muse_wide_dr1_position_samples=muse, jwst_direct_image=obs, jwst_direct_dateobs=direct_metadata["date_utc"],
        scope="Native Chandra mosaic apertures without PSF deblending; archival HST mosaics are static context",
        overlap_sha256=text_digest(overlap_path), jwst_products_sha256=text_digest(jwst_path),
        script_sha256=text_digest(__file__),
    ))


if __name__ == "__main__":
    main()
