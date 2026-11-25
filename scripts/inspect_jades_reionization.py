"""Cross-check Lyα screening responses against published data and observed images."""

import json
import re
from html.parser import HTMLParser
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import urlencode, urljoin

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from download_jades_reionization_sample import LEDGER, OUT, ROOT
from inspect_jades_cosmic_evolution import DENSITY_UNIT, receipt
from matplotlib.figure import Figure
from measure_jades_reionization import C_KMS, LYA, WIDTHS, line_models, spectral_arrays

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.visual_evidence import image_panel, read_images, save_figure

PUBLISHED_TABLE = "https://arxiv.org/html/2409.06405"
PUBLISHED_ID = re.compile(r"G([NS])\+(\d+\.\d+)([+−-])(\d+\.\d+)$")
PUBLISHED_FLUX = re.compile(r"(\d+(?:\.\d+)?)\s*±\s*(\d+(?:\.\d+)?)")


class LiteratureTable(HTMLParser):
    """Read the coordinate labels and the R1000 flux column in the author table."""

    def __init__(self):
        super().__init__()
        self.cells = None
        self.value = None
        self.in_annotation = False
        self.rows = []

    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.cells = []
        elif tag == "td" and self.cells is not None:
            self.value = ""
        elif tag == "annotation":
            self.in_annotation = True

    def handle_data(self, data):
        if self.value is not None and not self.in_annotation:
            self.value += " " + data

    def handle_endtag(self, tag):
        if tag == "annotation":
            self.in_annotation = False
        if tag == "td" and self.value is not None:
            self.cells.append(" ".join(self.value.split()))
            self.value = None
        if tag == "tr" and self.cells is not None:
            if len(self.cells) == 9 and (name := PUBLISHED_ID.fullmatch(self.cells[0])):
                flux = PUBLISHED_FLUX.match(self.cells[7])
                self.rows.append(dict(name=name[0].replace("−","-"), field=f"G{name[1]}", ra=float(name[2]),
                                      dec=float(name[4]) * (1 if name[3]=="+" else -1),
                                      z=float(self.cells[1]), grating_line_reported=flux is not None,
                                      literature_flux=float(flux[1]) if flux else None,
                                      literature_error=float(flux[2]) if flux else None))
            self.cells = None


def literature_matches(cache, rows, observations):
    raw = cache.fetch(PUBLISHED_TABLE, "metadata/jones_2409_06405.html", survey="JADES Jones et al.",
                      version="arXiv HTML", kind="published_literature_table")
    parser = LiteratureTable()
    parser.feed(raw.read_text(encoding="utf8"))
    if len(parser.rows) != 150:
        raise ValueError(f"Expected 150 published catalogue rows, found {len(parser.rows)}")
    sky = SkyCoord([r["ra"] for r in parser.rows], [r["dec"] for r in parser.rows], unit="deg")
    result = []
    for row in rows:
        delta = SkyCoord(row["ra"], row["dec"], unit="deg").separation(sky).arcsec
        matches = [dict(p, separation_arcsec=float(sep)) for sep, p in zip(delta, parser.rows, strict=True)
                   if sep <= .20 and row["field"] == p["field"] and abs(row["z"]-p["z"])/(1+p["z"]) <= .002]
        fitted = observations.loc[row["dr4_id"]]
        result.append(dict(dr4_id=row["dr4_id"], searched=True, above_screen=bool(fitted.above_screen),
                           matches=matches))
    return receipt(raw, PUBLISHED_TABLE), result


def inspect_images(cache, entry, row, image_axis, slit_axis):
    cutout = "https://grizli-cutout.herokuapp.com/thumb?" + urlencode(dict(
        ra=entry["ra"], dec=entry["dec"], size=3, filters="f150w-clear,f277w-clear,f444w-clear", output="fits_weight"
    ))
    image_path = cache.fetch(cutout, f"images/jades_lya_{entry['dr4_id']}.fits", survey="DAWN JWST Archive",
                             version="DJA-cutout-live", kind="observed_image_cutout")
    images, rgb, details = read_images(image_path, entry["ra"], entry["dec"])
    band = next(iter(images), None) if rgb is None else None
    image_panel(image_axis, images, rgb, SimpleNamespace(**entry), band=band)
    image_axis.set_title("JWST NIRCam RGB" if rgb is not None else band or "Unavailable", fontsize=9)
    listing = ROOT / entry["listing"]["path"]
    if digest(listing) != entry["listing"]["sha256"]:
        raise ValueError("Directory receipt mismatch")
    name = Path(entry["spectrum"]["path"]).name.replace("_x1d.fits", "_s2d.fits")
    if name not in re.findall(r'href="([^"]+)"', listing.read_text()):
        raise ValueError(f"2D counterpart missing from official directory for {entry['dr4_id']}")
    url = urljoin(entry["listing"]["url"], name)
    path = cache.fetch(url, f"independent/{name}", survey="JADES", version="DR4", kind="observed_spectral_image")
    yvel = LYA*(1+row.z)
    with fits.open(path) as h:
        header = h[0].header
        sep = SkyCoord(header["RA"], header["DEC"], unit="deg").separation(
            SkyCoord(entry["ra"], entry["dec"], unit="deg")
        ).arcsec
        if sep > .2 or int(header["PROPOSID"]) != entry["header_pid"] or int(header["TARGET"]) != entry["header_target_number"]:
            raise ValueError(f"2D product failed header checks for {entry['dr4_id']}")
        w = h["WAVELENGTH"]
        velocity = ((w.data*u.Unit(w.header["BUNIT"])).to_value(u.AA)/yvel-1)*C_KMS
        window = (velocity >= -1800) & (velocity <= 3300)
        plane = h["FLUX"]
        scale = u.Unit(plane.header["BUNIT"]).to(DENSITY_UNIT)/1e-20
        err = h["FLUX_ERR"].data[:,window]*scale
        z = plane.data[:,window]*scale
        v = np.where(np.isfinite(err) & (err > 0), z/np.where(err > 0, err, np.nan), np.nan)
        spatial = plane.header["CRVAL2"] + (np.arange(len(z))+1-plane.header["CRPIX2"])*plane.header["CDELT2"]
        slit_axis.pcolormesh(velocity[window], spatial, v, shading="nearest", cmap="RdBu_r", vmin=-4, vmax=4)
        slit_axis.set_ylim(min(spatial),max(spatial))
    slit_axis.axhline(0,color=".3",linestyle=":",linewidth=.7)
    slit_axis.axvline(row.best_velocity_kms,color=".3",linestyle=":",linewidth=.7)
    slit_axis.set(title=f"{entry['dr4_id']} · z={row.z:.3f}", ylabel="Slit position (shutters)", xlim=(-1800,3300))
    tests = []
    for a,n in ((5,3),(3,3),(5,2)):
        s = read_dr4_spectrum(ROOT/entry["spectrum"]["path"], aperture=a, nods=n, wavelength_range=(.65,5.4))
        models = line_models(spectral_arrays(s),yvel,WIDTHS[0])
        fixed = next((p for p in models if p["offset_kms"] == row.best_velocity_kms), None)
        tests.append(dict(aperture=a,nods=n,snr=fixed["flux"]/fixed["error"] if fixed else None))
    return dict(dr4_id=entry["dr4_id"], image=receipt(image_path,cutout), wcs_diagnostics=details,
                spectral_image=receipt(path,url), header_offset_arcsec=float(sep), extraction_checks=tests)


def atlas(cache, ids, targets, observed, stem, title):
    fig = Figure(figsize=(10., 2.4*len(ids)+.6), layout="constrained")
    axes = fig.subplots(len(ids),2,squeeze=False,gridspec_kw=dict(width_ratios=(1,3)))
    images = []
    for sid,(left,right) in zip(ids,axes,strict=True):
        images.append(inspect_images(cache,targets[sid],observed.loc[sid],left,right))
        print(sid,images[-1]["extraction_checks"],flush=True)
    axes[-1,1].set_xlabel("Velocity relative to systemic Lyα (km/s)")
    fig.suptitle(f"{title}\n2D colours show signed flux divided by released pixel error (±4)",fontsize=13)
    save_figure(fig,ROOT/f"figures/{stem}")
    return images


def main():
    sample = json.loads(LEDGER.read_text())
    summary = OUT/"summary.json"
    data = json.loads(summary.read_text())
    observed = pd.read_csv(OUT/"measurements.csv").set_index("dr4_id")
    targets = {r["dr4_id"]:r for r in sample["records"]}
    cache = DownloadCache(ROOT, 28_000*1024**2)
    reference, matches = literature_matches(cache,sample["records"],observed)
    images = atlas(cache, data["example_ids"], targets, observed, "jades_reionization_examples",
                   "Observed JADES galaxies and 2D G140M spectra")
    unseen = [row["dr4_id"] for row in matches if row["above_screen"] and not row["matches"]]
    unmatched = atlas(cache, unseen, targets, observed, "jades_reionization_unmatched",
                      "Lyα responses without a match in the Jones et al. table") if unseen else []
    write_json(OUT/"crosschecks.json",dict(
        executed_utc=utcnow(), summary_sha256=text_digest(summary), script_sha256=text_digest(__file__),
        published_table=reference, literature_matches=matches, examples=images, catalogue_unmatched=unmatched,
        literature_scope="0.20 arcsec and normalized redshift agreement within 0.002 in Jones et al. (2025).",
    ))


if __name__ == "__main__":
    main()
