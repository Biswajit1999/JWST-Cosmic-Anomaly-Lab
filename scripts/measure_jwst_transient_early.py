"""Check F150W images from the six days before the transient search visit."""

import json
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from download_jwst_cross_detector import API, OUT, ROOT, download
from download_jwst_field_b_candidate import STORAGE_MIB
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import format_flux, pixel_scale, show_stamp
from vet_jwst_field_b_candidate import FIGURES, vet_mosaic

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

VISITS = [
    ("jw01180013001", "NRCA2"),
    ("jw04549005001", "NRCA2"),
    ("jw04549006001", "NRCA2"),
    ("jw04549006001", "NRCA4"),
    ("jw08060002001", "NRCA2"),
    ("jw12577002001", "NRCA2"),
]


def main():
    lead_path = OUT / "field_b_transient_vetting.json"
    lead = json.loads(lead_path.read_text())
    sky = tuple(lead["position_degrees"])
    flt = "F150W-CLEAR"
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    url = API + "/assoc_mosaic?" + urlencode(dict(coords=f"{sky[0]:.8f},{sky[1]:.8f},.02", filters=flt, output="csv"))
    association_list = download(cache, url, "transient_f150w_associations.csv", "archival_multiepoch_photometry")
    entries = pd.read_csv(ROOT / association_list["path"])
    records, mosaics, rows = [], [], []
    for visit, detector in VISITS:
        found = entries.query("visit == @visit and detector == @detector and tile == 672 and status == 2")
        if len(found) != 1:
            raise ValueError(f"Expected one F150W product for {visit} {detector}")
        row = found.iloc[0]
        product = dict(band=flt, detector=detector, visit=visit, mjd=float(row.expstart), dateobs=row.dateobs,
                       exposure_count=int(row.nexp), assoc_name=row.assoc_name)
        for kind in ("sci", "wht"):
            url = row.file.replace("_sci.fits.gz", f"_{kind}.fits.gz")
            product[kind] = download(cache, url, url.rsplit("/", 1)[1], "observed_candidate_epoch_"+kind)
        records.append(product)
        mosaic, measurements = vet_mosaic(product, sky)
        mosaics.append(mosaic)
        rows.append(measurements)
        print(row.dateobs, detector, format_flux(measurements["apertures"][-1]), flush=True)
    display_samples = [im[round(y)-7:round(y)+8, round(x)-7:round(x)+8].ravel()
                       for im, _, x, y, _, _ in mosaics]
    limit = float(np.nanpercentile(np.concatenate(display_samples), 99))
    figure = Figure(figsize=(8.0, 5.7), layout="constrained")
    for ax, (im, wcs, x, y, _, _), row in zip(figure.subplots(2, 3).flat, mosaics, rows, strict=True):
        show_stamp(ax, im, x, y, pixel_scale(wcs), (-.08, limit), "gray", .24, half_arcsec=.60)
        date = row["dateobs"][:16]
        ax.set_title(f"{date} · {row['detector']}\n{format_flux(row['apertures'][-1])}", fontsize=8.5)
    figure.suptitle("Early F150W imaging of the GOODS-S transient\n"
                    "Observed cutouts · UTC dates · fixed 0.24″ aperture", fontsize=13)
    save_figure(figure, FIGURES / "field_b_transient_f150w")
    write_json(OUT / "field_b_transient_f150w.json", {
        "executed_utc": utcnow(), "candidate": lead["candidate"], "position_degrees": sky,
        "products": records, "association_list": association_list, "measurements": rows,
        "main_vetting_sha256": text_digest(lead_path), "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()
