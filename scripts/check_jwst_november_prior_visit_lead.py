"""Check a November GOODS-S lead in two additional JWST epochs."""

import json

from astropy.coordinates import SkyCoord
from download_jwst_cross_detector import OUT, ROOT
from matplotlib.figure import Figure
from search_jwst_november_prior_visits import STEM
from vet_jwst_field_b_candidate import FIGURES
from vet_jwst_november_catalogues import STORAGE_MIB, inspect

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

CANDIDATE = "GSN-53.184076-27.806633"
PAIRS = (("08060-001", "12577-001"), ("12577-001", "12577-004"))
OUTPUT = f"{STEM}_october_check"


def plot_flux_changes(pairs):
    fig = Figure(figsize=(8.2, 3.8), layout="constrained")
    ax = fig.subplots()
    filters = ("F150W", "F200W", "F277W", "F356W")
    for pair, offset in zip(pairs, (-.19, 0., .19), strict=True):
        samples = []
        for i, band in enumerate(filters):
            profile = pair["bands"][band].get("change")
            if profile and None not in (profile[-1]["flux_njy"], profile[-1]["local_mad_sigma_njy"]):
                change = profile[-1]
                samples.append((i+offset, change["flux_njy"], change["local_mad_sigma_njy"]))
        if samples:
            x, y, err = zip(*samples, strict=True)
            ax.errorbar(x, y, yerr=err, fmt="o", capsize=3,
                        label=f"{pair['older']} → {pair['newer']}")
    ax.axhline(0., color="0.45", linewidth=1, zorder=0)
    ax.set_xticks(range(len(filters)), filters)
    ax.set(xlabel="JWST NIRCam filter", ylabel="Flux change in 0.24″ aperture (nJy)", title=CANDIDATE)
    ax.legend(frameon=False, fontsize=8, ncols=3, loc="best")
    save_figure(fig, FIGURES/f"{OUTPUT}_flux_changes")


def main():
    vetting = OUT/f"{STEM}_vetting.json"
    previous = next(row for row in json.loads(vetting.read_text())["candidates"] if row["candidate"] == CANDIDATE)
    sky = SkyCoord(*previous["icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    comparisons = []
    for old, new in PAIRS:
        lead = {key: previous[key] for key in ("candidate", "icrs_deg", "catalogue_bands", "sign")}
        lead.update(older=old, newer=new, image_basename=f"{OUTPUT}_{old}_{new}")
        lead["bands"] = inspect(cache, sky, lead)
        comparisons.append(lead)
    plot_flux_changes((previous, *comparisons))
    write_json(OUT/f"{OUTPUT}.json", dict(
        executed_utc=utcnow(), candidate=CANDIDATE, comparisons=comparisons,
        vetting_sha256=text_digest(vetting), script_sha256=text_digest(__file__),
        measurement_method_sha256=text_digest(ROOT/"scripts/vet_jwst_november_catalogues.py"),
        scope="Same-position forced photometry in previously unmeasured visits covering the candidate.",
    ))


if __name__ == "__main__":
    main()
