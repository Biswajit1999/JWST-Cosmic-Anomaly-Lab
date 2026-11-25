"""Show real NIRCam images and remeasure DR4 spectra for four sample galaxies."""

import json
import re
from types import SimpleNamespace
from urllib.parse import urlencode, urljoin

import numpy as np
from analyse_jades_cosmic_evolution import FIELD_NAMES, OUT, ROOT
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
from matplotlib.figure import Figure

from cosmic_anomaly_lab.balmer_census import signed_core_integral
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.visual_evidence import image_panel, read_images, save_figure

LINES_VACUUM_A = {"Hβ": 4862.64, "[O III] 4959": 4960.29, "[O III] 5007": 5008.24}
DENSITY_UNIT = u.erg / u.s / u.cm**2 / u.AA
DISPERSER = "f170lp-g235m"


def receipt(path, url):
    return dict(path=path.relative_to(ROOT).as_posix(), url=url, sha256=digest(path))


def spectrum_for_target(cache, row):
    tier, number = row.dr4_id.rsplit("_", 1)
    directory = f"https://jades.herts.ac.uk/DR4/{tier[:7]}/spectra/{DISPERSER}/{tier}/"
    listing = cache.fetch(directory, f"metadata/dr4_{tier}_g235m_f170lp.html", survey="JADES",
                          version="DR4-listing", kind="archive_directory")
    links = re.findall(r'href="([^"]+)"', listing.read_text())
    matches = [link for link in links if link.endswith("_x1d.fits")
               and (m := re.search(r"-(\d+)_" + DISPERSER + "_", link))
               and int(m.group(1)) == int(number)]
    if len(matches) != 1:
        raise ValueError(f"{row.dr4_id} has {len(matches)} listed G235M 1D products")
    url = urljoin(directory, matches[0])
    path = cache.fetch(url, f"independent/{matches[0]}", survey="JADES", version="DR4", kind="observed_spectrum")
    spec = read_dr4_spectrum(path, aperture=5, nods=3)
    off = SkyCoord(spec.metadata["RA"], spec.metadata["DEC"], unit="deg").separation(
        SkyCoord(row.ra, row.dec, unit="deg")
    ).arcsec
    if off > .2 or int(spec.metadata["PROPOSID"]) != row.pid or int(spec.metadata["TARGET"]) != int(number):
        raise ValueError(f"{row.dr4_id} failed the coordinate, programme or ID check")
    return spec, dict(listing=receipt(listing, directory), spectrum=receipt(path, url), header_offset_arcsec=float(off),
                      header_pid=row.pid, header_target_number=int(number))


def plot_spectrum(ax, spec, row):
    rest = spec.wavelength.to_value(u.AA) / (1 + row.z)
    flux = spec.flux.to_value(DENSITY_UNIT, equivalencies=u.spectral_density(spec.wavelength)) / 1e-20
    err = spec.error.to_value(DENSITY_UNIT, equivalencies=u.spectral_density(spec.wavelength)) / 1e-20
    keep = (rest >= 4800) & (rest <= 5080)
    y = np.where(spec.good, flux, np.nan)
    ax.step(rest[keep], y[keep], where="mid", color="#236c7a", linewidth=1.0)
    ax.fill_between(rest[keep], (y-err)[keep], (y+err)[keep], step="mid", color="#236c7a", alpha=.22, linewidth=0)
    ax.axhline(0, color="0.4", linewidth=.6)
    ax.set(xlim=(4800, 5080), ylabel="Fλ (10⁻²⁰ cgs/Å)")
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom, top+.13*(top-bottom))
    lines = {}
    for name, wavelength in LINES_VACUUM_A.items():
        velocity = (rest / wavelength - 1) * 299792.458
        lines[name] = signed_core_integral(velocity, flux, err, spec.good, wavelength * (1+row.z), 600)
        ax.axvline(wavelength, linewidth=.6, linestyle=":", color="#a75134")
        ax.text(wavelength, .95, name, color="#8f422b", fontsize=8, ha="center", va="top",
                transform=ax.get_xaxis_transform())
    return lines


def main():
    summary = OUT / "summary.json"
    targets = [SimpleNamespace(**row) for row in json.loads(summary.read_text())["example_targets"]]
    cache = DownloadCache(ROOT, 28_000*1024**2)
    fig = Figure(figsize=(10.2, 10.6), layout="constrained")
    axes = fig.subplots(len(targets), 2, gridspec_kw={"width_ratios": (1, 2.9)})
    products = []
    for i, (row, (image_ax, ax)) in enumerate(zip(targets, axes, strict=True)):
        spec, record = spectrum_for_target(cache, row)
        url = "https://grizli-cutout.herokuapp.com/thumb?" + urlencode(dict(
            ra=row.ra, dec=row.dec, size=3, filters="f150w-clear,f277w-clear,f444w-clear", output="fits_weight"
        ))
        path = cache.fetch(url, f"images/jades_evolution_{row.dr4_id}.fits", survey="DAWN JWST Archive",
                           version="DJA-cutout-live", kind="observed_image_cutout")
        images, rgb, details = read_images(path, row.ra, row.dec)
        band = next(iter(images), None) if rgb is None else None
        if not images:
            raise ValueError(f"No observed image coverage for {row.dr4_id}")
        image_panel(image_ax, images, rgb, row, band=band)
        image_ax.set_title("JWST NIRCam RGB" if rgb is not None else band, fontsize=8)
        lines = plot_spectrum(ax, spec, row)
        age = Planck18.age(row.z).to_value(u.Gyr)
        ax.set_title(f"{FIELD_NAMES[row.field]} · DR5 {row.nircam_id} · z={row.z:.3f} · {age:.2f} Gyr",
                     fontsize=11, loc="left")
        if i == len(targets)-1:
            ax.set_xlabel("Rest wavelength (vacuum Å)")
        products.append(dict(dr4_id=row.dr4_id, **record, image=receipt(path, url), image_bands=list(images),
                             image_display="F444W/F277W/F150W RGB" if rgb is not None else band,
                             image_diagnostics=details, signed_600_kms_line_cores=lines))
        print(row.dr4_id, list(images), lines["Hβ"].get("status"), flush=True)
    fig.suptitle("Observed JADES galaxies near each field's median line ratio\n"
                 "Endpoint examples · NIRSpec G235M, public 5-pixel extraction", fontsize=13)
    save_figure(fig, ROOT / "figures/jades_cosmic_evolution_observed")
    write_json(OUT / "observed_examples.json", dict(
        executed_utc=utcnow(), input_sha256=text_digest(summary), script_sha256=text_digest(__file__),
        products=products, sample_scope="Illustrative detected sources selected by distance to their field/bin median.",
        measurement="Fixed ±600 km/s line cores after sideband continuum fit; catalogue and spectra share observations.",
    ))


if __name__ == "__main__":
    main()
