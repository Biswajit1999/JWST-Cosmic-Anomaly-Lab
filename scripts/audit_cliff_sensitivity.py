"""Test whether The Cliff's empirical tail preference survives masking its core."""

from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.io import fits

from cosmic_anomaly_lab.line_profiles import C_KMS, fit_complex
from cosmic_anomaly_lab.provenance import digest, utcnow, write_json
from cosmic_anomaly_lab.spectroscopy import read_spectrum

root = Path(__file__).resolve().parents[1]
path = root / "data/raw/spectra/rubies-uds31-v4_g395m-f290lp_4233_154183.spec.fits"
spec = read_spectrum(path)
z = 3.546018
wave = spec.wavelength.to_value(u.AA) / (1 + z)
window = (wave >= 6420) & (wave <= 6800)
conversion = (np.ones(len(wave)) * u.uJy).to_value(
    u.erg / u.s / u.cm**2 / u.AA, equivalencies=u.spectral_density(spec.wavelength)
) / 1e-20
with fits.open(root / "data/raw/instrument/g395m_disp.fits") as hdus:
    nominal = float(np.interp(0.656461 * (1 + z), hdus[1].data["WAVELENGTH"], hdus[1].data["R"]))
components = [
    [(6564.61, 1)],
    [(6549.86, 1 / 3), (6585.27, 1)],
    [(6679.995, 1)],
    [(6718.29, 1)],
    [(6732.67, 1)],
]
records = []
for core_mask in [False, True]:
    good = spec.good[window].copy()
    if core_mask:
        good &= abs(C_KMS * (wave[window] / 6564.61 - 1)) > 400
    for lsf in [1, 1.3, 1.7]:
        fits_by_family = [
            fit_complex(
                wave[window],
                (spec.flux.value * conversion)[window],
                (spec.error.value * conversion)[window],
                good,
                main_line=6564.61,
                components=components,
                resolution=nominal * lsf,
                family=family,
            )
            for family in ["gaussian", "exponential"]
        ]
        g, e = fits_by_family
        record = dict(
            mask_core_400_km_s=core_mask,
            lsf_scale=lsf,
            n_pixels=g["n"],
            delta_bic_gaussian_minus_exponential=g["bic"] - e["bic"],
            delta_bic_doubled_errors=(g["chi2"] - e["chi2"]) / 4
            + (g["k"] - e["k"]) * np.log(g["n"]),
            gaussian_reduced_chi2=g["reduced_chi2"],
            exponential_reduced_chi2=e["reduced_chi2"],
            parameters_at_bound_gaussian=g["parameter_indices_near_bound"],
            parameters_at_bound_exponential=e["parameter_indices_near_bound"],
        )
        records.append(record)
        print(record, flush=True)
write_json(
    root / "results/cliff_sensitivity.json",
    dict(
        executed_utc=utcnow(),
        raw_sha256=digest(path),
        redshift=z,
        reference="https://arxiv.org/abs/2503.16600",
        cases=records,
        interpretation="Post-screen published-benchmark diagnostic; no independent physical classification",
        limitations=[
            "No absorption model",
            "Fixed archive redshift with free local centroids",
            "Diagonal errors and Gaussian LSF",
            "Core masking can remove narrow-component information",
            "No multiple-testing discovery probability or independent reduction",
        ],
    ),
)
