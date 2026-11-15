"""Conditional local Halpha models; not a physical-class or mass estimator."""

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares

from .spectroscopy import pixel_edges

C_KMS = 299792.458
PARAMETERS = [
    "continuum",
    "slope",
    "narrow_area",
    "narrow_sigma",
    "narrow_velocity",
    "nii6585_area",
    "broad1_area",
    "broad1_sigma",
    "broad2_area",
    "broad2_sigma",
    "broad_velocity",
]


def gaussian(v, area, sigma, center):
    """Integral over velocity is 1000 * area; no luminosity conversion."""
    return area * 1000 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * ((v - center) / sigma) ** 2)


class BalmerModel:
    """Apply broad-only absorption BEFORE Gaussian LSF and pixel integration.

    Fixed covering fraction one: tau is an effective local shape parameter.
    Narrow Halpha and [NII], and a linear continuum, remain unabsorbed.
    Gaussian LSF is constant in this local velocity coordinate, not calibrated
    to this source's slit illumination. The [NII] 6585/6549 area ratio is 2.96.
    """

    def __init__(self, velocity, resolution, family="absorbed", step=5):
        self.velocity = np.asarray(velocity, float)
        if family not in ["co_centered", "absorbed", "split_emission"]:
            raise ValueError("Unknown family")
        if not np.isfinite(resolution) or resolution <= 0 or step <= 0:
            raise ValueError("Positive resolution and step required")
        self.family = family
        self.edges = pixel_edges(self.velocity)
        self.sigma_lsf = C_KMS / (2.354820045 * resolution)
        self.step = step
        margin = max(8 * self.sigma_lsf, 100)
        self.grid = np.arange(self.edges[0] - margin, self.edges[-1] + margin + step, step)

    @property
    def names(self):
        return (
            PARAMETERS
            + (
                {
                    "absorbed": ["tau", "abs_sigma", "abs_velocity"],
                    "split_emission": ["broad2_velocity"],
                    "co_centered": [],
                }[self.family]
            )
        )

    def predict(self, p):
        v = self.grid
        broad = gaussian(v, p[6], p[7], p[10]) + gaussian(
            v, p[8], p[9], p[11] if self.family == "split_emission" else p[10]
        )
        if self.family == "absorbed":
            broad *= np.exp(-p[11] * np.exp(-0.5 * ((v - p[13]) / p[12]) ** 2))
        line = broad + gaussian(v, p[2], p[3], p[4])
        for rest, ratio in [(6585.27, 1), (6549.86, 1 / 2.96)]:
            line += gaussian(v, p[5] * ratio, p[3], C_KMS * (rest / 6564.61 - 1) + p[4])
        blurred = gaussian_filter1d(line, self.sigma_lsf / self.step, mode="constant", truncate=8)
        integral = cumulative_trapezoid(blurred, v, initial=0)
        binned = np.diff(np.interp(self.edges, v, integral)) / np.diff(self.edges)
        # Exact linear-continuum average on each (possibly nonuniform) pixel.
        return binned + p[0] + p[1] * (self.edges[:-1] + self.edges[1:]) / 10000


def fit_balmer(velocity, flux, error, good, resolution, family, *, starts=4, step=5):
    velocity, flux, error, good = map(np.asarray, (velocity, flux, error, good))
    if any(x.shape != velocity.shape for x in [flux, error, good]):
        raise ValueError("Matching shapes required")
    mask = good & np.isfinite(flux + error) & (error > 0)
    model = BalmerModel(velocity, resolution, family, step)
    if mask.sum() <= len(model.names) + 5:
        raise ValueError("Insufficient valid pixels")
    scale = max(float(np.percentile(abs(flux[mask]), 80)), float(np.median(error[mask])))
    y, e = flux / scale, error / scale
    lower = [-5, -10, 0, 20, -400, 0, 0, 200, 0, 200, -1200]
    upper = [10, 10, 100, 300, 400, 100, 100, 3000, 100, 4000, 1200]
    if family == "absorbed":
        lower += [0, 25, -500]
        upper += [20, 500, 500]
    elif family == "split_emission":
        lower += [-1200]
        upper += [1200]
    attempts = []
    solutions = []
    for index in range(starts):
        initial = [
            float(np.clip(np.percentile(y[mask], 15), -4, 9)),
            0,
            0.5,
            65,
            0,
            0.02,
            2,
            [450, 700, 1000, 300][index % 4],
            2,
            1500,
            0,
        ]
        if family == "absorbed":
            initial += [[1, 3, 0.2, 6][index % 4], [100, 150, 70, 200][index % 4], 0]
        elif family == "split_emission":
            initial[10] = [-300, 300, -700, 700][index % 4]
            initial += [-initial[10]]
        result = least_squares(
            lambda p: ((model.predict(p) - y) / e)[mask],
            initial,
            bounds=(lower, upper),
            x_scale="jac",
            max_nfev=1200,
            ftol=1e-8,
            xtol=1e-8,
            gtol=1e-8,
        )
        chi2 = float(np.sum(result.fun**2))
        attempts.append(dict(success=bool(result.success), nfev=result.nfev, chi2=chi2))
        if result.success and np.isfinite(chi2):
            solutions.append(result)
    if not solutions:
        raise ValueError("No converged fit")
    best = min(solutions, key=lambda s: np.sum(s.fun**2))
    p = best.x.copy()
    for i in [0, 1, 2, 5, 6, 8]:
        p[i] *= scale
    prediction = model.predict(p)
    chi2 = float(np.sum(((flux[mask] - prediction[mask]) / error[mask]) ** 2))
    n, k = int(mask.sum()), len(p)
    near = (best.x - lower < 1e-4 * (np.array(upper) - lower)) | (
        upper - best.x < 1e-4 * (np.array(upper) - lower)
    )
    return dict(
        family=family,
        resolution=resolution,
        chi2=chi2,
        n=n,
        k=k,
        bic=chi2 + k * float(np.log(n)),
        reduced_chi2=chi2 / (n - k),
        parameters=dict(zip(model.names, p.tolist(), strict=True)),
        parameter_vector=p.tolist(),
        parameters_near_bounds=np.array(model.names)[near].tolist(),
        attempts=attempts,
        model=prediction.tolist(),
        fit_mask=mask.tolist(),
        flux_scale=scale,
        grid_step_kms=step,
    )
