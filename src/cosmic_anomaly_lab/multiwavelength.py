"""Published photometry conversion without discarding negative measurements."""

import numpy as np


def ab25_flux(flux, error, weight):
    """Convert 3D-HST ZP=25 total-flux units to microJy; invalid inputs stay null.

    -99 is conservatively reserved as a missing-data sentinel. Other signed
    fluxes are retained. A catalogue use_phot flag must be checked separately.
    """
    if not np.isfinite([flux, error, weight]).all() or error <= 0 or weight <= 0 or flux == -99:
        return None, None
    return float(flux * .3631), float(error * .3631)
