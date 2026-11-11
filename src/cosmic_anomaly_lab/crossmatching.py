"""Spherical candidate matching; ambiguous associations remain unresolved."""

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky


def candidates(left, right, radius_arcsec=0.3):
    if not np.isfinite(radius_arcsec) or radius_arcsec <= 0:
        raise ValueError("Matching radius must be positive")
    if len(left) == 0 or len(right) == 0:
        return pd.DataFrame(columns=["left_index", "right_index", "separation_arcsec", "ambiguous"])
    a = SkyCoord(np.asarray(left.ra_deg) * u.deg, np.asarray(left.dec_deg) * u.deg)
    b = SkyCoord(np.asarray(right.ra_deg) * u.deg, np.asarray(right.dec_deg) * u.deg)
    i, j, sep, _ = search_around_sky(a, b, radius_arcsec * u.arcsec)
    result = pd.DataFrame({"left_index": i, "right_index": j, "separation_arcsec": sep.arcsec})
    result["ambiguous"] = result.left_index.duplicated(False) | result.right_index.duplicated(False)
    return result.sort_values(["left_index", "separation_arcsec"]).reset_index(drop=True)
