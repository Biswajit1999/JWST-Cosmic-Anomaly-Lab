import numpy as np

from cosmic_anomaly_lab.multiwavelength import ab25_flux


def test_ab25_keeps_negative_measurement_and_rejects_bad_errors():
    assert ab25_flux(-2, 1, .5) == (-.7262, .3631)
    assert ab25_flux(1, -2, .5) == (None, None)
    assert ab25_flux(1, 2, 0) == (None, None)
    assert ab25_flux(-99, 1, 1) == (None, None)
    assert ab25_flux(np.nan, 1, 1) == (None, None)
