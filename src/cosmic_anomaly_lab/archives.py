"""Metadata-first public MAST searches, programme inventory and bounded product lists."""

from astroquery.mast import Observations


def search_mast(programme, product_type="spectrum"):
    if not str(programme).isdigit():
        raise ValueError("Programme ID must be numeric")
    Observations.TIMEOUT = 45
    return Observations.query_criteria(
        obs_collection="JWST",
        proposal_id=str(programme),
        dataproduct_type=product_type,
        dataRights="PUBLIC",
        intentType="science",
    )


def product_inventory(observations, max_observations=3):
    if max_observations < 1:
        raise ValueError("A positive observation budget is required")
    if len(observations) == 0:
        return None
    return Observations.get_product_list(observations[:max_observations])


def processing_level(calib_level):
    return {
        0: "raw",
        1: "raw_or_minimally_processed",
        2: "calibrated",
        3: "high_level_science_product",
        4: "high_level_science_product",
    }.get(int(calib_level), "unknown")
