import sys

sys.path.insert(0, "src")
from cosmic_anomaly_lab.archives import product_inventory, search_mast
from cosmic_anomaly_lab.provenance import utcnow, write_json

results = []
for survey, pid in [("JADES", 1210), ("CEERS", 1345), ("UNCOVER", 2561)]:
    try:
        obs = search_mast(pid)
        cols = [
            x
            for x in [
                "obsid",
                "obs_id",
                "proposal_id",
                "instrument_name",
                "filters",
                "calib_level",
                "dataRights",
                "s_ra",
                "s_dec",
                "t_exptime",
            ]
            if x in obs.colnames
        ]
        obs[cols].write(f"data/catalogues/mast_{pid}.ecsv", format="ascii.ecsv", overwrite=True)
        products = product_inventory(obs, 1)
        if products is not None:
            products.write(
                f"data/catalogues/mast_{pid}_products.ecsv", format="ascii.ecsv", overwrite=True
            )
        result = dict(
            survey=survey,
            pid=pid,
            status="queried",
            observations=len(obs),
            products_first_observation=len(products) if products is not None else 0,
        )
    except Exception as exc:
        result = dict(survey=survey, pid=pid, status="failed", error=str(exc))
    print(result, flush=True)
    results.append(result | {"accessed_utc": utcnow()})
    write_json("data/manifests/mast_queries.json", results)
