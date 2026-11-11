"""Version-specific DJA, JADES, and UNCOVER adapters; preserve missing values."""

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

from .provenance import write_json
from .registry import Registry, source_id

DJA_BASE = "https://s3.amazonaws.com/msaexp-nirspec/extractions"
COLUMNS = [
    "file",
    "srcid",
    "ra",
    "dec",
    "grating",
    "filter",
    "dataset",
    "root",
    "version",
    "exptime",
    "wmin",
    "wmax",
    "sn50",
    "flux50",
    "err50",
    "objid",
    "z_best",
    "zline",
    "grade",
    "valid",
    "comment",
    "z_phot",
    "phot_dr",
    "file_phot",
    "id_phot",
    "line_ha",
    "line_ha_err",
    "line_ha_nii",
    "line_ha_nii_err",
    "line_hb",
    "line_hb_err",
    "line_oiii",
    "line_oiii_err",
    "line_oiii_5007",
    "line_oiii_5007_err",
    "line_oiii_4959",
    "line_oiii_4959_err",
    "line_oii",
    "line_oii_err",
    "eqw_ha_nii",
    "eqw_oiii",
]


def survey_field(root):
    if root.startswith("ceers"):
        return "CEERS", "EGS"
    if root.startswith("uncover"):
        return "UNCOVER", "Abell2744"
    if root.startswith("rubies"):
        return "RUBIES", "UDS" if "uds" in root else "EGS"
    if root.startswith(("jades", "goodsn", "gds-deep")):
        return "JADES", "GOODS-N" if "gdn" in root or "goodsn" in root else "GOODS-S"
    return None, None


def clean(value):
    if pd.isna(value):
        return None
    return value.item() if isinstance(value, np.generic) else value


def ingest(root):
    root = Path(root)
    frame = pd.read_csv(root / "data/raw/dja_v4.4.csv.gz", usecols=COLUMNS, low_memory=False)
    total_rows = len(frame)
    mapped = frame.root.map(survey_field)
    frame["survey"] = [a for a, b in mapped]
    frame["field"] = [b for a, b in mapped]
    frame = frame.loc[frame.survey.notna()].copy()
    frame["programme"] = frame.dataset.str.extract(r"jw(\d{5})", expand=False).astype(int)
    frame["source_id"] = frame.objid.map(lambda x: source_id("DJA-v4.4-objid", str(x)))
    frame["url"] = DJA_BASE + "/" + frame.root + "/" + frame.file
    frame["redshift_method"] = "archive_spectroscopic_fit_visual_grade"
    frame.loc[frame.z_best < 0, "z_best"] = np.nan
    frame.loc[frame.z_phot < 0, "z_phot"] = np.nan
    registry = Registry(root / "data/processed/registry.sqlite")
    # Use a robust graded observation preferentially as the source representative.
    ordered = frame.sort_values(
        ["grade", "sn50", "file"], ascending=[False, False, True], na_position="last"
    )
    source_rows = []
    for objid, group in ordered.groupby("objid", sort=True):
        row = group.iloc[0]
        sid = registry.add_source(
            "DJA-v4.4-objid",
            str(objid),
            row.ra,
            row.dec,
            astrometric_uncertainty_arcsec=None,
            redshift=clean(row.z_best),
            redshift_uncertainty=None,
            redshift_method="archive_spectroscopic_fit_visual_grade",
            redshift_quality=clean(row.grade),
            physical_class="UNKNOWN_OR_UNRESOLVED",
            status="INSUFFICIENT_EVIDENCE",
            lensing_assessment="cluster_field_not_source_confirmation"
            if row.field == "Abell2744"
            else "not_assessed",
            namespace_note="DJA upstream association; not independently deblended or source-plane deduplicated",
            aliases=sorted(set(group.file)),
            surveys=sorted(set(group.survey)),
            fields=sorted(set(group.field)),
            programmes=sorted(set(int(x) for x in group.programme)),
        )
        source_rows.append(
            dict(
                source_id=sid,
                dja_objid=objid,
                ra_deg=row.ra,
                dec_deg=row.dec,
                survey=";".join(sorted(set(group.survey))),
                field=row.field,
                z_spec=clean(row.z_best),
                z_quality=clean(row.grade),
                n_observations=len(group),
                assessment="INSUFFICIENT_EVIDENCE",
                physical_class="UNKNOWN_OR_UNRESOLVED",
            )
        )
    for row in frame.itertuples(index=False):
        registry.add_observation(
            row.file,
            row.source_id,
            row.survey,
            row.field,
            row.url,
            programme=row.programme,
            instrument="NIRSpec",
            grating=row.grating,
            filter=row.filter,
            reduction="DJA-msaexp-v4.4",
            product_level="high_level_science_product",
            redshift=clean(row.z_best),
            line_fit_redshift=clean(row.zline),
            quality=clean(row.grade),
            exptime_s=clean(row.exptime),
            spectral_min_um=clean(row.wmin),
            spectral_max_um=clean(row.wmax),
            photometric_match_separation_arcsec=clean(row.phot_dr),
            upstream_comment=clean(row.comment),
        )
    registry.close()
    sources = pd.DataFrame(source_rows)
    sources.to_csv(root / "data/catalogues/sources.csv.gz", index=False)
    frame.to_csv(root / "data/catalogues/observations.csv.gz", index=False)
    jades = Table.read(root / "data/raw/jades_dr4.fits", hdu="Obs_info").to_pandas()
    for column in jades.select_dtypes(include=["object"]).columns:
        jades[column] = jades[column].map(
            lambda x: x.decode().strip() if isinstance(x, bytes) else x
        )
    jades = jades[
        [
            "Unique_ID",
            "PID",
            "RA_TARG",
            "Dec_TARG",
            "Field",
            "z_phot",
            "z_Spec",
            "z_Spec_flag",
            "oddball_source",
        ]
    ]
    jades = jades.rename(columns={"RA_TARG": "ra_deg", "Dec_TARG": "dec_deg"})
    jades.loc[jades.z_Spec < 0, "z_Spec"] = np.nan
    jades.loc[jades.z_phot < 0, "z_phot"] = np.nan
    jades.to_csv(root / "data/catalogues/jades_dr4_reference.csv.gz", index=False)
    uncover = Table.read(
        root / "data/raw/uncover_dr4_1.dat", format="ascii.commented_header"
    ).to_pandas()
    uncover = uncover.rename(columns={"ra": "ra_deg", "dec": "dec_deg"})
    for col in ["z_spec", "z_spec16", "z_spec50", "z_spec84"]:
        uncover.loc[uncover[col] < 0, col] = np.nan
    uncover.to_csv(root / "data/catalogues/uncover_dr4_1_reference.csv", index=False)
    summary = dict(
        dja_snapshot_rows=total_rows,
        ingested_observations=len(frame),
        ingested_source_records=len(sources),
        independent_physical_source_count=None,
        jades_reference_rows=len(jades),
        uncover_reference_rows=len(uncover),
        by_survey=frame.groupby("survey")
        .agg(observations=("file", "size"), source_records=("source_id", "nunique"))
        .reset_index()
        .to_dict("records"),
        fields=sorted(frame.field.unique().tolist()),
        note="Reference catalogue rows overlap DJA; never add these counts together.",
    )
    write_json(root / "results/catalogue_summary.json", summary)
    return frame, sources


def select_pilot(frame, seed=17775, config=None):
    """Deterministic PRISM selection, unrelated to labels or infrared colour cuts."""
    config = config or {}
    seed = config.get("selection_seed", seed)
    spectral = config.get("spectroscopy", {})
    grating = spectral.get("grating", "PRISM")
    grade = spectral.get("quality_grade", 3)
    min_sn = spectral.get("min_archive_median_pixel_sn", 5)
    bins = spectral.get("redshift_bins", [[0.5, 2], [2, 4], [4, 7], [7, 12]])
    counts = spectral.get("per_survey_per_bin", [2, 1, 1, 1])
    if len(bins) != len(counts) or any(n < 0 for n in counts):
        raise ValueError("Invalid bin/count configuration")
    if any(lo >= hi for lo, hi in bins):
        raise ValueError("Redshift bin lower bounds must precede upper bounds")
    eligible = frame.loc[
        (frame.grating == grating) & (frame.grade == grade) & (frame.sn50 >= min_sn)
    ].copy()
    eligible["order"] = eligible.file.map(
        lambda f: hashlib.sha256(f"{seed}:{f}".encode()).hexdigest()
    )
    eligible = eligible.sort_values("order").drop_duplicates("source_id")
    selections = []
    for survey in sorted(eligible.survey.unique()):
        subset = eligible[eligible.survey == survey]
        for (low, high), count in zip(bins, counts):
            candidates = subset[(subset.z_best >= low) & (subset.z_best < high)]
            for _, row in candidates.head(count).iterrows():
                record = row.to_dict()
                record.update(
                    sample="B" if low == 0.5 else "C",
                    selection=f"grade={grade}; {grating}; sn50>={min_sn}; {low}<=z<{high}; sha256 order seed={seed}",
                    population_label="low_z_comparison_candidate"
                    if low == 0.5
                    else "archive_selected_unclassified",
                    classification_reference=None,
                )
                selections.append(record)
    benchmark = frame[
        (frame.root == "uncover-v4") & (frame.srcid == 45924) & (frame.grating == "PRISM")
    ]
    for _, row in benchmark.head(1).iterrows():
        record = row.to_dict()
        record.update(
            sample="A",
            selection="Published LRD/AGN benchmark; excluded from blind selection",
            population_label="published_LRD_AGN_candidate",
            classification_reference="https://arxiv.org/abs/2412.04557",
        )
        selections.append(record)
    if not selections:
        return frame.iloc[:0].copy()
    return (
        pd.DataFrame(selections)
        .drop_duplicates("source_id", keep="last")
        .drop(columns="order", errors="ignore")
    )
