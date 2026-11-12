"""Regenerate two A1 conference posters from verified result files.

Run with Python providing reportlab and pypdf (the Codex bundled runtime has both).
Plots are embedded as PDF pages to preserve vector text and curves. Actual image
pixels remain raster observations. No generated astronomical imagery is used.
"""

import hashlib
import json
from datetime import date
from pathlib import Path

from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.lib.colors import HexColor, white
from reportlab.lib.pagesizes import A1, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/pdf"
OUT.mkdir(parents=True, exist_ok=True)
TMP = ROOT / "data/processed/posters"
TMP.mkdir(parents=True, exist_ok=True)
W, H = landscape(A1)
INK = HexColor("#14273A")
TEAL = HexColor("#1D6872")
PALE = HexColor("#EAF1F2")
MUTED = HexColor("#536472")


def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


summary = read("results/expanded/summary.json")
cliff = read("results/cliff_sensitivity.json")
inspection = read("results/expanded/inspection_audit.json")
snapshot_date = date.fromisoformat(summary["executed_utc"][:10])
snapshot_label = snapshot_date.strftime("%d %B %Y")
inspected_count = len(inspection["sources"])
medium_count = sum(
    r["medium_resolution_status"] == "CONDITIONAL_PROFILE_FITS" for r in inspection["sources"]
)
full = next(r for r in cliff["cases"] if not r["mask_core_400_km_s"] and r["lsf_scale"] == 1.3)
masked = next(r for r in cliff["cases"] if r["mask_core_400_km_s"] and r["lsf_scale"] == 1.3)


def paragraph(c, text, x, top, width, *, font=26, leading=35, color=INK, max_height=None):
    style = ParagraphStyle(
        "body", fontName="Helvetica", fontSize=font, leading=leading, textColor=color
    )
    p = Paragraph(text, style)
    _, height = p.wrap(width, H)
    if max_height is not None and height > max_height:
        raise ValueError(f"Poster text overflow ({height}>{max_height}): {text[:60]}")
    p.drawOn(c, x, top - height)
    return top - height


def heading(c, text, x, y, width=700):
    c.setFillColor(TEAL)
    c.setFont("Helvetica-Bold", 28)
    c.drawString(x, y, text)
    c.setStrokeColor(HexColor("#BFCFD2"))
    c.setLineWidth(1.2)
    c.line(x, y - 16, x + width, y - 16)


def header(c, line1, line2, number):
    c.setFillColor(INK)
    c.rect(0, H - 275, W, 275, fill=1, stroke=0)
    c.setFillColor(white)
    c.setFont("Helvetica-Bold", 74)
    c.drawString(70, H - 95, line1)
    c.drawString(70, H - 177, line2)
    c.setFont("Helvetica", 25)
    c.drawString(
        74,
        H - 233,
        f"Biswajit Jana  |  JWST Cosmic Anomaly Lab  |  Research progress: {snapshot_label}",
    )
    c.setFont("Helvetica-Bold", 27)
    c.drawRightString(W - 70, H - 70, f"CONFERENCE POSTER {number:02d}")
    c.setFont("Helvetica", 24)
    c.drawRightString(W - 70, H - 112, "CEERS / JADES / RUBIES / UNCOVER")
    c.drawRightString(W - 70, H - 153, "Public observations. Reproducible analysis.")
    c.drawRightString(W - 70, H - 194, "No new black-hole-star identification.")


def footer(c, text):
    c.setStrokeColor(HexColor("#BFCFD2"))
    c.line(70, 84, W - 70, 84)
    paragraph(c, text, 70, 67, W - 140, font=18, leading=23, max_height=54, color=MUTED)
    c.linkURL(
        "https://github.com/Biswajit1999/JWST-Cosmic-Anomaly-Lab", (70, 15, W - 70, 75), relative=0
    )


def embed_and_save(base, destination, plots):
    page = PdfReader(base).pages[0]
    for path, x, y, width, height in plots:
        plot = PdfReader(ROOT / path).pages[0]
        pw, ph = float(plot.mediabox.width), float(plot.mediabox.height)
        scale = min(width / pw, height / ph)
        dx, dy = x + (width - pw * scale) / 2, y + (height - ph * scale) / 2
        page.merge_transformed_page(plot, Transformation().scale(scale).translate(dx, dy))
    writer = PdfWriter()
    writer.add_page(page)
    writer.add_metadata(
        {
            "/Author": "Biswajit Jana",
            "/Title": destination.stem.replace("_", " "),
            "/Subject": "JWST observational research progress; no new physical classification",
        }
    )
    with destination.open("wb") as f:
        writer.write(f)
    assert len(PdfReader(destination).pages) == 1


# Poster 1: scope, purpose, fixed-sample search and actual observed evidence.
base = TMP / "purpose_base.pdf"
c = canvas.Canvas(str(base), pagesize=(W, H))
header(c, "SEARCHING FOR", "BLACK-HOLE STARS", 1)
heading(c, "01 / THE RESEARCH QUESTION", 70, 1350)
paragraph(
    c,
    "Can unusual JWST sources be explained by stars, ordinary or obscured active galactic nuclei, "
    "dense accretion envelopes, or a mixture of these? The aim is to distinguish these possibilities "
    "with observations and competing predictions.",
    70,
    1306,
    700,
    max_height=220,
)
heading(c, "02 / A DEFINED DATASET", 70, 1050)
paragraph(
    c,
    f"<b>{summary['analysed']:,} / {summary['eligible_denominator']:,} eligible source records analysed.</b><br/>"
    "PRISM spectra; archive grade 3; median pixel S/N at least 5; 0.5 &lt;= z &lt; 12. "
    "One highest-S/N observation per upstream source ID. No colour or physical-class filter.",
    70,
    1006,
    700,
    max_height=230,
)
paragraph(
    c,
    "This is 13.1% of the parent catalogue's 20,144 source records, not a complete JWST census. "
    "Targeting, depth and redshift success affect selection.",
    70,
    770,
    700,
    font=23,
    leading=30,
    max_height=125,
    color=MUTED,
)
heading(c, "03 / FROM ARCHIVE TO EVIDENCE", 70, 600)
paragraph(
    c,
    "<b>1.</b> Verify units, quality masks and product checksums.<br/>"
    "<b>2.</b> Measure continuum windows and empirical line blends.<br/>"
    "<b>3.</b> Vary instrument resolution and assumed errors.<br/>"
    "<b>4.</b> Inspect images, 2D spectra and higher-resolution data.<br/>"
    "<b>5.</b> Check published benchmarks and competing explanations.",
    70,
    556,
    700,
    font=25,
    leading=34,
    max_height=340,
)
c.setFillColor(PALE)
c.rect(70, 116, 700, 100, fill=1, stroke=0)
paragraph(
    c,
    f"<b>{summary['joint_followup_priorities']} follow-up priorities.</b><br/>"
    "A screening pass is not a physical identification.",
    90,
    197,
    660,
    font=25,
    leading=32,
    max_height=75,
)
heading(c, "ACTUAL TELESCOPE DATA", 900, 1350, 1400)
footer(
    c,
    "Data credit: NASA / ESA / CSA / STScI; CEERS, JADES, RUBIES, UNCOVER; DAWN JWST Archive / msaexp. "
    "Images are false-colour NIRCam cutouts; PSFs are unmatched. Full methods, data hashes and code: "
    "github.com/Biswajit1999/JWST-Cosmic-Anomaly-Lab",
)
c.save()
embed_and_save(
    base,
    OUT / "jwst_project_purpose_methods_A1.pdf",
    [
        ("figures/observed_image_comparison.pdf", 960, 860, 1240, 465),
        ("figures/expanded_search.pdf", 865, 115, 1450, 705),
    ],
)

# Poster 2: a focused, falsifiable benchmark result and its limitations.
base = TMP / "evidence_base.pdf"
c = canvas.Canvas(str(base), pagesize=(W, H))
header(c, "FROM A CANDIDATE", "TO A PHYSICAL CLAIM", 2)
heading(c, "THE CLIFF / RUBIES-UDS-154183", 70, 1350, 1440)
heading(c, "01 / WHAT THE DATA PREFER", 1590, 1350, 715)
paragraph(
    c,
    "The Cliff is a <b>published candidate</b>, added separately as a benchmark. "
    "Our G395M fit compares a Gaussian broad component with a Gaussian convolved with symmetric "
    "exponential tails. Narrow nuisance lines and a local continuum are included.",
    1590,
    1306,
    715,
    max_height=245,
)
c.setFillColor(TEAL)
c.setFont("Helvetica-Bold", 65)
c.drawString(1590, 1010, f"{full['delta_bic_gaussian_minus_exponential']:.1f}")
c.drawString(1970, 1010, f"{masked['delta_bic_gaussian_minus_exponential']:.1f}")
paragraph(c, "Full-window delta BIC", 1590, 986, 340, font=22, leading=28, max_height=60)
paragraph(c, "With four core pixels masked", 1970, 986, 335, font=22, leading=28, max_height=60)
heading(c, "02 / THE SENSITIVITY TEST", 1590, 875, 715)
paragraph(
    c,
    "Positive delta BIC favours exponential tails. Across LSF tests, the full-window value is "
    "56.8-69.6, but it falls to 3.5-9.4 after excluding the central +/-400 km/s. With doubled "
    "errors, the masked comparison no longer favours the exponential family.",
    1590,
    831,
    715,
    max_height=230,
)
heading(c, "03 / WHAT THIS DOES NOT SHOW", 1590, 565, 715)
paragraph(
    c,
    "The preference depends on the line core, where absorption is unmodelled. "
    "Empirical shapes cannot uniquely identify an accretion envelope. No black-hole mass, "
    "density or physical-class probability is inferred here.",
    1590,
    521,
    715,
    max_height=195,
)
heading(c, "04 / THE NEXT PHYSICAL TEST", 1590, 265, 715)
paragraph(
    c,
    "Fit absorption and emission together; validate extraction and calibration; compare "
    "stellar, AGN and dense-envelope predictions across multiple independent observables.",
    1590,
    221,
    715,
    font=24,
    leading=32,
    max_height=130,
)
paragraph(
    c,
    "<b>A missed-benchmark check:</b> the PRISM continuum screen finds The Cliff, but the fixed-width "
    f"broad-blend screen misses it. The {summary['joint_followup_priorities']}-source queue is therefore incomplete for published "
    f"candidates. {inspected_count} sources received detailed inspection; {medium_count} yielded conditional "
    "medium-resolution fits. No new source is established as a black-hole star.",
    70,
    245,
    1440,
    font=25,
    leading=34,
    max_height=145,
)
footer(
    c,
    "Benchmark: de Graaff et al. (2025), A&A 701, A168, arXiv:2503.16600. Observations: RUBIES / JWST; "
    "reduction: DJA / msaexp. This project performs a conditional reanalysis, not an independent reduction. "
    "Results and executable methods: github.com/Biswajit1999/JWST-Cosmic-Anomaly-Lab",
)
c.save()
embed_and_save(
    base,
    OUT / "jwst_evidence_limits_A1.pdf",
    [
        ("figures/expanded_sources/CAL-ebc2a6a997b3535d_medium.pdf", 70, 285, 1440, 1020),
    ],
)
manifest = dict(
    format="A1 landscape, 841 x 594 mm",
    audience="research/conference",
    date=snapshot_date.isoformat(),
    author="Biswajit Jana",
    build_script="scripts/build_research_posters.py",
    input_hash_convention="UTF-8 text with universal-newline normalization",
    result_inputs={
        name: hashlib.sha256((ROOT / name).read_text(encoding="utf-8").encode("utf-8")).hexdigest()
        for name in [
            "results/expanded/summary.json",
            "results/cliff_sensitivity.json",
            "results/expanded/inspection_audit.json",
        ]
    },
    outputs=[p.relative_to(ROOT).as_posix() for p in sorted(OUT.glob("*A1.pdf"))],
    interpretation="Research progress snapshot; no new black-hole-star identification",
)
(OUT / "poster_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
print(json.dumps(manifest, indent=2))
