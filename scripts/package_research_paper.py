"""Package the compiled preprint with its exact LaTeX and figure inputs."""

import json
import re
import zipfile
from pathlib import Path

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
source = root / "paper/jwst_cosmic_anomaly_lab.tex"
pdf = root / "output/pdf/jwst_cosmic_anomaly_lab.pdf"
if not pdf.is_file() or not pdf.read_bytes().startswith(b"%PDF-"):
    raise ValueError("Compile the LaTeX PDF before packaging")
figures = [root / "figures" / name for name in re.findall(r"\\fig\{([^}]+)\}", source.read_text())]
assert len(figures) == 7 and all(path.is_file() for path in figures)
files = [source, root / "paper/README.md", root / "paper/linkedin_draft.md", *figures]
manifest = dict(created_utc=utcnow(), research_snapshot="2026-11-16", compiler="Tectonic 0.17.0",
    prior_analysis_commit="6c9a7144fcafae72f397d209ad319a01a3400225",
    status="RESEARCH_PREPRINT_NOT_PEER_REVIEWED", physical_classification="NOT_ESTABLISHED",
    figure_count=len(figures), input_hash_convention="Raw bytes for PDFs; canonical UTF-8/LF for text",
    input_sha256={p.relative_to(root).as_posix(): (digest(p) if p.suffix == ".pdf" else text_digest(p)) for p in files},
    output=dict(path=pdf.relative_to(root).as_posix(), bytes=pdf.stat().st_size, sha256=digest(pdf)),
    analysis_summary_sha256={p: text_digest(root / p) for p in [
        "results/multiwavelength/dja_summary.json", "results/multiwavelength/context_summary.json"]})
write_json(root / "paper/build_manifest.json", manifest)
bundle = root / "output/jwst_cosmic_anomaly_lab_latex_source.zip"
with zipfile.ZipFile(bundle, "w", compression=zipfile.ZIP_DEFLATED) as archive:
    for path in [*files, pdf, root / "paper/build_manifest.json"]:
        archive.write(path, path.relative_to(root).as_posix())
with zipfile.ZipFile(bundle) as archive:
    assert archive.testzip() is None
print(json.dumps(dict(bundle=str(bundle), bytes=bundle.stat().st_size,
                      bundled_files=len(files)+2, pdf_bytes=pdf.stat().st_size)))
