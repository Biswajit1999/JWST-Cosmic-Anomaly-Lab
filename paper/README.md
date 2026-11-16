# Research paper, snapshot 16 November 2026

This is a research preprint/project methods report, **not peer reviewed**. It
consolidates the executed project, including the new cross-reduction and
Hubble/Spitzer/Chandra catalogue checks. It does not announce a newly discovered
black-hole star.

The editable source is `jwst_cosmic_anomaly_lab.tex`. References are embedded in
the LaTeX source; no separate BibTeX installation is required. Figures are real
observations, derived measurements or explicitly identified synthetic controls.

## Build

From this directory with Tectonic 0.17.0 available:

```powershell
tectonic --keep-logs --outdir ../output/pdf jwst_cosmic_anomaly_lab.tex
```

Alternatively, run `pdflatex` twice from this directory with a standard TeX Live
installation. The first Tectonic run downloads its TeX resources; later runs can
reuse the cache. The downloaded executable in `.tools/` is intentionally not
version controlled. Its official release is:
https://github.com/tectonic-typesetting/tectonic/releases/tag/tectonic%400.17.0

The source bundle preserves the `paper/` and `figures/` relative paths. Extract
the ZIP, enter `paper/`, and use the command above after creating `output/pdf/`.
The bundled build manifest records the exact input/output hashes.

## Scientific scope

Prior programme results are pinned by commit `6c9a714`. New execution records
are in `experiments/multiwavelength_validation.yaml`; raw retrieval provenance
is in `data/manifests/downloads.json`. The paper reports both successful and
unsuccessful checks. A local test-suite pass does not constitute external
scientific review.

For LinkedIn, attach the PDF and use or edit `linkedin_draft.md`. No social-media
post has been sent automatically. Historical conference posters describe their
dated results; the current paper supersedes their shortlist interpretation.
