# Observed visual evidence: methods

The atlas extends image coverage of the existing 20-source pilot. It does not
rerun its population search or upgrade any source's physical classification.
The original `results/run_summary.json` describes the original pilot, including
its three image sources. `results/visual_atlas_run.json` describes this separate
visual/imaging extension and contains the new aperture measurements.

## Inputs and reproducibility

Run from the repository root after installation and pilot spectrum download:

```sh
python scripts/download_cutouts.py --all
python scripts/build_visual_atlas.py
cosmic-lab verify
```

The first command uses the existing 512 MiB single-writer download cache.
`data/manifests/downloads.json` records exact URLs, dates and SHA256 hashes.
The cutout service is live, not a frozen survey release. To reproduce these
specific products, use `scripts/replay_manifest.py`, which checks the stored
hashes and rejects upstream changes. The atlas manifest identifies the selected
catalogue, renderer and individual spectrum/image hashes.

The selection comes from `data/catalogues/pilot_selection.csv`, unchanged by
this extension. Images are DJA cutouts of public NIRCam mosaics at those catalogue
coordinates. A matching coordinate does not prove an unblended counterpart.
Of 20 cutouts, 18 provide usable F150W, F277W and F444W images; JADES 34958
(`CAL-8e9bf6bf1be95f7a`) contains only F444W in this service response. This is
missing coverage, not a blue-band upper limit or a measured dropout. JADES 3038
(`CAL-b1486215c5645b42`) returns an F277W extension without sufficient valid
background pixels; its other usable bands remain visible, but RGB is withheld.

## Image rendering

SCI pixels are converted to µJy/pixel using each HDU's PHOTFNU (Jy per native
unit). Nonfinite SCI/WHT or WHT ≤ 0 pixels are masked. The display background is
the sigma-clipped median in a 1.5–2.5 arcsecond annulus. This is an approximate
local sky estimate; crowding and diffuse emission can bias it.

All pixel centres must agree between band WCS grids within 0.05 pixel before
composition. No interpolation or sharpening is applied. The images retain the
native archive orientation. Each cyan circle is centred on the catalogue WCS
position with 0.4 arcsecond radius; it is a locator, not a measured source size
or photometric aperture. Scale bars span one arcsecond using the WCS pixel scale.

False colour assigns F444W to red, F277W to green and F150W to blue. A common
intensity-dependent asinh transform with Q=8 preserves positive channel ratios
at each pixel. Display scale is the 99.7th percentile of the positive mean
channel intensity inside one arcsecond of the target. Saturated pixels divide
all channels by the same factor. Negative values are clipped only for RGB
display; invalid RGB pixels are grey. No missing channel is synthesized.

Each source has its own intensity scale, so displayed brightness cannot be
compared across sources. Monochrome panels use independent asinh stretches and
retain signed background-subtracted values; their apparent brightness cannot
be compared across filters. The manifest records calibration, sky, mask coverage
and RGB scale. False colour is not human-eye visible colour.

## Spectral evidence

The 1D plot uses the released SPEC1D flux and error columns with their FITS units
and existing quality mask, with no extra flux rescaling. Gaps remain gaps. The
shaded uncertainty is conditional statistical ±1σ; calibration, extraction,
redshift and correlated-error systematics are not included.

The continuum detail uses the existing rest-frame 0.34–0.36 and 0.40–0.42 µm
windows and reports their red/blue fν ratio only when coverage and significance
requirements pass. The ratio is not a physical-model likelihood or a uniquely
stellar Balmer break. Archive redshifts are used without propagated uncertainty.

The lower panel displays the actual released 2D SCI extension in µJy, masking
nonfinite or nonpositive-weight pixels. It uses the SPEC1D wavelength centres
and their midpoint bin edges, and YTRACE for the spatial origin. Pixel areas
follow the nonuniform wavelength grid. The asinh display and clipping are for
visibility; this is not a pixel significance map. It is a rectified archive
image, not raw detector counts. It makes neighbours, residuals and trace
structure inspectable, but no automated 2D extraction validation is claimed.

## Aperture diagnostics and limits

The extension also runs the existing 0.3 arcsecond aperture measurement in every
available band. Its sky annulus is 0.8–1.2 arcseconds, distinct from the image
display background. Signed fluxes, diagonal errors, flags and unavailable
measurements are retained in the atlas manifest. These are neither PSF-matched
total photometry nor physical SED fits. Drizzle covariance, aperture corrections,
deblending, lens magnification and differential magnification remain unmodelled.

The images and 2D spectra can guide further checks; they cannot establish a
dense accretion envelope, electron-scattering wings or a black-hole mass.
