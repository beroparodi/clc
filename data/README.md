# Light-curve data bundled with the tool

Source: COBS Comet Observation Database, https://cobs.si (Observation List API, retrieved
2026 November 16). Licence: **CC BY-NC-SA 4.0** — this folder and the light curves built into the
page are shared under that licence, unlike the MIT-licensed code. Credit: COBS Comet Observation
Database – CC BY-NC-SA 4.0. We thank the observers worldwide who contributed to COBS.

`cobs_lightcurves.txt` — one section per comet (`#key`), lines `date  magnitude  error`:

* `hb`, `hyak`, `neo`, `p12`, `a3`, `i3`, `t5`, `y1`, `sw1`: total magnitudes from visual
  estimates (methods S, M, B, I), CCD visual-equivalent (Z), CCD V (V) and unfiltered CCD (C)
  observations, median per bin (15 days for Hale–Bopp, 12 for 12P and Tsuchinshan–ATLAS, 9 for
  C/2025 Y1, 8 for C/2024 T5, 5 for 3I/ATLAS, 3 for Hyakutake and NEOWISE, 52 for 29P), at least
  2 per bin, dated at the median time. Error: sqrt(0.15² + (1.4826 MAD)²/N); bins with an error
  above 0.6 mag are dropped. No aperture correction: the tool compares them with a fixed
  aperture of 20 000 km (j = 0), and K in each preset is fitted to them (`test/fit_cobs.js`).
* `q2`, `k2`, `t2`, `p1`, `e3`, `l3`, `l4`, `a6`, `v2j`, `v2z`: the second ten comets of the paper's
  forecast test, chosen by a rule fixed before their forecasts were computed (the ten comets not yet in
  the sample with the most COBS observations in one apparition among comets with perihelion in
  2010–2025; census of 2026 November 17). All records of the apparition (methods S, M, B, I, Z, V, C,
  retrieved 2026 November 17) binned by one rule: bins of max(3 d, span/50) counted from the first
  record, median of at least 2 values dated at the median time, the same error and the same 0.6 mag cut
  as above (`code/bin_cobs.py`; a `%` line before each section gives the record count and bin width).
* `un271`: total CCD magnitudes (V, R and V-calibrated unfiltered) moved to the paper's system,
  r′ in a 6″ radius aperture: the stated photometric aperture is read as a diameter (it equals
  the coma diameter in 301 of 305 records), the nucleus is removed, the coma flux is scaled by
  6″/ρ, and V − r′ = 0.22, r′ − R = 0.18; monthly medians of at least 3 values (code:
  `cobs_un271.py` in the paper's supplementary material). Not used in the paper's fit; the two
  professional points of the paper are kept separately in the tool.

`un271_cobs_raw.txt` — the individual C/2014 UN271 records: date, magnitude, method, aperture
size, unit (m = arcmin, s = arcsec), shape, observer code, telescope (cm), reference catalogue,
type (C = CCD, V = visual).
