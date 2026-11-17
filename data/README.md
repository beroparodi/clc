# Light-curve data bundled with the tool

Source: COBS Comet Observation Database, https://cobs.si (Observation List API, retrieved
2026 November 17; C/2014 UN271 2026 November 16). Licence: **CC BY-NC-SA 4.0** — this folder and the light curves built into the
page are shared under that licence, unlike the MIT-licensed code. Credit: COBS Comet Observation
Database – CC BY-NC-SA 4.0. We thank the observers worldwide who contributed to COBS.

`cobs_lightcurves.txt` — one section per comet (`#key`), lines `date  magnitude  error`:

* all comets except `un271` (the nine other comets of the tool and the fifteen added to the paper's
  forecast test: `q2`, `k2`, `t2`, `p1`, `e3`, `l3`, `l4`, `a6`, `v2j`, `v2z`, the ten most observed
  comets with perihelion in 2010–2025, and `us10`, `k1`, `g2`, `e2`, `n2`, the five most observed of
  that period first seen beyond 5 au): total magnitudes from visual estimates (methods S, M, B, I),
  CCD visual-equivalent (Z), CCD V (V) and unfiltered CCD (C) observations of the apparition (12P from
  2020, 29P 2010–2026), retrieved 2026 November 17 and binned by one rule (`code/bin_cobs.py`):
  records between the 1st and 99th percentile of their dates, bins of max(3 d, span/50), median of at
  least 2 values dated at the median time, error sqrt(0.15² + (1.4826 MAD)²/N), bins with an error
  above 0.6 mag dropped. A `%` line before each section gives the record count and bin width. No
  aperture correction: the tool compares them with a fixed aperture of 20 000 km (j = 0).
* `un271`: total CCD magnitudes (V, R and V-calibrated unfiltered) moved to the paper's system,
  r′ in a 6″ radius aperture: the stated photometric aperture is read as a diameter (it equals
  the coma diameter in 301 of 305 records), the nucleus is removed, the coma flux is scaled by
  6″/ρ, and V − r′ = 0.22, r′ − R = 0.18; monthly medians of at least 3 values (code:
  `cobs_un271.py` in the paper's supplementary material). Not used in the paper's fit; the two
  professional points of the paper are kept separately in the tool.

`un271_cobs_raw.txt` — the individual C/2014 UN271 records: date, magnitude, method, aperture
size, unit (m = arcmin, s = arcsec), shape, observer code, telescope (cm), reference catalogue,
type (C = CCD, V = visual).
