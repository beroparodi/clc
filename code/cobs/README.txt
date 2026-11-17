COBS data used in the paper

Source: COBS Comet Observation Database, https://cobs.si (retrieved 2026 November 17; the records of
C/2014 UN271 2026 November 16).
Licence: CC BY-NC-SA 4.0 (unlike the CC BY code and text of this supplementary material).
Credit: "COBS Comet Observation Database - CC BY-NC-SA 4.0". We thank the observers worldwide
who contributed to COBS.

un271_cobs.csv
  The individual CCD records of comet C/2014 UN271 used by cobs_un271.py (Sec. IV A, Fig. 3).

cobs_lightcurves.txt
  Binned light curves used by forecast_test.py (Sec. III C), one section per comet (#key), lines
  "date  magnitude  error"; a "%" line before each section gives the number of records and the bins.
  - The 24 comets other than C/2014 UN271 (keys in forecast_test.py: the nine other comets of the
    interactive tool, the ten most observed comets with perihelion in 2010-2025 and the five most
    observed of that period first seen beyond 5 au): total magnitudes from visual estimates
    (methods S, M, B, I), CCD visual-equivalent (Z), CCD V (V) and unfiltered CCD (C) observations
    of the apparition
    (12P from 2020, 29P 2010-2026), retrieved 2026 November 17, binned by bin_cobs.py with one rule:
    records between the 1st and 99th percentile of their dates, bins of max(3 d, span/50) counted
    from the first, median of at least 2 values dated at the median time, error
    sqrt(0.15^2 + (1.4826 MAD)^2/N), bins with an error above 0.6 mag dropped. No aperture
    correction: they are compared with a fixed aperture of 20 000 km (j = 0).
  - un271: monthly medians (at least 3 values) of the C/2014 UN271 records, moved to r' in a
    6" radius aperture as in cobs_un271.py.
