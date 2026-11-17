COBS data used in the paper

Source: COBS Comet Observation Database, https://cobs.si (retrieved 2026 November 16; the records
of the second ten comets 2026 November 17).
Licence: CC BY-NC-SA 4.0 (unlike the CC BY code and text of this supplementary material).
Credit: "COBS Comet Observation Database - CC BY-NC-SA 4.0". We thank the observers worldwide
who contributed to COBS.

un271_cobs.csv
  The individual CCD records of comet C/2014 UN271 used by cobs_un271.py (Sec. IV A, Fig. 3).

cobs_lightcurves.txt
  Binned light curves used by forecast_test.py (Sec. III C), one section per comet (#key), lines
  "date  magnitude  error":
  - hb (C/1995 O1), hyak (C/1996 B2), neo (C/2020 F3), p12 (12P), a3 (C/2023 A3), i3 (3I/ATLAS),
    t5 (C/2024 T5), y1 (C/2025 Y1), sw1 (29P): total magnitudes from visual estimates (methods S,
    M, B, I), CCD visual-equivalent (Z), CCD V (V) and unfiltered CCD (C) observations, median per
    bin (15 days for Hale-Bopp, 12 for 12P and C/2023 A3, 9 for C/2025 Y1, 8 for C/2024 T5, 5 for
    3I/ATLAS, 3 for Hyakutake and NEOWISE, 52 for 29P), at least 2 per bin, dated at the median
    time. Error: sqrt(0.15^2 + (1.4826 MAD)^2/N); bins with an error above 0.6 mag are dropped.
    No aperture correction: they are compared with a fixed aperture of 20 000 km (j = 0).
  - q2 (C/2014 Q2), k2 (C/2017 K2), t2 (C/2017 T2), p1 (C/2009 P1), e3 (C/2022 E3), l3 (C/2019 L3),
    l4 (C/2011 L4), a6 (C/2025 A6), v2j (C/2015 V2), v2z (C/2020 V2): the second ten comets of the
    forecast test, chosen by a rule fixed before their forecasts were computed (the comets not yet
    in the sample with the most COBS observations in one apparition among comets with perihelion in
    2010-2025). All records of the apparition, same methods, binned by bin_cobs.py with one rule:
    bins of max(3 d, span/50) counted from the first record, otherwise as above. A "%" line before
    each section gives the number of records and the bin width. bin_cobs.py --fetch downloads the
    records again (internet needed; COBS keeps adding records, so the bins may change slightly).
  - un271: monthly medians (at least 3 values) of the C/2014 UN271 records, moved to r' in a
    6" radius aperture as in cobs_un271.py.
