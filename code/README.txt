Supplementary material: "Forecasting comet light curves from sublimation physics"

Python 3 code (numpy and scipy; tested with Python 3.13, numpy 2.5, scipy 1.18) that reproduces the
numbers, tables and figures of the paper. Keep all files in one folder and run each script from
there. EXERCISES.txt suggests how to use the code and the data in a course.

  un271_model.py           Energy balance (Eqs. 8-11) at the subsolar point, on the sunlit
                           hemisphere of a slow rotator and for an isothermal body; nucleus +
                           coma model (Eqs. 2-7); analysis of comet C/2014 UN271 (Sec. IV A),
                           worked examples of the notes, Monte Carlo slope k
  un271_orbit.py           Two-body orbit of C/2014 UN271 and the Earth (note on the orbit)
  backtest.py              Backtest (Sec. III): SOHO/SWAN exponents, baselines (fixed p, fixed
                           band), medians by perihelion distance for long- and short-period
                           comets, criterion (ii) for all comets with known radii, dust mantle;
                           comet C/2017 K2 (incl. the CRIRES+ water ratio), comet NEOWISE,
                           3I/ATLAS, long-period comets (Lacerda et al. 2025; Holt et al. 2026)
  swan_data.py             SWAN water-production power laws of Combi et al. (2019) and later
                           papers (2015-2025, incl. 3I/ATLAS); SWAN rates of C/2017 K2 and
                           C/2022 E3 (Combi et al. 2025) and our fits to them
  forecast_test.py         Forecasts of Eq. (14) against the empirical law, Eq. (1), for 25
                           COBS light curves (Sec. III C and its note), with the switch-on flag
                           and the comparison of flagged comets with the rest: fit on the first half,
                           forecast rms on the second half, against Eq. (1) with n fitted and
                           with n fixed at 3.5 and 4; criterion (iii); hindsight fits; forecasts
                           with subsolar water imposed (C/1995 O1, C/2023 A3)
  cobs_un271.py            Test of the C/2014 UN271 forecast with COBS photometry (Sec. IV A,
                           Fig. 4): conversion to r' in 6", offsets for all cases, observers,
                           coma profile from the apertures, time trend of the residuals, the
                           filtered calibration points for both conversions
  predictions.py           Forecasts: C/2014 UN271 (all cases, robustness of the CO2 inference:
                           subsolar vs hemisphere, albedo, emissivity and gas-drag grain speed;
                           pre-discovery test and coma shares in 6" and 1"; one-line formula,
                           Eq. 14), C/2025 Y1 and C/2024 T5 (JPL orbits; check of Barker's
                           parabolic equation)
  earlier_model_revisited.py  Exploration, not used in the paper: the six comets of the earlier
                           model, Parodi (2013) (standalone)
  cobs/un271_cobs.csv      The COBS records of C/2014 UN271 (COBS Comet Observation Database,
                           https://cobs.si; licensed CC BY-NC-SA 4.0, not CC BY)
  bin_cobs.py              The binning rule of the forecast test; --fetch downloads the COBS
                           records of all its comets but C/2014 UN271 (internet needed)
  cobs/cobs_lightcurves.txt  Binned COBS light curves of the 25 comets (median per bin;
                           CC BY-NC-SA 4.0, see cobs/README.txt)

Usage (running times on a laptop):
       python un271_model.py          (seconds)
       python backtest.py             (about 20 s)
       python backtest.py fig2        (coordinates of Fig. 3, the p-q diagram)
       python backtest.py fign        (coordinates of Fig. 2, the activity index)
       python predictions.py          (seconds; prints the coordinates of the hemisphere curve of Fig. 4)
       python cobs_un271.py           (seconds)
       python forecast_test.py        (about a minute)
       python earlier_model_revisited.py

The binned light curves were made from COBS records retrieved on 2026 November 17 (C/2014 UN271:
November 16; see cobs/README.txt). The raw records are not included: bin_cobs.py --fetch downloads
them and bin_cobs.py bins them with the rule of the paper (COBS keeps adding records, so the bins
may change slightly).

License: CC BY 4.0 (code and text); the COBS data in cobs/ are CC BY-NC-SA 4.0.
