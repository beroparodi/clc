# How bright will a comet get? A physical magnitude law

The interactive tool and code accompanying

> B. R. Parodi, *How bright will a comet get? A physical magnitude law* (in preparation)

## Light Up a Comet (interactive)

Open **https://beroparodi.github.io/clc/** in any browser, or download `index.html` and open it locally (no installation needed).
Pick a comet, the ice that drives it (CO, CO<sub>2</sub> or water), or two ices on separate patches, and how the Sun’s heat spreads over the nucleus, from the subsolar point (⟨cos ζ⟩ = 1) to an isothermal, fast-rotating body (⟨cos ζ⟩ = ¼). The page solves the power balance of the standard model (Eqs. 8–10) for the share η(r) of the absorbed sunlight that sublimates the ice (Eq. 11), gives the production index *n*<sub>Q</sub> and the activity index of the coma *n*<sub>C</sub> = 3/2 + *n*<sub>Q</sub> (Eqs. 7, 12), and the magnitude in a fixed angular or physical aperture with the physical law, Eq. (14):

*m* = *m*<sub>1</sub> + 5 log Δ + 8.75 log *r* − 2.5 log Φ(*r*), Φ = Σ *k*<sub>i</sub> η<sub>i</sub>(*r*) + ν *r*<sup>3/2</sup>.

In an aperture of fixed angle the coma fades as Δ<sup>−1</sup>, not as Δ<sup>−2</sup>.

- **Light curve**: model, nucleus alone, your photometry, perihelion and today; drag the date line to read the model at any date (*r*, Δ, phase angle, magnitude, nucleus share, surface temperature of each ice, Σ*k*η, *n*<sub>Q</sub> and *n*, Afρ).
- **Afρ growth**: the growth of Afρ from a reference date for all ices; works with *q*, *e* and *T* alone (as for C/2025 Y1).
- **Fit and compare**: fit the coma scale *K* (and optionally *H*<sub>N</sub>) to your photometry; with two ices, also their dust share *k* (Sec. IV of the paper). “Compare ices” fits all five single-ice cases and three pairs of ices at once and lists χ²<sub>ν</sub>.
- **Before and after perihelion**: fit both laws separately to the two arcs, as in Sec. V of the paper; the fit reports how the active area (and, for two ices, *k*) changed at perihelion.
- **Forecasts against the empirical law**: the light curve also shows the empirical law *m* = *m*<sub>0</sub> + 5 log Δ + 2.5 *n* log *r* (Eq. 2), with *m*<sub>0</sub> and *n* fitted to the same points. With a “fit only before” date, both are fitted to the earlier points only and the later ones test the two forecasts (rms in the legend), as in Fig. 4 of the paper. In our test of 24 COBS light curves, the physical law beats the empirical law with fitted *n* for 7 of the 9 comets whose ice switches on within the fitted arc, against 3 of the other 15; a fixed *n* = 3.5, the law’s value for a fully active ice, does about as well: the gain comes mostly from knowing that a steep rise cannot last. The fifteen comets added to that test by fixed rules are built in with the split used there.
- **Physics panels**: η(*r*) for the three ices and the coma index *n*<sub>C</sub>(*r*) (Fig. 2 of the paper), for two ices also Σ*k*η and its index, with the empirical fit of Holt et al. (2026) and the measured index of C/2014 UN271.
- **Backtests**: the power budget of a nucleus (NEOWISE, 3I/ATLAS or your own numbers) and the SOHO/SWAN water-production exponents against the constant-area band.

Built-in comets: from the paper, C/2014 UN271 (Bernardinelli–Bernstein) with the photometry used there (Fig. S2b of the supplementary material), C/1995 O1 (Hale–Bopp) fitted before and after perihelion (Fig. S1b of the supplementary material; forecast of Fig. 4a), 12P/Pons–Brooks with the forecast of Fig. 4b, and C/2019 L3 (ATLAS) with CO and water ice on separate patches (Fig. S1c of the supplementary material); further, C/2024 T5 (ATLAS), C/2025 Y1 (ATLAS), 1P/Halley (return of 2061), 2P/Encke (2027), 29P/Schwassmann–Wachmann, 67P/Churyumov–Gerasimenko (2028), C/1996 B2 (Hyakutake), C/2020 F3 (NEOWISE), C/2023 A3 (Tsuchinshan–ATLAS) and the interstellar 3I/ATLAS; and the fifteen comets added to our forecast test, C/2014 Q2 (Lovejoy), C/2017 K2 (PANSTARRS), C/2017 T2 (PANSTARRS), C/2009 P1 (Garradd), C/2022 E3 (ZTF), C/2019 L3 (ATLAS), C/2011 L4 (PANSTARRS), C/2025 A6 (Lemmon), C/2015 V2 (Johnson), C/2020 V2 (ZTF), C/2013 US10 (Catalina), C/2012 K1 (PANSTARRS), C/2021 G2 (ATLAS), C/2022 E2 (ATLAS) and C/2022 N2 (PANSTARRS); and two famous disappointments, C/1973 E1 (Kohoutek) and C/2012 S1 (ISON), run through the same forecast rules (they are not among the 24 comets); with orbits from the JPL Small-Body Database (Hale–Bopp: NASA NSSDC). For these, the coma scale K (and H<sub>N</sub> where no nucleus size is known) is a placeholder to be fitted to photometry. “Your comet” takes any orbit from the JPL Small-Body Database or the MPC.

## Code (`code/`)

Python 3 (numpy, scipy). Keep the files in one folder and run each script from there. These scripts were written for an earlier, research version of the paper; their section, equation and note numbers refer to that version.

| file | purpose |
|---|---|
| `un271_model.py` | energy balance (Eqs. 8–11), nucleus + coma model (Eqs. 2–7), analysis of C/2014 UN271 (Sec. IV A), worked examples of the notes, Monte Carlo slope *k* |
| `un271_orbit.py` | two-body orbit of C/2014 UN271 and the Earth |
| `backtest.py` | backtest (Sec. III): SOHO/SWAN exponents and normalization, with and without a dust mantle, C/2017 K2, NEOWISE, 3I/ATLAS, long-period comets (Lacerda et al. 2025; Holt et al. 2026) |
| `swan_data.py` | SWAN water-production power laws of Combi et al. (2019) and later papers |
| `predictions.py` | forecasts: C/2014 UN271 (all cases, pre-discovery test, Eq. 14), C/2025 Y1, C/2024 T5 |
| `cobs_un271.py` | test of the C/2014 UN271 forecast with COBS photometry (Sec. IV A, Fig. 3) |
| `forecast_test.py` | forecasts of Eq. (14) against the empirical law for 25 COBS light curves, with the switch-on flag (Sec. III C) |
| `bin_cobs.py` | fetches COBS records and bins them with the rule of the forecast test |
| `cobs/` | the COBS data these two scripts read (CC BY-NC-SA 4.0, see `cobs/README.txt`) |
| `earlier_model_revisited.py` | the six comets of the earlier model, Parodi (2013) |

## Page sources (`src/`)

`index.html` is assembled by `python build.py` from `src/page.html` (markup and style), `src/model.js` (power balance, Eqs. 8–10, the physical law, Eq. 14, orbits), `src/swan.js` (SWAN data) and `src/ui.js` (page logic). `node test/test_model.js` checks the JavaScript model against the numbers of the Python code.

## Data

C/2014 UN271 photometry: Kelley et al. 2022 (ApJL 933, L44) and Roth et al. 2025 (ApJL 986, L22); orbit from the JPL Small-Body Database. Orbits of the other built-in comets: JPL Small-Body Database (via TheSkyLive and Wikipedia), Hale–Bopp from NASA’s NSSDC; nucleus sizes of 1P/Halley (Giotto), 2P/Encke (Fernández et al. 2000), 29P (Schambeau et al. 2015), 67P (Rosetta), Hale–Bopp (Szabó et al. 2012, infrared) and Hyakutake (Lisse et al. 1999, infrared), each with the albedo stated in the tool. SWAN exponents: Combi et al. 2019 (Icarus 317, 610), Combi et al. 2021 (six comets 2017–2020, Icarus 365, 114509) and Combi et al. 2025 (Icarus 438, 116645; K2 and E3 exponents fitted to their tables). Activity index fit: Holt et al. 2026. NEOWISE: Combi et al. 2021; 3I/ATLAS: Tan et al. 2026, Hui et al. 2026. Light curves of C/2014 UN271, C/2024 T5, C/2025 Y1, Hale–Bopp, Hyakutake, NEOWISE, 12P, Tsuchinshan–ATLAS, 3I/ATLAS and 29P: binned magnitudes from the COBS Comet Observation Database (https://cobs.si), retrieved 2026 November 16; all except C/2014 UN271 re-retrieved 2026 November 17 and binned by `bin_cobs.py` with one rule; see `data/README.md` for the binning and the conversion of the UN271 magnitudes. **The COBS data (folders `data/` and `code/cobs/` and the light curves built into the page) are licensed CC BY-NC-SA 4.0**, unlike the MIT-licensed code. Elements of C/2025 Y1: JPL Small-Body Database (solution of 2026 September 8); of the fifteen added comets: JPL Small-Body Database, retrieved 2026 November 17. Please cite the original papers when using these data, and credit COBS: “COBS Comet Observation Database – CC BY-NC-SA 4.0”.

## How to cite

See [`CITATION.cff`](CITATION.cff).

## License

Code: MIT (see `LICENSE`). COBS data in `data/` and in the built page: CC BY-NC-SA 4.0 (see `data/README.md`).
