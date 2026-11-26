// ---- Page logic for "Light Up a Comet"
const $ = (id) => document.getElementById(id);
const DAY = 86400000;
const ICE_NAME = { CO: "CO", CO2: "CO₂", H2O: "H₂O" };
const ICE_VAR = { CO: "--co", CO2: "--co2", H2O: "--h2o" };
const SCEN = [
  { x: "CO", c: 1, lab: "CO, subsolar" }, { x: "CO2", c: 1, lab: "CO₂, subsolar" }, { x: "CO2", c: 0.25, lab: "CO₂, isothermal" },
  { x: "H2O", c: 1, lab: "H₂O, subsolar" }, { x: "H2O", c: 0.25, lab: "H₂O, isothermal" }];
// two ices on separate patches (Eq. 14), both at the subsolar point; the share k is fitted
const PAIRS = [
  { x: "CO", y: "H2O", c: 1, lab: "CO + H₂O, subsolar" }, { x: "CO2", y: "H2O", c: 1, lab: "CO₂ + H₂O, subsolar" },
  { x: "CO", y: "CO2", c: 1, lab: "CO + CO₂, subsolar" }];
const keyOf = (sc) => sc.x + (sc.y || "") + sc.c;
// ranking of fits: chi2_nu; a fit with as many constants as points cannot win, and a pair must do clearly better
const rankOf = (f, sc) => (f.nu > 0 ? f.score : Infinity) * (sc && sc.y ? 1.001 : 1);
// the ices of the model: one ice, or two on separate patches with dust shares k and 1 - k
function mixFrom(x, y, k) { return !y || y === "none" || y === x ? [{ x, w: 1 }] : [{ x, w: k }, { x: y, w: 1 - k }]; }
function mixOf(P) { return mixFrom(P.ice, P.ice2, P.k); }
function mixName(P) { return P.ice2 && P.ice2 !== "none" && P.ice2 !== P.ice ? `${ICE_NAME[P.ice]} + ${ICE_NAME[P.ice2]}` : ICE_NAME[P.ice]; }

const PRESETS = {
  un271: {
    orbit: { q: "10.94525581", e: "0.99869477", i: "95.55335892", Om: "190.0846361196", w: "326.33161880945", T: fmtDT(msOf(2462892.9068)) },
    ice: "CO2", cz: 1, HN: 8.40, beta: 0.04, lk: Math.log10(4.110e6), j: 1, ap: 6, msun: -26.93, y0: 2010, y1: 2036, ref: "2024-02-15", view: "mag", cut: "2024-06-01",
    obs: "# r′ magnitude in a 6″ radius aperture\n2021-12-15 18.99 0.05  # LOOK, Oct 2021–Feb 2022 (Kelley et al. 2022)\n2024-02-15 17.38 0.03  # LCO, Feb 2024 (Roth et al. 2025)",
    nbox: { r0: 16.70, r1: 19.4, n: 9.6, s: 0.3 },        // the coma index measured between the two points (Sec. VI)
    note: "Nucleus 120–140 km across (H<sub>N</sub> = 8.40 in r′). K is anchored on the 2024 point with CO₂ at the subsolar point (μ = 271), as in Fig. 4c of the paper. The later points are monthly medians of amateur CCD photometry from COBS (2024 Oct – 2026 Nov), not fitted: each magnitude, measured in a 8–73″ aperture, is converted to r′ in 6″ by removing the nucleus, scaling the coma flux by 6″/ρ (a steady coma) and taking V − r′ = 0.22. Converted this way they follow the CO₂ forecast (median offset 0.00 mag), while the empirical law through the two professional points, with a constant n, is about 0.45 mag too bright. The conversion is uncertain by about half a magnitude: the main observer’s apertures suggest a coma flux growing as ρ<sup>0.6</sup>, which would make the points 0.65 mag brighter and favor the constant n. Clear “fit only before” to fit all points."
  },
  t5: {
    orbit: { q: "3.84251371", e: "1.00240368", i: "52.38681511", Om: "100.67985066075", w: "352.41723883275", T: "2027-05-06 02:38" },
    ice: "CO2", cz: 1, HN: 15, beta: 0.04, lk: 4.997, j: 0, ap: 20000, msun: -26.76, y0: 2025, y1: 2029, ref: "2026-11-10", view: "afr", obs: "",
    note: "Perihelion 3.84 au on 2027 May 6, at the outer edge of the water zone. H<sub>N</sub> is a placeholder; K is fitted to the COBS total magnitudes (aperture 20 000 km). So far CO₂ or CO fits (χ²<sub>ν</sub> ≈ 2), water does not: Compare ices."
  },
  y1: {
    orbit: { q: "2.601924069", e: "0.999552254", i: "14.55698337", Om: "166.5518278289", w: "124.2978493120", T: fmtDT(msOf(2461692.5495)) },
    ice: "CO2", cz: 1, HN: 15, beta: 0.04, lk: 4.289, j: 0, ap: 20000, msun: -26.76, y0: 2025, y1: 2029, ref: "2026-12-19", view: "afr", obs: "",
    note: "Perihelion 2.6 au on 2027 Oct 14, where water switches on (JPL elements, 2026 Sep solution). H<sub>N</sub> is a placeholder; K is fitted to the COBS total magnitudes. So far CO₂ fits (χ²<sub>ν</sub> ≈ 1) and water does not; a brightening well above the CO₂ curve in 2027 would show water ice at work."
  }
};

const MORE = (orbit, o) => ({ orbit, beta: 0.04, j: 0, ap: 20000, msun: -26.76, cz: 1, lk: 4, view: "afr", obs: "", ...o });
Object.assign(PRESETS, {
  halley: MORE({ q: "0.57471580", e: "0.96794279", i: "162.18787119", Om: "59.114482936739", w: "112.25779200487", T: "2061-07-28 12:00" },
    { ice: "H2O", HN: 13.91, y0: 2050, y1: 2066, ref: "2060-10-28",
      note: "Now near aphelion at 35 au; next perihelion 2061 July 28 at about 0.6 au. Elements of the 1986 apparition (JPL) moved to the predicted 2061 perihelion, without planetary perturbations. H<sub>N</sub> from Giotto’s nucleus, 15 × 8 × 8 km (D ≈ 11 km seen side-on, p = 0.04); K is a placeholder." }),
  encke: MORE({ q: "0.33794828", e: "0.84774970", i: "11.41227811", Om: "334.19358460368", w: "187.13424636959", T: fmtDT(msOf(2461447.5650)) },
    { ice: "H2O", HN: 15.53, y0: 2026, y1: 2028, ref: "2026-01-01",
      note: "Period 3.3 years, perihelion 2027 Feb 11 at 0.34 au. Only a small part of its surface is active: try the dust mantle. H<sub>N</sub> from R<sub>N</sub> ≈ 2.4 km and p = 0.047; K is a placeholder." }),
  sw1: MORE({ q: "5.77512802", e: "0.04479814", i: "9.36477141", Om: "312.39091654503", w: "50.662700071797", T: fmtDT(msOf(2458589.9259)) },
    { lk: 6.080, ice: "CO", HN: 10.22, y0: 2010, y1: 2040, ref: "2026-01-01",
      note: "Nearly circular orbit between 5.8 and 6.3 au, always active and known for outbursts. Compare CO, which works at any distance, with CO₂ and water. H<sub>N</sub> from R<sub>N</sub> ≈ 30 km and p = 0.04; K is fitted to the COBS magnitudes 2010–2026, above which the outbursts stand out." }),
  p67: MORE({ q: "1.210", e: "0.64989", i: "3.8719", Om: "36.33", w: "22.15", T: "2028-04-09 12:00" },
    { ice: "H2O", HN: 15.99, y0: 2026, y1: 2031, ref: "2027-01-01",
      note: "Rosetta’s comet; next perihelion 2028 Apr 9 at 1.21 au. H<sub>N</sub> from the mean radius 1.65 km and p<sub>V</sub> = 0.065; K is a placeholder." }),
  p12: MORE({ q: "0.78086475", e: "0.95456023", i: "74.19083929", Om: "255.85505059528", w: "198.98800614657", T: fmtDT(msOf(2460421.6316)) },
    { lk: 5.124, ice: "H2O", HN: 13, y0: 2022, y1: 2026, ref: "2023-07-09", cut: "2024-01-01",
      note: "Halley-type comet with a 71-year period and repeated outbursts before its perihelion of 2024 Apr 21 at 0.78 au. As in Fig. 4b of the paper, both laws are fitted to the rise up to the end of 2023 (2.1 au): the empirical law gives <i>n</i> = 8.8 and misses the rest by 5.0 mag, water ice by 0.5 mag. The outbursts of 2023 stand above the water curve. Its nucleus is not measured: H<sub>N</sub> is a placeholder; K is fitted to the COBS magnitudes before the cut." }),
  a3: MORE({ q: "0.39145906", e: "1.00012475", i: "139.11502467", Om: "21.558609111275", w: "308.4872848828", T: fmtDT(msOf(2460581.2353)) },
    { lk: 4.571, ice: "CO2", HN: 15, y0: 2023, y1: 2026, ref: "2024-01-22",
      note: "The bright comet of October 2024; perihelion 2024 Sep 27 at 0.39 au on a nearly parabolic orbit. H<sub>N</sub> is a placeholder; K is fitted to the COBS magnitudes. CO₂ fits the whole arc far better than water (χ²<sub>ν</sub> 19 against 120): it was active at 7 au." }),
  neo: MORE({ q: "0.295", e: "0.99921", i: "128.94", Om: "61.010", w: "37.279", T: "2020-07-03 12:00" },
    { lk: 4.379, ice: "H2O", HN: 15, y0: 2019, y1: 2022, ref: "2019-11-03",
      note: "The bright comet of July 2020; perihelion 2020 Jul 3 at 0.29 au. Its peak water production needs a nucleus of at least 3.5 km radius: see the power budget below. H<sub>N</sub> is a placeholder; K is fitted to the COBS magnitudes. The faint points in April 2020 are small-aperture CCD magnitudes." }),
  hb: MORE({ q: "0.9141", e: "0.9951", i: "89.43", Om: "282.47", w: "130.59", T: "1997-04-01 03:19" },
    { lk: 7.324, split: true, post: { lk: 7.2716, k: 1 }, ice: "CO", HN: 9.77, y0: 1995, y1: 2000, ref: "1995-07-23",
      note: "The Great Comet of 1997, discovered on 1995 Jul 23 at 7.1 au, far beyond the water zone; perihelion 1997 Apr 1 at 0.91 au. H<sub>N</sub> for a nucleus 37 km in radius (Szabó et al. 2012; earlier estimates 30 ± 10 km) and p = 0.04. As in Fig. S1b of the paper's supplementary material, CO and the empirical law are fitted separately before and after perihelion: rms 0.31 and 0.24 mag, <i>n</i> = 3.1 and 3.5, and the active area after perihelion is 0.89 times that before. Water alone does not fit (Compare ices). For the forecast of Fig. 4a, untick “before and after perihelion”, enter 1996-03-15 under “fit only before” and Fit: the empirical law (<i>n</i> = 4.6) misses the rest by 1.2 mag, CO by 0.4 mag." }),
  hyak: MORE({ q: "0.23022353", e: "0.99989165", i: "124.92204939", Om: "188.04513199216", w: "130.1751209781", T: fmtDT(msOf(2450204.8941)) },
    { lk: 4.804, ice: "H2O", HN: 15.71, y0: 1996, y1: 1997, ref: "1996-01-30", view: "mag", cursor: "1996-03-25",
      note: "The Great Comet of 1996. Discovered on 1996 Jan 30, it passed 0.10 au from the Earth on Mar 25, five weeks before perihelion at 0.23 au: the magnitude view shows the close approach. H<sub>N</sub> from the nucleus seen in the infrared (Lisse et al. 1999: D ≈ 4.8 km, with p = 0.04; radar echoes suggested 1–3 km); K is fitted to the COBS magnitudes. Its SWAN water exponents are the points at q = 0.23 au below." }),
  i3: MORE({ q: "1.35638454", e: "6.13941774", i: "175.11310480", Om: "322.15684249105", w: "128.01051366675", T: fmtDT(msOf(2460977.9814)) },
    { lk: 3.999, ice: "H2O", HN: 17.04, y0: 2025, y1: 2027, ref: "2025-07-16",
      note: "The third interstellar object, on a hyperbola with e = 6.14; perihelion 2025 Oct 29 at 1.36 au. H<sub>N</sub> from R<sub>N</sub> = 1.3 km and p = 0.04; K is fitted to the COBS magnitudes. Its water production exceeds what this nucleus can deliver: see the power budget below." }),
});

// comets of our forecast test of 24 light curves (supplementary material): ice and K from the first half, cut at the split
Object.assign(PRESETS, {
  q2: MORE({ q: "1.290358238881392", e: "0.9977738520356274", i: "80.30130181438155", Om: "94.97505909428631", w: "12.39526184803142", T: fmtDT(msOf(2457052.569249606)) },
    { lk: 5.188, ice: "H2O", cz: 0.25, HN: 15, y0: 2014, y1: 2016, ref: "2014-10-03", cut: "2015-04-08",
      note: "The bright green comet of January 2015; perihelion 2015 Jan 30 at 1.29 au. In our forecast test of 24 comets, water on an isothermal nucleus fits the points before the cut; the forecast rms is 4.19 mag, against 5.51 for the empirical law with fitted <i>n</i> (<i>n</i> = 10.4) and 1.36 for <i>n</i> = 3.5: the chosen ice switches on within the fitted half. H<sub>N</sub> is a placeholder." }),
  k2: MORE({ q: "1.797281961376049", e: "1.000616562373503", i: "87.55487509681637", Om: "88.24109127986824", w: "236.1923489989567", T: fmtDT(msOf(2459933.227066954)) },
    { lk: 5.956, ice: "CO", cz: 1, HN: 15, y0: 2018, y1: 2026, ref: "2018-03-22", cut: "2021-12-04",
      note: "Active at 16 au; perihelion 2022 Dec 19 at 1.80 au. In our forecast test of 24 comets, CO at the subsolar point fits the points before the cut; the forecast rms is 1.65 mag, against 1.22 for the empirical law with fitted <i>n</i> (<i>n</i> = 2.9) and 1.65 for <i>n</i> = 3.5: the chosen ice does not switch on within the fitted half, so the model is the law with <i>n</i> ≈ 3.5. H<sub>N</sub> is a placeholder." }),
  t2: MORE({ q: "1.615095955888523", e: "0.9996765980850288", i: "57.23142385337199", Om: "64.37877215685057", w: "92.99331253710214", T: fmtDT(msOf(2458974.449996825)) },
    { lk: 4.649, ice: "CO", cz: 1, HN: 15, y0: 2018, y1: 2021, ref: "2018-11-12", cut: "2019-12-30",
      note: "Perihelion 2020 May 4 at 1.62 au. In our forecast test of 24 comets, CO at the subsolar point fits the points before the cut; the forecast rms is 0.38 mag, against 0.42 for the empirical law with fitted <i>n</i> (<i>n</i> = 3.2) and 0.38 for <i>n</i> = 3.5: the chosen ice does not switch on within the fitted half, so the model is the law with <i>n</i> ≈ 3.5. H<sub>N</sub> is a placeholder." }),
  p1: MORE({ q: "1.550535584223489", e: "1.000997435576857", i: "106.1772916768743", Om: "325.9972285528648", w: "90.74706541748093", T: fmtDT(msOf(2455919.177555371)) },
    { lk: 5.554, ice: "CO", cz: 1, HN: 15, y0: 2010, y1: 2013, ref: "2010-11-08", cut: "2011-10-16",
      note: "C/2009 P1 (Garradd), perihelion 2011 Dec 23 at 1.55 au and bright for months afterwards. In our forecast test of 24 comets, CO at the subsolar point fits the points before the cut; the forecast rms is 0.50 mag, against 0.20 for the empirical law with fitted <i>n</i> (<i>n</i> = 2.6) and 0.50 for <i>n</i> = 3.5: the chosen ice does not switch on within the fitted half, so the model is the law with <i>n</i> ≈ 3.5. H<sub>N</sub> is a placeholder." }),
  e3: MORE({ q: "1.11224437022534", e: "1.000301905819192", i: "109.169480756749", Om: "302.5550197168474", w: "145.81492879", T: fmtDT(msOf(2459957.28519883)) },
    { lk: 4.062, ice: "H2O", cz: 1, HN: 15, y0: 2022, y1: 2024, ref: "2022-04-13", cut: "2022-10-22",
      note: "The green comet of early 2023; perihelion 2023 Jan 12 at 1.11 au. In our forecast test of 24 comets, water at the subsolar point fits the points before the cut; the forecast rms is 0.30 mag, against 0.48 for the empirical law with fitted <i>n</i> (<i>n</i> = 4.6) and 0.60 for <i>n</i> = 3.5: the chosen ice switches on within the fitted half. H<sub>N</sub> is a placeholder." }),
  l3: MORE({ q: "3.554492702308672", e: "1.001405275505796", i: "48.36374698279935", Om: "290.7888867057587", w: "171.6137051016775", T: fmtDT(msOf(2459589.132623246)) },
    { lk: 6.035, ice: "CO", ice2: "H2O", k: 0.04, split: true, post: { lk: 6.0265, k: 0.897 }, cz: 1, HN: 15, y0: 2020, y1: 2024, ref: "2020-01-06",
      note: "Active from 8 au; perihelion 2022 Jan 9 at 3.55 au. As in Fig. S1c of the paper's supplementary material, CO and water ice on separate patches are fitted before perihelion (rms 0.26 mag, k = 0.04). No single ice fits the steep rise (<i>n</i> = 7.8): with 2022-01-09 under “fit only before”, Compare ices gives χ²<sub>ν</sub> = 2.7 for the pair and 22 or more for any single ice. Far out, CO supplied the dust; as water switched on, its share grew and the index rose. After perihelion the comet faded slowly (<i>n</i> = 3.6) and the fit gives k = 0.90, as if its water patches had stopped working. In our forecast test of 24 comets (cut 2021-10-04), CO₂ on an isothermal nucleus fits the points before the cut; the forecast rms is 1.77 mag, against 1.53 for the empirical law with fitted <i>n</i> (<i>n</i> = 7.5) and 1.86 for <i>n</i> = 3.5. H<sub>N</sub> is a placeholder." }),
  l4: MORE({ q: "0.3015442238164768", e: "1.000032767108972", i: "84.20819702923397", Om: "65.6658863915034", w: "333.6516376633611", T: fmtDT(msOf(2456361.669960002)) },
    { lk: 5.056, ice: "CO", cz: 1, HN: 15, y0: 2012, y1: 2014, ref: "2012-05-22", cut: "2013-03-05",
      note: "C/2011 L4 (PANSTARRS), a naked-eye comet in March 2013; perihelion 2013 Mar 10 at 0.30 au. In our forecast test of 24 comets, CO at the subsolar point fits the points before the cut; the forecast rms is 0.96 mag, against 0.92 for the empirical law with fitted <i>n</i> (<i>n</i> = 2.9) and 0.96 for <i>n</i> = 3.5: the chosen ice does not switch on within the fitted half, so the model is the law with <i>n</i> ≈ 3.5. H<sub>N</sub> is a placeholder." }),
  a6: MORE({ q: "0.5299041302912412", e: "0.9955362293967411", i: "143.663707214442", Om: "108.0979873301086", w: "132.9688404751593", T: fmtDT(msOf(2460988.039387102)) },
    { lk: 4.196, ice: "CO", cz: 1, HN: 15, y0: 2025, y1: 2027, ref: "2025-08-26", cut: "2026-01-17",
      note: "C/2025 A6 (Lemmon), bright in October 2025; perihelion 2025 Nov 8 at 0.53 au. In our forecast test of 24 comets, CO at the subsolar point fits the points before the cut; the forecast rms is 1.37 mag, against 1.47 for the empirical law with fitted <i>n</i> (<i>n</i> = 3.3) and 1.37 for <i>n</i> = 3.5: the chosen ice does not switch on within the fitted half, so the model is the law with <i>n</i> ≈ 3.5. H<sub>N</sub> is a placeholder." }),
  v2j: MORE({ q: "1.636950053071025", e: "1.001702614530096", i: "49.87481167153636", Om: "69.8521003063327", w: "164.8980978916873", T: fmtDT(msOf(2457916.847458616)) },
    { lk: 4.625, ice: "CO", cz: 1, HN: 15, y0: 2016, y1: 2019, ref: "2016-02-10", cut: "2017-01-26",
      note: "C/2015 V2 (Johnson), perihelion 2017 Jun 12 at 1.64 au. In our forecast test of 24 comets, CO at the subsolar point fits the points before the cut; the forecast rms is 0.62 mag, against 0.24 for the empirical law with fitted <i>n</i> (<i>n</i> = 2.5) and 0.62 for <i>n</i> = 3.5: the chosen ice does not switch on within the fitted half, so the model is the law with <i>n</i> ≈ 3.5. H<sub>N</sub> is a placeholder." }),
  v2z: MORE({ q: "2.227959873700067", e: "1.00089953037377", i: "131.6110418896687", Om: "212.3712351547342", w: "162.4243172377264", T: fmtDT(msOf(2460073.056589952)) },
    { lk: 5.168, ice: "CO2", cz: 0.25, HN: 15, y0: 2021, y1: 2026, ref: "2021-03-09", cut: "2023-03-28",
      note: "C/2020 V2 (ZTF), active beyond 8 au; perihelion 2023 May 8 at 2.23 au. In our forecast test of 24 comets, CO₂ on an isothermal nucleus fits the points before the cut; the forecast rms is 0.76 mag, against 0.64 for the empirical law with fitted <i>n</i> (<i>n</i> = 4.1) and 0.46 for <i>n</i> = 3.5: the chosen ice switches on within the fitted half. H<sub>N</sub> is a placeholder." }),
  us10: MORE({ q: "0.8229694514684749", e: "1.000280139519764", i: "148.8782967756639", Om: "186.1445797203693", w: "340.358754276709", T: fmtDT(msOf(2457342.221641829)) },
    { lk: 5.003, ice: "CO", cz: 1, HN: 15, y0: 2014, y1: 2017, ref: "2014-06-19", cut: "2015-08-22",
      note: "C/2013 US10 (Catalina), first seen by amateurs at 8.2 au; perihelion 2015 Nov 15 at 0.82 au. In our forecast test of 24 comets, CO at the subsolar point fits the points before the cut; the forecast rms is 1.21 mag, against 1.11 for the empirical law with fitted <i>n</i> (<i>n</i> = 3.1) and 1.21 for <i>n</i> = 3.5: the chosen ice does not switch on within the fitted half, so the model is the law with <i>n</i> ≈ 3.5. H<sub>N</sub> is a placeholder." }),
  k1: MORE({ q: "1.054597098294", e: "1.000152915493971", i: "142.4282581726916", Om: "317.7376491689354", w: "203.1066008825931", T: fmtDT(msOf(2456897.160363006)) },
    { lk: 4.816, ice: "CO", cz: 1, HN: 15, y0: 2013, y1: 2016, ref: "2013-03-24", cut: "2014-04-28",
      note: "C/2012 K1 (PANSTARRS), first seen by amateurs at 8.7 au; perihelion 2014 Aug 27 at 1.05 au. In our forecast test of 24 comets, CO at the subsolar point fits the points before the cut; the forecast rms is 0.58 mag, against 0.45 for the empirical law with fitted <i>n</i> (<i>n</i> = 3.0) and 0.58 for <i>n</i> = 3.5: the chosen ice does not switch on within the fitted half, so the model is the law with <i>n</i> ≈ 3.5. H<sub>N</sub> is a placeholder." }),
  g2: MORE({ q: "4.982360288479666", e: "1.000603552652205", i: "48.46954508656729", Om: "221.095294828915", w: "343.2672327577748", T: fmtDT(msOf(2460562.634007078)) },
    { lk: 5.195, ice: "CO2", cz: 1, HN: 15, y0: 2023, y1: 2027, ref: "2023-04-13", cut: "2025-05-20",
      note: "C/2021 G2 (ATLAS), a distant comet that came no closer than 4.98 au (perihelion 2024 Sep 9). In our forecast test of 24 comets, CO₂ at the subsolar point fits the points before the cut; the forecast rms is 0.33 mag, against 0.37 for the empirical law with fitted <i>n</i> (<i>n</i> = 4.0) and 0.28 for <i>n</i> = 3.5: the chosen ice does not switch on within the fitted half, so the model is the law with <i>n</i> ≈ 3.5. H<sub>N</sub> is a placeholder." }),
  e2: MORE({ q: "3.666295808424358", e: "1.001060418168188", i: "137.1316685887328", Om: "125.3761151010601", w: "41.72559864492467", T: fmtDT(msOf(2460567.63388719)) },
    { lk: 4.762, ice: "CO", cz: 1, HN: 15, y0: 2023, y1: 2027, ref: "2023-01-18", cut: "2025-02-01",
      note: "C/2022 E2 (ATLAS), first seen by amateurs at 7.8 au; perihelion 2024 Sep 14 at 3.67 au. In our forecast test of 24 comets, CO at the subsolar point fits the points before the cut; the forecast rms is 0.86 mag, against 0.92 for the empirical law with fitted <i>n</i> (<i>n</i> = 3.2) and 0.86 for <i>n</i> = 3.5: the chosen ice does not switch on within the fitted half, so the model is the law with <i>n</i> ≈ 3.5. H<sub>N</sub> is a placeholder." }),
  n2: MORE({ q: "3.825529163602361", e: "1.003236448782429", i: "5.503022530288324", Om: "319.7406957022052", w: "75.38965400185295", T: fmtDT(msOf(2460888.270675301)) },
    { lk: 4.354, ice: "CO2", cz: 1, HN: 15, y0: 2023, y1: 2027, ref: "2023-10-17", cut: "2025-08-08",
      note: "C/2022 N2 (PANSTARRS), first seen by amateurs at 7.0 au; perihelion 2025 Jul 31 at 3.83 au. In our forecast test of 24 comets, CO₂ at the subsolar point fits the points before the cut; the forecast rms is 0.32 mag, against 0.33 for the empirical law with fitted <i>n</i> (<i>n</i> = 3.9) and 0.32 for <i>n</i> = 3.5: the chosen ice does not switch on within the fitted half, so the model is the law with <i>n</i> ≈ 3.5. H<sub>N</sub> is a placeholder." }),
});

// two famous disappointments, run through the rules of the forecast test (not among its 24 comets)
Object.assign(PRESETS, {
  koh: MORE({ q: "0.1424250322974156", e: "1.000007809885518", i: "14.3041449344601", Om: "258.4894909953879", w: "37.79756893411028", T: fmtDT(msOf(2442044.930683474144)) },
    { lk: 4.8026, ice: "CO", cz: 1, HN: 15, y0: 1973, y1: 1975, ref: "1973-10-06", cut: "1974-01-08",
      note: "The “comet of the century” of 1973 that disappointed. Discovered on 1973 Mar 18 at 4.6 au; extrapolations with a fixed index (<i>n</i> = 4 to 6) promised magnitude −10 or brighter, but it peaked near −3 at perihelion on 1973 Dec 28, 0.14 au from the Sun and too close to it for ground observers. COBS holds almost no magnitudes from before October 1973. With the rules of our forecast test (it is not one of the 24 comets), CO at the subsolar point fits the points before the cut; the forecast rms is 1.55 mag, against 1.46 for the empirical law with fitted <i>n</i> (<i>n</i> = 2.8) and 1.55 for <i>n</i> = 3.5. No law foresees how fast it faded after perihelion: with one μ before and one after, μ<sub>post</sub>/μ<sub>pre</sub> = 0.27. H<sub>N</sub> is a placeholder." }),
  ison: MORE({ q: "0.01246672509220626", e: "1.00000509491957", i: "62.16289620932959", Om: "295.6865414817095", w: "345.5411966905318", T: fmtDT(msOf(2456625.264531166622)) },
    { lk: 4.1954, ice: "CO", cz: 1, HN: 15, y0: 2012, y1: 2014, ref: "2012-12-17", cut: "2013-05-06",
      note: "The sungrazer of 2013, which prompted Parodi (2013), the forerunner of this model. Discovered on 2012 Sep 21 at 6.3 au; early extrapolations promised a comet brighter than the full Moon. Its brightness stalled in the first half of 2013 and rose again after its conjunction with the Sun; it reached about magnitude −2 on Nov 27 and disintegrated at perihelion on 2013 Nov 28, 0.012 au from the Sun’s centre. The last COBS bin is Nov 18. With the rules of our forecast test (it is not one of the 24 comets), CO at the subsolar point fits the points before the cut; the forecast rms is 1.29 mag, against 3.67 for the empirical law with fitted <i>n</i> (<i>n</i> = 0.1, the stall) and 1.30 for <i>n</i> = 3.5. The fitted index, frozen in the stall, misses the rise by magnitudes; the physical law does not. It forecasts about magnitude −0.5 for Nov 27 and −10 for perihelion: the comet broke up before. H<sub>N</sub> is a placeholder (the nucleus was less than 1 km across)." }),
});

// COBS light curves (data/cobs_lightcurves.txt, CC BY-NC-SA 4.0), appended to the presets
const COBS_HEAD = {
  un271: "# COBS CCD photometry moved to r′ in 6″ (monthly medians, see the note; open circles in Fig. 4c of the paper), not fitted",
  other: "# total magnitudes from COBS (median per bin; CC BY-NC-SA 4.0)" };
if (typeof COBS === "object") for (const [k, txt] of Object.entries(COBS)) {
  const p = PRESETS[k]; if (!p) continue;
  p.obs = (p.obs ? p.obs + "\n" : "") + (COBS_HEAD[k] || COBS_HEAD.other) + "\n" + txt;
  if (k !== "un271") p.view = "mag";
}

const S = { view: "mag", cursor: Date.now(), cmp: null, nbox: null, series: null, post: null };

// ---------- small helpers
function num(s) { const v = parseFloat(String(s).replace("−", "-").replace(",", ".")); return Number.isFinite(v) ? v : NaN; }
function pad(n) { return String(n).padStart(2, "0"); }
function fmtD(ms) { const d = new Date(ms); return `${d.getUTCFullYear()}-${pad(d.getUTCMonth() + 1)}-${pad(d.getUTCDate())}`; }
function fmtDT(ms) { const d = new Date(ms); return `${fmtD(ms)} ${pad(d.getUTCHours())}:${pad(d.getUTCMinutes())}`; }
function parseDate(s) {
  s = String(s).trim(); if (!s) return NaN;
  let m = s.match(/^(\d{4})-(\d{1,2})-(\d{1,2})(?:[ T](\d{1,2}):(\d{2})(?::(\d{2}))?)?$/);
  if (m) return Date.UTC(+m[1], +m[2] - 1, +m[3], +(m[4] || 0), +(m[5] || 0), +(m[6] || 0));
  const v = num(s); if (!Number.isFinite(v)) return NaN;
  if (v > 2e6) return msOf(v);                                   // Julian date
  if (v > 1000 && v < 3000) { const y = Math.floor(v), a = Date.UTC(y, 0, 1), b = Date.UTC(y + 1, 0, 1); return a + (v - y) * (b - a); }
  return NaN;
}
function yearOf(ms) { const y = new Date(ms).getUTCFullYear(), a = Date.UTC(y, 0, 1), b = Date.UTC(y + 1, 0, 1); return y + (ms - a) / (b - a); }
function msOfYear(v) { const y = Math.floor(v), a = Date.UTC(y, 0, 1), b = Date.UTC(y + 1, 0, 1); return a + (v - y) * (b - a); }
function mi(s) { return String(s).replace(/-/g, "−"); }
function sci(v, d = 3) {
  if (!Number.isFinite(v)) return "—"; if (v === 0) return "0";
  const e = Math.floor(Math.log10(Math.abs(v)));
  if (e >= d - 1 && e <= 4) return mi(Number(v.toPrecision(d)).toLocaleString("en-US"));
  if (e >= -2 && e < d - 1) return mi(Number(v.toPrecision(d)).toString());
  return `${mi((v / 10 ** e).toFixed(d - 1))} × 10<sup>${mi(e)}</sup>`;
}
function css(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
function alpha(color, a) {
  const c = color.replace("#", ""); if (c.length !== 6) return color;
  return `rgba(${parseInt(c.slice(0, 2), 16)},${parseInt(c.slice(2, 4), 16)},${parseInt(c.slice(4, 6), 16)},${a})`;
}

// ---------- reading the controls
function readOrbit() {
  const q = num($("oq").value), e = num($("oe").value), i = num($("oi").value), Om = num($("oO").value), w = num($("ow").value);
  const T = parseDate($("oT").value);
  const ok = q > 0 && e >= 0 && Number.isFinite(T);
  return { ok, q, e, i, Om, w, Tp: jdOf(T), Tms: T, full: ok && [i, Om, w].every(Number.isFinite) };
}
function params() {
  const i2 = document.querySelector('input[name="ice2"]:checked');
  return { ice: document.querySelector('input[name="ice"]:checked').value, ice2: i2 ? i2.value : "none", k: +$("kmix").value, split: $("split").checked,
    cz: +$("cz").value, HN: +$("hn").value, beta: +$("beta").value,
    lk: +$("lk").value, h: +$("hm").value >= 2.995 ? Infinity : 10 ** +$("hm").value, j: +$("j").value, ap: Math.max(1e-3, num($("ap").value) || 6), msun: Number.isFinite(num($("msun").value)) ? num($("msun").value) : -26.76 };
}
function windowMs() {
  let y0 = num($("y0").value), y1 = num($("y1").value);
  if (!Number.isFinite(y0)) y0 = 2020; if (!Number.isFinite(y1) || y1 <= y0) y1 = y0 + 1;
  return [Date.UTC(Math.round(y0), 0, 1), Date.UTC(Math.round(y1), 0, 1)];
}
function cutMs() { const v = parseDate($("cut").value); return Number.isFinite(v) ? v : Infinity; }
// empirical law m = m0 + 5 lg D + 2.5 n lg r, weighted least squares in (m0, n)
function fitEmp(pts) {
  let S0 = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
  for (const p of pts) { const w = 1 / p.s ** 2, x = 2.5 * Math.log10(p.r), y = p.m - 5 * Math.log10(p.D); S0 += w; Sx += w * x; Sy += w * y; Sxx += w * x * x; Sxy += w * x * y; }
  const det = S0 * Sxx - Sx * Sx; if (pts.length < 2 || !(Math.abs(det) > 1e-12 * S0 * S0)) return null;
  const n = (S0 * Sxy - Sx * Sy) / det; return { n, m0: (Sy - n * Sx) / S0 };
}
function empMag(f, r, D) { return f.m0 + 5 * Math.log10(D) + 2.5 * f.n * Math.log10(r); }
function parseObs() {
  const out = [], bad = [];
  $("obs").value.split(/\n/).forEach((line, k) => {
    const t = line.replace(/#.*/, "").trim(); if (!t) return;
    const p = t.split(/[\s,;]+/); let ms, rest;
    if (p.length >= 3 && /^\d{1,2}:\d{2}/.test(p[1])) { ms = parseDate(p[0] + " " + p[1]); rest = p.slice(2); } else { ms = parseDate(p[0]); rest = p.slice(1); }
    const m = num(rest[0]), s = rest.length > 1 ? num(rest[1]) : 0.1;
    if (Number.isFinite(ms) && Number.isFinite(m)) out.push({ ms, m, s: s > 0 ? s : 0.1 }); else bad.push(k + 1);
  });
  return { obs: out, bad };
}

// ---------- model series over the window
// the model piece in use at time ms: the sliders, or after perihelion the separate fit (S.post)
function pieceAt(el, P, ms) {
  if (P.split && S.post && ms >= el.Tms) { const K = 10 ** S.post.lk; return { mix: mixFrom(P.ice, P.ice2, S.post.k), K, mu: muFrom(K, P), post: true }; }
  const K = 10 ** P.lk; return { mix: mixOf(P), K, mu: muFrom(K, P), post: false };
}
function modelAt(el, P, g, ms) { const pc = pieceAt(el, P, ms), et = etaMix(g.r, pc.mix, P.cz, P.h); return { ...magnitude(g.r, g.D, g.alpha, et, P, pc.mu), et, pc }; }
function compute() {
  const el = readOrbit(), P = params(), [t0, t1] = windowMs();
  const N = 520, t = [], r = [], D = [], al = [], m = [], mN = [], A = [];
  const K = 10 ** P.lk, mu = muFrom(K, P);
  if (el.ok) for (let k = 0; k <= N; k++) {
    const ms = t0 + (t1 - t0) * k / N, g = ephem(el, ms), pc = pieceAt(el, P, ms), et = etaMix(g.r, pc.mix, P.cz, P.h);
    t.push(ms); r.push(g.r); D.push(g.D); al.push(g.alpha); A.push(afrho(pc.K, g.r, et));
    if (el.full) { const mg = magnitude(g.r, g.D, g.alpha, et, P, pc.mu); m.push(mg.m); mN.push(mg.mN); } else { m.push(NaN); mN.push(NaN); }
  }
  let emp = null, empPost = null, mE = [];
  if (el.full && $("emp").checked) {
    const cut = cutMs(), pts = parseObs().obs.filter((o) => o.ms < cut).map((o) => { const g = ephem(el, o.ms); return { ...o, r: g.r, D: g.D }; });
    if (P.split) { emp = fitEmp(pts.filter((o) => o.ms < el.Tms)); empPost = fitEmp(pts.filter((o) => o.ms >= el.Tms)); }
    else emp = fitEmp(pts);
    if (emp || empPost) mE = r.map((rv, k) => { const f = t[k] >= el.Tms && empPost ? empPost : emp || empPost; return empMag(f, rv, D[k]); });
  }
  S.series = { el, P, t0, t1, t, r, D, al, m, mN, A, K, mu, emp, empPost, mE };
  return S.series;
}
function empAt(s, ms) { return ms >= s.el.Tms && s.empPost ? s.empPost : s.emp || s.empPost; }

// ---------- canvas plumbing
function frame(cv, ar, minH, maxH) {
  const dpr = Math.min(2, window.devicePixelRatio || 1), w = Math.max(240, Math.floor(cv.clientWidth || cv.parentNode.clientWidth - 24));
  const h = Math.max(minH, Math.min(maxH || 1e9, Math.round(w * ar)));
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
  cv.style.height = h + "px";
  const g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, w, h);
  return { g, w, h };
}
function lin(d0, d1, p0, p1) { const f = (v) => p0 + (v - d0) / (d1 - d0) * (p1 - p0); f.inv = (p) => d0 + (p - p0) / (p1 - p0) * (d1 - d0); f.d = [d0, d1]; f.p = [p0, p1]; return f; }
function logs(d0, d1, p0, p1) { const a = Math.log10(d0), b = Math.log10(d1); const f = (v) => p0 + (Math.log10(v) - a) / (b - a) * (p1 - p0); f.inv = (p) => 10 ** (a + (p - p0) / (p1 - p0) * (b - a)); f.d = [d0, d1]; f.p = [p0, p1]; f.log = true; return f; }
function niceStep(span, n) { const s0 = span / n, mg = 10 ** Math.floor(Math.log10(s0)), e = s0 / mg; return (e >= 7.5 ? 10 : e >= 3.5 ? 5 : e >= 1.5 ? 2 : 1) * mg; }
function linTicks(d0, d1, n) { const lo = Math.min(d0, d1), hi = Math.max(d0, d1), st = niceStep(hi - lo, n), out = []; for (let v = Math.ceil(lo / st - 1e-9) * st; v <= hi + 1e-9 * st; v += st) out.push(+v.toFixed(10)); return out; }
function logTicks(d0, d1) {
  const out = [], a = Math.floor(Math.log10(d0)), b = Math.ceil(Math.log10(d1)), many = b - a > 4;
  for (let e = a; e <= b; e++) for (const f of many ? [1] : [1, 2, 5]) { const v = f * 10 ** e; if (v >= d0 * 0.999 && v <= d1 * 1.001) out.push(v); }
  return out;
}
function fmtTick(v) { if (v >= 1e4 || (v > 0 && v < 1e-2)) { const e = Math.round(Math.log10(v)); return Math.abs(v - 10 ** e) / v < 1e-6 ? `1e${e}` : v.toExponential(0); } return mi(String(+v.toPrecision(6))); }
const FONT = '12px "IBM Plex Mono", ui-monospace, Menlo, monospace', FONTS = '11px "IBM Plex Mono", ui-monospace, Menlo, monospace', FONTL = '13px "IBM Plex Sans Condensed", "Arial Narrow", system-ui, sans-serif';
function axes(g, X, Y, o) {
  const [l, r] = X.p, [b, t] = Y.p, ink = css("--mute"), rule = css("--rule");
  g.save(); g.lineWidth = 1; g.font = FONT; g.fillStyle = ink;
  for (const v of o.xt) { const x = X(v); g.strokeStyle = alpha(rule, 0.9); g.beginPath(); g.moveTo(x, t); g.lineTo(x, b); g.stroke(); g.textAlign = "center"; g.textBaseline = "top"; g.fillText(o.xf ? o.xf(v) : fmtTick(v), x, b + 5); }
  for (const v of o.yt) { const y = Y(v); g.strokeStyle = alpha(rule, 0.9); g.beginPath(); g.moveTo(l, y); g.lineTo(r, y); g.stroke(); g.textAlign = "right"; g.textBaseline = "middle"; g.fillText(o.yf ? o.yf(v) : fmtTick(v), l - 6, y); }
  g.strokeStyle = rule; g.strokeRect(l, t, r - l, b - t);
  g.font = FONTL; g.fillStyle = ink;
  if (o.xl) { g.textAlign = "center"; g.textBaseline = "bottom"; g.fillText(o.xl, (l + r) / 2, g.canvas.height / (Math.min(2, window.devicePixelRatio || 1)) - 2); }
  if (o.yl) { g.save(); g.translate(12, (t + b) / 2); g.rotate(-Math.PI / 2); g.textAlign = "center"; g.textBaseline = "middle"; g.fillText(o.yl, 0, 0); g.restore(); }
  g.restore();
}
function path(g, xs, ys, X, Y, style) {
  g.save(); g.strokeStyle = style.c; g.lineWidth = style.w || 1.5; g.setLineDash(style.dash || []); g.lineJoin = "round";
  g.beginPath(); let on = false;
  for (let k = 0; k < xs.length; k++) {
    const yv = ys[k]; if (!Number.isFinite(yv) || (Y.log && yv <= 0)) { on = false; continue; }
    const x = X(xs[k]), y = Y(yv); if (on) g.lineTo(x, y); else { g.moveTo(x, y); on = true; }
  }
  g.stroke(); g.restore();
}
function clip(g, X, Y) { g.save(); g.beginPath(); g.rect(X.p[0], Y.p[1], X.p[1] - X.p[0], Y.p[0] - Y.p[1]); g.clip(); }

// ---------- light curve
function drawLC() {
  const cv = $("c1"), { g, w, h } = frame(cv, 0.44, 280, 460), s = S.series, P = s.P;
  const box = { l: 52, r: w - 14, t: 16, b: h - 40 };
  const X = lin(s.t0, s.t1, box.l, box.r);
  const two = mixOf(P).length > 1, iceC = two ? css("--ink") : css(ICE_VAR[P.ice]), ink = css("--ink"), mute = css("--mute"), ghost = css("--ghost"), obsC = css("--obs"), acc = css("--accent");
  const span = (s.t1 - s.t0) / (365.25 * DAY), st = span > 30 ? 5 : span > 12 ? 2 : span > 6 ? 1 : 0.5;
  const xt = []; for (let y = Math.ceil(yearOf(s.t0) / st) * st; y <= yearOf(s.t1) + 1e-9; y += st) xt.push(msOfYear(y));
  const xf = (v) => st < 1 ? (Math.abs(yearOf(v) % 1) < 1e-6 ? String(Math.round(yearOf(v))) : fmtD(v).slice(0, 7)) : String(Math.round(yearOf(v)));
  let Y;
  const { obs } = parseObs();
  S.box1 = { X, box };
  if (!s.el.ok) { msg(g, w, h, "Enter q, e and T of the orbit to draw the light curve."); return; }
  if (S.view === "mag") {
    if (!s.el.full) { msg(g, w, h, "Magnitudes need the full orbit (i, Ω, ω). The Afρ growth view works with q, e and T alone."); return; }
    let lo = Infinity, hi = -Infinity;
    s.m.forEach((v) => { if (Number.isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } });
    obs.forEach((o) => { if (o.ms >= s.t0 && o.ms <= s.t1) { lo = Math.min(lo, o.m - o.s); hi = Math.max(hi, o.m + o.s); } });
    s.mN.forEach((v) => { if (Number.isFinite(v) && v < hi + 3) hi = Math.max(hi, v); });
    if (!Number.isFinite(lo)) { lo = 10; hi = 20; }
    const pad = Math.max(0.3, 0.06 * (hi - lo)); lo = Math.floor((lo - pad) * 2) / 2; hi = Math.ceil((hi + pad) * 2) / 2;
    Y = lin(hi, lo, box.b, box.t);
    axes(g, X, Y, { xt, xf, yt: linTicks(lo, hi, 6), yl: "magnitude" });
    clip(g, X, Y);
    path(g, s.t, s.mN, X, Y, { c: ghost, w: 1.5, dash: [5, 4] });
    if (s.emp || s.empPost) path(g, s.t, s.mE, X, Y, { c: ink, w: 1.7, dash: [7, 3, 2, 3] });
    path(g, s.t, s.m, X, Y, { c: iceC, w: 2.4 });
    const cut = cutMs(), test = [];
    g.fillStyle = obsC; g.strokeStyle = obsC; g.lineWidth = 1.2;
    for (const o of obs) {
      if (o.ms >= cut) test.push(o);
      if (o.ms < s.t0 || o.ms > s.t1) continue; const x = X(o.ms); g.beginPath(); g.moveTo(x, Y(o.m - o.s)); g.lineTo(x, Y(o.m + o.s)); g.stroke();
      g.beginPath(); if (o.ms < cut) { g.arc(x, Y(o.m), 3.6, 0, 2 * Math.PI); g.fill(); } else { g.fillStyle = css("--panel"); g.arc(x, Y(o.m), 3.4, 0, 2 * Math.PI); g.fill(); g.stroke(); g.fillStyle = obsC; }
    }
    g.restore();
    let tst = "";
    if (test.length && Number.isFinite(cut)) {
      const rms = (f) => Math.sqrt(test.reduce((a, o) => { const gg = ephem(s.el, o.ms), d = o.m - f(gg); return a + d * d; }, 0) / test.length);
      const rmsT = (f) => Math.sqrt(test.reduce((a, o) => { const gg = ephem(s.el, o.ms), d = o.m - f(gg, o.ms); return a + d * d; }, 0) / test.length);
      const rm = rmsT((gg, ms) => modelAt(s.el, P, gg, ms).m), hasE = s.emp || s.empPost, re = hasE ? rmsT((gg, ms) => empMag(empAt(s, ms), gg.r, gg.D)) : NaN;
      S.forecast = { n: test.length, model: rm, emp: re };
      tst = ` <span class="key">forecast rms: model ${rm.toFixed(2)}${hasE ? `, empirical ${re.toFixed(2)}` : ""} mag</span>`;
    } else S.forecast = null;
    const nTxt = s.emp && s.empPost ? `${mi(s.emp.n.toFixed(1))} | ${mi(s.empPost.n.toFixed(1))}` : (s.emp || s.empPost) ? mi((s.emp || s.empPost).n.toFixed(1)) : "";
    const splitTxt = P.split && S.post ? " (before | after perihelion)" : "";
    $("key1").innerHTML = `<span class="key"><span class="sw" style="border-top-color:${iceC}"></span>${mixName(P)} model${splitTxt}</span>` + (nTxt ? ` <span class="key"><span class="sw e"></span>empirical law, <i>n</i> = ${nTxt}</span>` : "") +
      ` <span class="key"><span class="sw d"></span>nucleus alone</span> <span class="key"><span class="dot"></span>observed</span>` + (test.length && Number.isFinite(cut) ? ` <span class="key"><span class="ring"></span>after the cut</span>` + tst : "");
  } else {
    const ref = refMs(s), r0 = ephem(s.el, ref).r;
    const series = SCEN.map((sc) => s.r.map((r) => eta(r, sc.x, sc.c, P.h) / eta(r0, sc.x, sc.c, P.h) * (r / r0) ** -1.5));
    const mx = mixOf(P), cur = s.r.map((r) => etaMix(r, mx, P.cz, P.h) / etaMix(r0, mx, P.cz, P.h) * (r / r0) ** -1.5);
    let lo = Infinity, hi = -Infinity;
    [cur, ...series].forEach((a) => a.forEach((v) => { if (v > 0 && Number.isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } }));
    lo = Math.max(lo, 1e-4); hi = Math.min(Math.max(hi, 2), 1e6);
    lo = 10 ** Math.floor(Math.log10(lo)); hi = 10 ** Math.ceil(Math.log10(hi));
    Y = logs(lo, hi, box.b, box.t);
    axes(g, X, Y, { xt, xf, yt: logTicks(lo, hi), yl: "Afρ / Afρ(ref)" });
    clip(g, X, Y);
    SCEN.forEach((sc, k) => path(g, s.t, series[k], X, Y, { c: alpha(css(ICE_VAR[sc.x]), 0.55), w: 1.2, dash: sc.c < 1 ? [5, 4] : [] }));
    path(g, s.t, cur, X, Y, { c: iceC, w: 2.6 });
    const xr = X(ref); g.strokeStyle = mute; g.setLineDash([2, 3]); g.beginPath(); g.moveTo(xr, box.t); g.lineTo(xr, box.b); g.stroke(); g.setLineDash([]);
    g.fillStyle = mute; g.font = FONTS; g.textAlign = "left"; g.fillText("ref", xr + 4, box.b - 6);
    g.restore();
    $("key1").innerHTML = `<span class="key"><span class="sw" style="border-top-color:${iceC}"></span>${two ? "chosen ices" : "chosen ice"}</span> <span class="key"><span class="sw" style="border-top-color:${ghost}"></span>other ices, subsolar</span> <span class="key"><span class="sw d"></span>isothermal</span>`;
  }
  // perihelion, today, date line
  const vline = (ms, col, dash, label, top) => {
    if (ms < s.t0 || ms > s.t1) return; const x = X(ms);
    g.save(); g.strokeStyle = col; g.lineWidth = 1; g.setLineDash(dash); g.beginPath(); g.moveTo(x, box.t); g.lineTo(x, box.b); g.stroke(); g.restore();
    if (label) { g.save(); g.fillStyle = col; g.font = FONTS; const tw = g.measureText(label).width; g.textAlign = x + tw + 8 > box.r ? "right" : "left"; g.textBaseline = "top"; g.fillText(label, g.textAlign === "left" ? x + 4 : x - 4, box.t + 3 + (top || 0)); g.restore(); }
  };
  vline(s.el.Tms, mute, [2, 3], "perihelion", 0);
  vline(Date.now(), ghost, [], "today", 14);
  if (S.view === "mag" && Number.isFinite(cutMs())) vline(cutMs(), css("--obs"), [4, 3], "fit | forecast", box.b - box.t - 20);
  const c = Math.min(Math.max(S.cursor, s.t0), s.t1); S.cursor = c;
  vline(c, acc, [], "", 0);
  g.save(); g.fillStyle = acc; g.beginPath(); const xc = X(c); g.moveTo(xc - 6, box.b); g.lineTo(xc + 6, box.b); g.lineTo(xc, box.b - 8); g.fill(); g.restore();
}
function msg(g, w, h, text) {
  g.save(); g.fillStyle = css("--mute"); g.font = '15px "Source Serif 4", Georgia, serif'; g.textAlign = "center"; g.textBaseline = "middle";
  const words = text.split(" "), lines = [""]; for (const wd of words) { const t = lines[lines.length - 1] + (lines[lines.length - 1] ? " " : "") + wd; if (g.measureText(t).width > w - 60) lines.push(wd); else lines[lines.length - 1] = t; }
  lines.forEach((ln, k) => g.fillText(ln, w / 2, h / 2 + (k - (lines.length - 1) / 2) * 20)); g.restore();
}
function refMs(s) { const v = parseDate($("ref").value); return Number.isFinite(v) ? v : s.t0; }

// ---------- eta and n panels
function curveR(r0, r1, n) { const out = []; for (let k = 0; k <= n; k++) out.push(r0 * (r1 / r0) ** (k / n)); return out; }
function drawEta() {
  const cv = $("c2"), { g, w, h } = frame(cv, 0.68, 230, 360), s = S.series, P = s.P;
  const box = { l: 46, r: w - 12, t: 12, b: h - 38 };
  const X = logs(0.3, 50, box.l, box.r), Y = lin(0, 100, box.b, box.t), rs = curveR(0.3, 50, 260);
  axes(g, X, Y, { xt: [0.3, 1, 3, 10, 30, 50].filter((v) => v <= 50), yt: [0, 25, 50, 75, 100], yf: (v) => v + "%", xl: "distance from the Sun r (au)" });
  clip(g, X, Y);
  const rr = s.r.filter(Number.isFinite);
  if (rr.length) { const a = Math.max(0.3, Math.min(...rr)), b = Math.min(50, Math.max(...rr)); g.fillStyle = alpha(css("--ghost"), 0.16); g.fillRect(X(a), box.t, Math.max(1, X(b) - X(a)), box.b - box.t); }
  const sel = css(ICE_VAR[P.ice]);
  const lo = rs.map((r) => 100 * eta(r, P.ice, 0.25, P.h)), hi = rs.map((r) => 100 * eta(r, P.ice, 1, P.h));
  g.fillStyle = alpha(sel, 0.13); g.beginPath(); rs.forEach((r, k) => k ? g.lineTo(X(r), Y(hi[k])) : g.moveTo(X(r), Y(hi[k]))); for (let k = rs.length - 1; k >= 0; k--) g.lineTo(X(rs[k]), Y(lo[k])); g.fill();
  for (const x of ["CO", "CO2", "H2O"]) {
    const ys = rs.map((r) => 100 * eta(r, x, P.cz, P.h)), col = css(ICE_VAR[x]);
    path(g, rs, ys, X, Y, { c: x === P.ice ? col : alpha(col, 0.7), w: x === P.ice ? 2.6 : 1.4 });
    const k = ys.findIndex((v) => v < 50);
    g.fillStyle = col; g.font = FONTL; g.textBaseline = "middle";
    if (k > 0) { g.textAlign = "left"; g.fillText(ICE_NAME[x], X(rs[k]) + 6, Y(50) - 10); }
    else { g.textAlign = "right"; g.fillText(ICE_NAME[x], box.r - 6, Y(ys[ys.length - 1]) + 12); }
  }
  const mx = pieceAt(s.el, P, S.cursor).mix, two = mx.length > 1;
  if (two) { path(g, rs, rs.map((r) => 100 * etaMix(r, mx, P.cz, P.h)), X, Y, { c: css("--ink"), w: 2.4 }); g.fillStyle = css("--ink"); g.font = FONTL; g.textAlign = "left"; g.textBaseline = "middle"; g.fillText("Σ kη", box.l + 6, Y(100 * etaMix(0.3, mx, P.cz, P.h)) + 12); }
  const rc = ephem(s.el, S.cursor).r;
  if (rc >= 0.3 && rc <= 50) { const x = X(rc); g.strokeStyle = css("--accent"); g.lineWidth = 1; g.beginPath(); g.moveTo(x, box.t); g.lineTo(x, box.b); g.stroke(); g.fillStyle = two ? css("--ink") : sel; g.beginPath(); g.arc(x, Y(100 * etaMix(rc, mx, P.cz, P.h)), 4, 0, 2 * Math.PI); g.fill(); }
  g.restore();
}
function drawN() {
  const cv = $("c3"), { g, w, h } = frame(cv, 0.68, 230, 360), s = S.series, P = s.P;
  const box = { l: 40, r: w - 12, t: 12, b: h - 38 };
  const X = logs(1, 50, box.l, box.r), Y = lin(0, 20, box.b, box.t), rs = curveR(1, 50, 220);
  axes(g, X, Y, { xt: [1, 2, 5, 10, 20, 50], yt: [0, 5, 10, 15, 20], xl: "distance from the Sun r (au)" });
  clip(g, X, Y);
  const mute = css("--mute"), ink = css("--ink");
  g.strokeStyle = mute; g.setLineDash([2, 3]); g.beginPath(); g.moveTo(box.l, Y(3.5)); g.lineTo(box.r, Y(3.5)); g.stroke(); g.setLineDash([]);
  g.fillStyle = mute; g.font = FONTS; g.textAlign = "left"; g.textBaseline = "top"; g.fillText("fully active 7/2", box.l + 4, Y(3.5) + 3);
  if (S.nbox) { const b = S.nbox, oc = css("--obs"); g.fillStyle = alpha(oc, 0.22); g.strokeStyle = oc; g.fillRect(X(b.r0), Y(b.n + b.s), X(b.r1) - X(b.r0), Y(b.n - b.s) - Y(b.n + b.s)); g.strokeRect(X(b.r0), Y(b.n + b.s), X(b.r1) - X(b.r0), Y(b.n - b.s) - Y(b.n + b.s)); g.fillStyle = oc; g.textAlign = "left"; g.textBaseline = "middle"; g.fillText("UN271", X(b.r1) + 5, Y(b.n)); }
  const holt = (r) => 0.40 * r * (Math.log(r) + 1) - 1.78;
  path(g, curveR(2.5, 10, 60), curveR(2.5, 10, 60).map(holt), X, Y, { c: ink, w: 1.4, dash: [6, 4] });
  path(g, curveR(10, 15, 20), curveR(10, 15, 20).map(holt), X, Y, { c: alpha(css("--ghost"), 1), w: 1.2, dash: [2, 3] });
  for (const x of ["CO", "CO2", "H2O"]) {
    const col = css(ICE_VAR[x]), ys = rs.map((r) => 2 + kSlope(r, x, P.cz, P.h));
    path(g, rs, ys, X, Y, { c: x === P.ice ? col : alpha(col, 0.7), w: x === P.ice ? 2.6 : 1.4 });
    let k = ys.findIndex((v) => v > 17); if (k < 0) k = ys.length - 1;
    g.fillStyle = col; g.font = FONTL; g.textAlign = "right"; g.textBaseline = "middle";
    const yy = Number.isFinite(ys[k]) ? Y(Math.min(ys[k], 17)) : Y(17); g.fillText(ICE_NAME[x], X(rs[k]) - 6, x === "CO" ? Y(3.5) - 12 : yy);
  }
  const mx = pieceAt(s.el, P, S.cursor).mix, two = mx.length > 1;
  if (two) path(g, rs, rs.map((r) => 2 + kSlopeMix(r, mx, P.cz, P.h)), X, Y, { c: ink, w: 2.4 });
  const rc = ephem(s.el, S.cursor).r;
  if (rc >= 1 && rc <= 50) { const x = X(rc), nv = 2 + kSlopeMix(rc, mx, P.cz, P.h); g.strokeStyle = css("--accent"); g.lineWidth = 1; g.beginPath(); g.moveTo(x, box.t); g.lineTo(x, box.b); g.stroke(); if (nv <= 20) { g.fillStyle = two ? ink : css(ICE_VAR[P.ice]); g.beginPath(); g.arc(x, Y(nv), 4, 0, 2 * Math.PI); g.fill(); } }
  g.restore();
}

// ---------- readouts and table
function readouts() {
  const s = S.series, P = s.P, el = s.el, cell = (b, v) => `<div><b>${b}</b><span>${v}</span></div>`;
  if (!el.ok) { $("read").innerHTML = cell("Date", "—"); return; }
  const g = ephem(el, S.cursor), pc = pieceAt(el, P, S.cursor), two = pc.mix.length > 1, bal = balance(g.r, P.ice, P.cz, P.h);
  const et = etaMix(g.r, pc.mix, P.cz, P.h), k = kSlopeMix(g.r, pc.mix, P.cz, P.h);
  let mtxt = "—", share = "—";
  if (el.full) { const mg = magnitude(g.r, g.D, g.alpha, et, P, pc.mu); mtxt = mg.m.toFixed(2); share = (100 * mg.nucShare).toFixed(mg.nucShare < 0.01 ? 2 : 1) + " %"; }
  const T2 = two ? balance(g.r, P.ice2, P.cz, P.h).T : NaN;
  $("read").innerHTML = cell("Date", fmtD(S.cursor)) + cell('<span class="m">r</span> (au)', g.r.toFixed(2)) + cell('<span class="m">Δ</span> (au)', el.full ? g.D.toFixed(2) : "—") +
    cell('Phase <span class="m">α</span>', el.full ? g.alpha.toFixed(1) + "°" : "—") + cell('Magnitude <span class="m">m</span>', mi(mtxt)) + cell("Nucleus share", share) +
    (two ? cell(`T ${ICE_NAME[P.ice]} / ${ICE_NAME[P.ice2]}`, `${bal.T.toFixed(0)} / ${T2.toFixed(0)} K`) : cell(Number.isFinite(P.h) ? "Surface / ice T" : "Surface T", Number.isFinite(P.h) ? `${bal.T.toFixed(0)} / ${bal.Ti.toFixed(0)} K` : bal.T.toFixed(0) + " K")) +
    cell(two ? '<span class="m">Σ kη</span>' : '<span class="m">η<sub>' + ICE_NAME[P.ice] + "</sub></span>", sci(et, 3)) +
    cell('<span class="m">n<sub>Q</sub></span> (<span class="m">n<sub>C</sub></span>)', Number.isFinite(k) ? `${(k + 0.5).toFixed(2)} (${(k + 2).toFixed(1)})` : "—") + cell('<span class="m">Afρ</span> (cm)', sci(afrho(pc.K, g.r, et), 3));
}
function table() {
  const s = S.series, P = s.P, el = s.el, tb = $("tab").tBodies[0];
  if (!el.ok) { tb.innerHTML = ""; return; }
  const g = ephem(el, S.cursor), ref = refMs(s), r0 = ephem(el, ref).r, rows = [...SCEN, ...PAIRS];
  let best = null, bs = Infinity; if (S.cmp) for (const sc of rows) { const c = S.cmp[keyOf(sc)]; if (c && (!best || rankOf(c, sc) < bs)) { best = c; bs = rankOf(c, sc); } }
  tb.innerHTML = rows.map((sc) => {
    const fitc = S.cmp && S.cmp[keyOf(sc)], kk = sc.y ? (fitc ? fitc.k : P.ice === sc.x && P.ice2 === sc.y ? P.k : 0.5) : 1;
    const mx = mixFrom(sc.x, sc.y, kk);
    const k = kSlopeMix(g.r, mx, sc.c, P.h), gr = etaMix(g.r, mx, sc.c, P.h) / etaMix(r0, mx, sc.c, P.h) * (g.r / r0) ** -1.5;
    let m = "—";
    if (el.full) { const PP = { ...P, HN: fitc ? fitc.HN : P.HN }, K = fitc ? 10 ** fitc.lk : s.K; m = magnitude(g.r, g.D, g.alpha, etaMix(g.r, mx, sc.c, P.h), PP, muFrom(K, PP)).m.toFixed(2); }
    const on = sc.x === P.ice && (sc.y || "none") === (P.ice2 === P.ice ? "none" : P.ice2) && Math.abs(sc.c - P.cz) < 1e-6, isBest = best && fitc === best;
    const chi = fitc ? (fitc.nu > 0 ? fitc.chi2 / fitc.nu : fitc.chi2).toFixed(fitc.chi2 < 100 ? 2 : 0) + (fitc.nu > 0 ? "" : " (χ²)") : "—";
    const chip = (x) => `<span class="chip" style="background:var(${ICE_VAR[x]})"></span>`;
    const lab = sc.y ? `${chip(sc.x)}${chip(sc.y)}${sc.lab}${fitc ? `, k = ${fitc.k.toFixed(2)}` : ""}` : `${chip(sc.x)}${sc.lab}`;
    return `<tr data-x="${sc.x}" data-y="${sc.y || "none"}" data-c="${sc.c}" data-key="${keyOf(sc)}" class="${on ? "on" : ""} ${isBest ? "best" : ""}" tabindex="0"><td>${lab}</td><td>${Number.isFinite(k) ? (k + 0.5).toFixed(2) : "—"}</td><td>${Number.isFinite(k) ? (k + 2).toFixed(1) : "—"}</td><td>× ${sci(gr, 3)}</td><td>${mi(m)}</td><td>${mi(chi)}</td></tr>`;
  }).join("");
}

// ---------- fitting
function fitOne(obs, el, P, x, c, fitH, y = "none", Hfix = null) {
  const two = y && y !== "none" && y !== x;
  const pts = obs.map((o) => { const g = ephem(el, o.ms); return { ...o, r: g.r, D: g.D, al: g.alpha, eA: eta(g.r, x, c, P.h), eB: two ? eta(g.r, y, c, P.h) : 0 }; });
  const chi = (lk, HN, k) => { const PP = { ...P, HN }, mu = muFrom(10 ** lk, PP); let s2 = 0; for (const p of pts) { const et = two ? k * p.eA + (1 - k) * p.eB : p.eA, d = (p.m - magnitude(p.r, p.D, p.al, et, PP, mu).m) / p.s; s2 += d * d; } return s2; };
  const rmsOf = (lk, HN, k) => { const PP = { ...P, HN }, mu = muFrom(10 ** lk, PP); let s2 = 0; for (const p of pts) { const et = two ? k * p.eA + (1 - k) * p.eB : p.eA, d = p.m - magnitude(p.r, p.D, p.al, et, PP, mu).m; s2 += d * d; } return Math.sqrt(s2 / pts.length); };
  const gr = (Math.sqrt(5) - 1) / 2;
  const golden = (f, a, b, n) => { for (let q = 0; q < n; q++) { const c1 = b - gr * (b - a), c2 = a + gr * (b - a); if (f(c1) < f(c2)) b = c2; else a = c1; } return (a + b) / 2; };
  const best1 = (HN, k) => {
    let bl = 0, bv = Infinity; for (let lk = -1; lk <= 10; lk += 0.1) { const v = chi(lk, HN, k); if (v < bv) { bv = v; bl = lk; } }
    const lk = golden((l) => chi(l, HN, k), bl - 0.1, bl + 0.1, 40); return { lk, v: chi(lk, HN, k) };
  };
  const bestK = (HN) => {                          // share k of the first ice: grid, then golden section
    if (!two) return { k: 1, ...best1(HN, 1) };
    let bk = 0, bv = Infinity; for (let k = 0; k <= 1.0001; k += 0.02) { const v = best1(HN, k).v; if (v < bv) { bv = v; bk = k; } }
    const k = Math.min(1, Math.max(0, golden((kk) => best1(HN, Math.min(1, Math.max(0, kk))).v, bk - 0.02, bk + 0.02, 20)));
    return { k, ...best1(HN, k) };
  };
  const np = 1 + (two ? 1 : 0);
  if (!fitH || Hfix !== null) { const HN = Hfix !== null ? Hfix : P.HN, r = bestK(HN); return { lk: r.lk, k: r.k, HN, chi2: r.v, rms: rmsOf(r.lk, HN, r.k), nu: pts.length - np, score: r.v / Math.max(1, pts.length - np) }; }
  const H0 = Math.max(2, P.HN - 3), H1 = Math.min(22, P.HN + 3);
  let k = two ? bestK(P.HN).k : 1, HN = P.HN;
  for (let it = 0; it < (two ? 3 : 1); it++) {     // with two ices, alternate between H_N and k
    let bH = HN, bv = Infinity; for (let h = H0; h <= H1 + 1e-9; h += 0.1) { const r = best1(h, k); if (r.v < bv) { bv = r.v; bH = h; } }
    HN = golden((h) => best1(h, k).v, bH - 0.1, bH + 0.1, 30);
    if (two) k = bestK(HN).k;
  }
  const r = best1(HN, k);
  return { lk: r.lk, k, HN, chi2: r.v, rms: rmsOf(r.lk, HN, k), nu: pts.length - np - 1, score: r.v / Math.max(1, pts.length - np - 1), edge: HN - H0 < 0.15 || H1 - HN < 0.15 };
}
function obsForFit() {
  const el = readOrbit(); const { obs, bad } = parseObs();
  if (!el.full) return { err: "Fitting needs the full orbit (i, Ω, ω) to know Δ and the phase angle." };
  if (!obs.length) return { err: "Add photometry first: one line per date with magnitude and error." };
  const cut = cutMs(), use = obs.filter((o) => o.ms < cut);
  if (!use.length) return { err: "No points before the “fit only before” date." };
  return { el, obs: use, bad, cut: Number.isFinite(cut) ? fmtD(cut) : null };
}
function fitText(r, n, two) {
  return `K = ${sci(10 ** r.lk, 3)} cm` + (two ? `, k = ${r.k.toFixed(2)}` : "") + ", " +
    (r.nu > 0 ? `χ²<sub>ν</sub> = ${(r.chi2 / r.nu).toFixed(2)}, rms ${r.rms.toFixed(2)} mag for ${n} points` : `χ² = ${r.chi2.toFixed(2)}: as many parameters as points`);
}
function doFit() {
  const f = obsForFit(); if (f.err) return score(f.err, "mid");
  const P = params(), fitH = $("fith").checked, two = mixOf(P).length > 1;
  if (fitH && f.obs.length < 2) return score("Fitting H_N as well needs at least two points.", "mid");
  const pre = P.split ? f.obs.filter((o) => o.ms < f.el.Tms) : f.obs, post = P.split ? f.obs.filter((o) => o.ms >= f.el.Tms) : [];
  if (!pre.length) return score("No points before perihelion to fit; untick “fit before and after perihelion separately”.", "mid");
  const r = fitOne(pre, f.el, P, P.ice, P.cz, fitH, P.ice2);
  $("lk").value = r.lk.toFixed(4); if (two) $("kmix").value = r.k.toFixed(3); if (fitH) $("hn").value = r.HN.toFixed(2);
  let rp = null; S.post = null;
  if (P.split && post.length) { rp = fitOne(post, f.el, P, P.ice, P.cz, false, P.ice2, r.HN); S.post = { lk: rp.lk, k: rp.k }; }
  const PP = params(), head = `${mixName(P)}, ⟨cos ζ⟩ = ${P.cz.toFixed(2)}${fitH ? `, H<sub>N</sub> = ${r.HN.toFixed(2)}` : ""}`;
  let txt;
  if (P.split) txt = `${head}. Before perihelion: ${fitText(r, pre.length, two)}. ` + (rp ? `After: ${fitText(rp, post.length, two)}; ` + (two ? `K changed by × ${(10 ** (rp.lk - r.lk)).toFixed(2)} and k from ${r.k.toFixed(2)} to ${rp.k.toFixed(2)}.` : `the active area changed by × ${(10 ** (rp.lk - r.lk)).toFixed(2)}.`) : "No points after perihelion: the curve after it uses the same constants.");
  else txt = `${head}: ${fitText(r, f.obs.length, two)}, μ = ${sci(muFrom(10 ** r.lk, PP), 3)}.`;
  txt += (f.cut ? ` Points before ${f.cut} only.` : "") + (r.edge ? ` H<sub>N</sub> ran to the edge of its range (±3 mag): these data do not constrain it, so untick “also H<sub>N</sub>”.` : "") + (f.bad.length ? ` Skipped line ${f.bad.join(", ")}.` : "");
  const bad = r.edge || (r.nu > 0 && r.chi2 / r.nu > 3) || (rp && rp.nu > 0 && rp.chi2 / rp.nu > 3);
  score(txt, bad ? "mid" : "good");
  update();
}
function doCompare() {
  const f = obsForFit(); if (f.err) return score(f.err, "mid");
  const P = params(), fitH = $("fith").checked && f.obs.length >= 2, rows = [...SCEN, ...PAIRS]; S.cmp = {};
  for (const sc of rows) S.cmp[keyOf(sc)] = fitOne(f.obs, f.el, P, sc.x, sc.c, fitH, sc.y || "none");
  let best = rows[0]; for (const sc of rows) if (rankOf(S.cmp[keyOf(sc)], sc) < rankOf(S.cmp[keyOf(best)], best)) best = sc;
  const b = S.cmp[keyOf(best)];
  score(`Best: ${best.lab}${best.y ? ` with k = ${b.k.toFixed(2)}` : ""}, ` + (b.nu > 0 ? `χ²<sub>ν</sub> = ${(b.chi2 / b.nu).toFixed(2)}` : `χ² = ${b.chi2.toFixed(2)}`) +
    `. A second ice adds one constant, k: it should lower χ²<sub>ν</sub> clearly to count.` + (P.split ? " (Compare fits all points together, not before and after perihelion.)" : "") + ` Click a row to use it.`, "good");
  update();
}
function score(html, cls) { const e = $("score"); e.innerHTML = html; e.className = "score " + (cls || ""); }

// ---------- power budget and SWAN backtest
const BPRE = { neo: { q: "5.27e30", r: "0.337", rn: "3.5" }, "3i": { q: "4.5e28", r: "1.49", rn: "1.3" } };
function budget() {
  const Q = num($("bq").value), r = num($("br").value), RN = num($("brn").value);
  if (!(Q > 0 && r > 0 && RN > 0)) { $("bout").innerHTML = "Enter Q, r and R<sub>N</sub>."; return; }
  const qm = qMax(RN, r, "H2O", 1), f = Q / qm, rmin = RN * Math.sqrt(f);
  $("bout").innerHTML = `Q<sub>max</sub> = ${sci(qm, 3)} s⁻¹ &nbsp;·&nbsp; f = ${f.toPrecision(3)} &nbsp;<span class="pill ${f <= 1 ? "ok" : "no"}">${f <= 1 ? "power suffices" : "needs another source"}</span><br>smallest nucleus for f ≤ 1: R<sub>N</sub> = ${rmin.toPrecision(3)} km`;
}
function drawSwan() {
  const cv = $("c4"), { g, w, h } = frame(cv, 0.5, 260, 420);
  const box = { l: 44, r: w - 12, t: 12, b: h - 38 };
  const X = logs(0.005, 2, box.l, box.r), Y = lin(-10, 4, box.b, box.t);
  axes(g, X, Y, { xt: [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2], xf: (v) => String(v), yt: [-10, -8, -6, -4, -2, 0, 2, 4], yf: (v) => mi(String(v)), xl: "perihelion distance q (au)" });
  clip(g, X, Y);
  const qs = curveR(0.005, 2, 70), lo = qs.map((q) => envelope(q, 0.25)), hi = qs.map((q) => envelope(q, 1)), mt = qs.map((q) => envelopeMantle(q)[1]);
  const band = (top, bot, col) => { g.fillStyle = col; g.beginPath(); qs.forEach((q, k) => k ? g.lineTo(X(q), Y(top[k])) : g.moveTo(X(q), Y(top[k]))); for (let k = qs.length - 1; k >= 0; k--) g.lineTo(X(qs[k]), Y(bot[k])); g.fill(); };
  band(mt, hi, alpha(css("--h2o"), 0.09)); band(hi, lo, alpha(css("--h2o"), 0.24));
  path(g, qs, hi, X, Y, { c: css("--h2o"), w: 1.2 }); path(g, qs, lo, X, Y, { c: css("--h2o"), w: 1.2, dash: [5, 4] }); path(g, qs, mt, X, Y, { c: css("--h2o"), w: 1.2, dash: [1.5, 3] });
  const good = css("--good"), mid = css("--warn"), bad = css("--obs"); let n = 0, ok = 0, okm = 0; S.swanPts = [];
  if (!S.swanEnv) { S.swanEnv = {}; for (const r of SWAN) S.swanEnv[r[1]] = [envelope(r[1], 0.25), envelope(r[1], 1), envelopeMantle(r[1])[1]]; }
  for (const [name, q, pp, ep, pq, eq] of SWAN) {
    const [L, H, HM] = S.swanEnv[q];
    [[pp, ep, 0], [pq, eq, 1]].forEach(([p, e, post]) => {
      if (p === null || e === null) return;
      const pass = p + 2 * e >= L && p - 2 * e <= H, passm = p + 2 * e >= L && p - 2 * e <= HM; n++; ok += pass; okm += passm;
      const col = pass ? good : passm ? mid : bad, x = X(q) + (post ? 3 : -3), y = Y(Math.max(-9.8, p));
      g.strokeStyle = col; g.fillStyle = col; g.lineWidth = 1;
      if (p >= -10) { g.beginPath(); g.moveTo(x, Y(Math.max(-10, p - e))); g.lineTo(x, Y(Math.min(4, p + e))); g.stroke(); }
      g.beginPath();
      if (p < -10) { g.moveTo(x, box.b - 2); g.lineTo(x - 4, box.b - 9); g.lineTo(x + 4, box.b - 9); g.closePath(); g.fill(); }
      else if (post) { g.rect(x - 3, y - 3, 6, 6); g.fill(); } else { g.arc(x, y, 3.3, 0, 2 * Math.PI); g.fill(); }
      S.swanPts.push({ x, y: p < -10 ? box.b - 6 : y, name, q, p, e, post, pass, passm, L, H, HM });
    });
  }
  // medians in bins of q with 68% bootstrap intervals (as in the backtest of the SWAN exponents)
  g.strokeStyle = css("--ink"); g.fillStyle = css("--panel"); g.lineWidth = 1.8;
  for (const [qm, m, p16, p84] of SWAN_MED) {
    const x = X(qm), y = Y(m);
    g.beginPath(); g.moveTo(x, Y(p16)); g.lineTo(x, Y(p84)); g.stroke();
    g.beginPath(); g.rect(x - 5, y - 5, 10, 10); g.fill(); g.stroke();
  }
  g.restore();
  $("swansum").innerHTML = `<span class="key"><span class="dot" style="background:${good}"></span>${ok} of ${n} pass</span> <span class="key"><span class="dot" style="background:${mid}"></span>${okm - ok} more with a mantle</span>`;
}
function swanHover(ev) {
  const cv = $("c4"), rc = cv.getBoundingClientRect(), x = ev.clientX - rc.left, y = ev.clientY - rc.top, tip = $("tip");
  let best = null, bd = 144; for (const p of S.swanPts || []) { const d = (p.x - x) ** 2 + (p.y - y) ** 2; if (d < bd) { bd = d; best = p; } }
  if (!best) { tip.hidden = true; return; }
  tip.hidden = false;
  tip.innerHTML = `${best.name}<br>q = ${best.q.toFixed(3)} au, ${best.post ? "after" : "before"} perihelion<br>p = ${mi(best.p.toFixed(1))} ± ${best.e} · band ${mi(best.L.toFixed(2))} … ${mi(best.H.toFixed(2))}<br>${best.pass ? "passes" : best.passm ? `passes with a dust mantle (top ${mi(best.HM.toFixed(2))})` : best.p > best.HM ? "flatter even with a mantle" : "steeper: the active area must grow"}`;
  const fr = cv.parentNode.getBoundingClientRect(), tx = ev.clientX - fr.left + 14, ty = ev.clientY - fr.top + 14;
  tip.style.left = Math.min(tx, fr.width - tip.offsetWidth - 8) + "px"; tip.style.top = Math.min(ty, fr.height - tip.offsetHeight - 8) + "px";
}

// ---------- wiring
function labels() {
  const P = params();
  $("czo").textContent = P.cz.toFixed(2) + (P.cz === 1 ? " subsolar" : P.cz === 0.25 ? " isothermal" : "");
  const D = 1329 / Math.sqrt(0.04) * 10 ** (-P.HN / 5);
  $("hno").textContent = `${P.HN.toFixed(2)} · D ≈ ${D < 10 ? D.toFixed(1) : D.toFixed(0)} km`;
  $("betao").textContent = P.beta.toFixed(3) + " mag/°";
  $("hmo").textContent = Number.isFinite(P.h) ? `h = ${P.h < 1 ? P.h.toPrecision(2) : P.h.toPrecision(3).replace(/\.0+$/, "")} W m⁻² K⁻¹` : "none (bare ice)";
  $("lko").innerHTML = sci(10 ** P.lk, 3) + " cm";
  const two = mixOf(P).length > 1; $("kbox").hidden = !two;
  if (two) $("kmixo").textContent = `${P.k.toFixed(2)} / ${(1 - P.k).toFixed(2)}`;
  const mu = muFrom(10 ** P.lk, P);
  $("munote").innerHTML = `μ = ${sci(mu, 3)} for this aperture, so m<sub>1</sub> = H<sub>N</sub> − 2.5 log μ = ${(P.HN - 2.5 * Math.log10(mu)).toFixed(2).replace("-", "−")} (Eq. 14). D for p = 0.04.`;
  const { obs, bad } = parseObs(); $("nobs").textContent = `${obs.length} point${obs.length === 1 ? "" : "s"}${bad.length ? `, ${bad.length} unreadable` : ""}`;
  document.querySelectorAll(".seg label").forEach(() => {});
  $("vmag").setAttribute("aria-pressed", S.view === "mag"); $("vafr").setAttribute("aria-pressed", S.view === "afr");
  $("vmag").disabled = !readOrbit().full;
}
let raf = 0;
function update() { if (raf) return; raf = requestAnimationFrame(() => { raf = 0; labels(); compute(); drawLC(); drawEta(); drawN(); readouts(); table(); }); }
function applyPreset(key) {
  const p = PRESETS[key];
  $("orbitbox").open = key === "custom" || key === "y1";
  if (!p) { $("cometnote").innerHTML = "Enter the orbit from the JPL Small-Body Database or the MPC, set H<sub>N</sub>, the band and the aperture of your photometry, paste it below and fit."; S.nbox = null; S.cmp = null; S.post = null; update(); return; }
  const o = p.orbit; $("oq").value = o.q; $("oe").value = o.e; $("oi").value = o.i; $("oO").value = o.Om; $("ow").value = o.w; $("oT").value = o.T;
  $("ice-" + p.ice).checked = true; $("ice2-" + (p.ice2 || "none")).checked = true; $("kmix").value = p.k ?? 0.5;
  $("split").checked = !!p.split; S.post = p.post || null; $("cz").value = p.cz; $("hm").value = 3; $("hn").value = p.HN; $("beta").value = p.beta; $("lk").value = p.lk.toFixed(4);
  $("j").value = String(p.j); $("ap").value = p.ap; $("msun").value = p.msun; $("y0").value = p.y0; $("y1").value = p.y1; $("ref").value = p.ref;
  $("obs").value = p.obs; $("cut").value = p.cut || ""; $("cometnote").innerHTML = p.note; S.view = p.view; S.nbox = p.nbox || null; S.cmp = null;
  const now = Date.now(), [t0, t1] = windowMs(); S.cursor = p.cursor ? parseDate(p.cursor) : p.view === "mag" && now >= t0 && now <= t1 ? now : parseDate(o.T);
  score(""); update();
}
function cursorFrom(ev) { const b = S.box1; if (!b) return; const rc = $("c1").getBoundingClientRect(); const x = Math.min(Math.max(ev.clientX - rc.left, b.box.l), b.box.r); S.cursor = b.X.inv(x); update(); }
function init() {
  $("comet").addEventListener("change", (e) => applyPreset(e.target.value));
  ["oq", "oe", "oi", "oO", "ow", "oT", "ap", "msun", "y0", "y1", "ref", "j"].forEach((id) => $(id).addEventListener("input", () => { S.cmp = null; S.post = null; update(); }));
  ["cz", "hn", "beta", "lk", "hm"].forEach((id) => $(id).addEventListener("input", update));
  document.querySelectorAll('input[name="ice"], input[name="ice2"]').forEach((e) => e.addEventListener("change", () => { S.post = null; update(); }));
  $("kmix").addEventListener("input", update);
  $("split").addEventListener("change", () => { if (!$("split").checked) score(""); update(); });
  $("obs").addEventListener("input", () => { S.cmp = null; S.post = null; update(); });
  $("cut").addEventListener("input", () => { S.cmp = null; S.post = null; update(); }); $("emp").addEventListener("change", update);
  $("fit").addEventListener("click", doFit); $("cmp").addEventListener("click", doCompare);
  $("vmag").addEventListener("click", () => { S.view = "mag"; update(); }); $("vafr").addEventListener("click", () => { S.view = "afr"; update(); });
  const c1 = $("c1"); let drag = false;
  c1.addEventListener("pointerdown", (e) => { drag = true; c1.setPointerCapture(e.pointerId); cursorFrom(e); });
  c1.addEventListener("pointermove", (e) => { if (drag) cursorFrom(e); });
  c1.addEventListener("pointerup", () => { drag = false; }); c1.addEventListener("pointercancel", () => { drag = false; });
  c1.addEventListener("keydown", (e) => { if (e.key !== "ArrowLeft" && e.key !== "ArrowRight") return; e.preventDefault(); const [t0, t1] = windowMs(); S.cursor += (e.key === "ArrowRight" ? 1 : -1) * (t1 - t0) / (e.shiftKey ? 20 : 200); update(); });
  $("tab").tBodies[0].addEventListener("click", (e) => pick(e.target.closest("tr")));
  $("tab").tBodies[0].addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); pick(e.target.closest("tr")); } });
  $("bpre").addEventListener("change", (e) => { const p = BPRE[e.target.value]; if (p) { $("bq").value = p.q; $("br").value = p.r; $("brn").value = p.rn; } budget(); });
  ["bq", "br", "brn"].forEach((id) => $(id).addEventListener("input", () => { $("bpre").value = "own"; budget(); }));
  $("c4").addEventListener("pointermove", swanHover); $("c4").addEventListener("pointerleave", () => { $("tip").hidden = true; });
  const redrawAll = () => { update(); requestAnimationFrame(drawSwan); };
  if (window.ResizeObserver) { let wd = 0; new ResizeObserver((en) => { const nw = en[0].contentRect.width; if (Math.abs(nw - wd) > 1) { wd = nw; redrawAll(); } }).observe(document.querySelector(".wrap")); } else window.addEventListener("resize", redrawAll);
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redrawAll);
  new MutationObserver(redrawAll).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(redrawAll);
  $("bq").value = BPRE.neo.q; $("br").value = BPRE.neo.r; $("brn").value = BPRE.neo.rn; budget();
  applyPreset("un271"); requestAnimationFrame(drawSwan);
}
function pick(tr) {
  if (!tr) return; const x = tr.dataset.x, y = tr.dataset.y || "none", c = +tr.dataset.c;
  $("ice-" + x).checked = true; $("ice2-" + y).checked = true; $("cz").value = c; S.post = null;
  const f = S.cmp && S.cmp[tr.dataset.key]; if (f) { $("lk").value = f.lk.toFixed(4); $("hn").value = f.HN.toFixed(2); if (y !== "none") $("kmix").value = f.k.toFixed(3); }
  update();
}
init();
