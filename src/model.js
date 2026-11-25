// ---- Model core: energy balance (Eqs. 8–11 of the paper), the physical law (Eq. 15), Afrho, orbits. Pure functions; also runs in Node for tests.
const AU_CM = 1.495978707e13, ARCSEC = Math.PI / 180 / 3600, DEG = Math.PI / 180;
const S0 = 1361.0, SIG = 5.670374e-8, KB = 1.380649e-23, AMU = 1.66054e-27;
const A_B = 0.01, EPS = 0.9, KG = 0.01720209895;
// A_X [Pa], B_X [K], molar mass [g/mol]  (Table I of the paper)
const ICES = { H2O: [3.56e12, 6141.667, 18.015], CO2: [1.23e12, 3167.8, 44.01], CO: [1.26e9, 764.16, 28.01] };

function latent(x) { const [, B, mu] = ICES[x]; return B * KB / (mu * AMU); }
function zflux(T, x) { const [A, B, mu] = ICES[x]; const m = mu * AMU; return A * Math.exp(-B / T) * Math.sqrt(m / (2 * Math.PI * KB * T)); }

// Surface temperature T, ice temperature Ti, mass flux Z [kg m^-2 s^-1] and energy share eta for ice x
// at r [au], <cos zeta> = c, under a dust mantle of conductance h [W m^-2 K^-1] (h = Infinity: bare ice,
// Eq. 8; otherwise surface: fin = eps sigma T^4 + h (T - Ti), ice front: h (T - Ti) = Z(Ti) L)
function balance(r, x, c, h = Infinity) {
  const fin = (1 - A_B) * S0 * c / (r * r), L = latent(x), bare = !Number.isFinite(h);
  let lo = 5, hi = 400;
  for (let k = 0; k < 100; k++) {
    const Ti = 0.5 * (lo + hi), zl = zflux(Ti, x) * L, T = bare ? Ti : Ti + zl / h;
    if (fin - EPS * SIG * T ** 4 - zl > 0) lo = Ti; else hi = Ti;
    if (hi - lo < 1e-12 * Ti) break;
  }
  const Ti = 0.5 * (lo + hi), Z = zflux(Ti, x), zl = Z * L;
  return { T: bare ? Ti : Ti + zl / h, Ti, Z, eta: zl / fin };
}
const _eta = new Map();
function eta(r, x, c, h = Infinity) {
  const key = x + c.toFixed(4) + r.toPrecision(10) + h;
  let v = _eta.get(key);
  if (v === undefined) { v = balance(r, x, c, h).eta; if (_eta.size > 300000) _eta.clear(); _eta.set(key, v); }
  return v;
}
// local slope k = 3/2 - dln(eta)/dln(r); activity index n = 2 + k = 3/2 + n_Q (Eqs. 12, 13)
function kSlope(r, x, c, hm = Infinity) {
  const h = 1e-3, a = eta(r * (1 + h), x, c, hm), b = eta(r * (1 - h), x, c, hm);
  if (!(a > 1e-300 && b > 1e-300)) return NaN;
  return 1.5 - (Math.log(a) - Math.log(b)) / (Math.log(1 + h) - Math.log(1 - h));
}

// Several ices on separate patches (Sec. IV of the paper): mix = [{x, w}] with shares w summing to 1.
// etaMix is the bracket sum_i k_i eta_i of Eq. (15); kSlopeMix its local slope, n = 2 + k.
function etaMix(r, mix, c, h = Infinity) { let v = 0; for (const p of mix) if (p.w > 0) v += p.w * eta(r, p.x, c, h); return v; }
function kSlopeMix(r, mix, c, hm = Infinity) {
  if (mix.length === 1) return kSlope(r, mix[0].x, c, hm);
  const h = 1e-3, a = etaMix(r * (1 + h), mix, c, hm), b = etaMix(r * (1 - h), mix, c, hm);
  if (!(a > 1e-300 && b > 1e-300)) return NaN;
  return 1.5 - (Math.log(a) - Math.log(b)) / (Math.log(1 + h) - Math.log(1 - h));
}

// Eq. (15), written with the nucleus as reference: m1 = HN - 2.5 lg mu, nu = 10^(-0.4 beta alpha) / mu.
// s = {HN, beta, j, ap, msun}; ap in arcsec (j = 1) or km (j = 0); K = Afrho scale [cm]
function muFrom(K, s) {
  const f = Math.pow(10, 0.4 * (s.HN - s.msun));
  return s.j === 1 ? K * s.ap * ARCSEC / (4 * AU_CM) * f : K * s.ap * 1e5 / (4 * AU_CM * AU_CM) * f;
}
function magnitude(r, D, alpha, et, s, mu) {
  const ph = alpha * s.beta, mN = s.HN + 5 * Math.log10(r * D) + ph;
  const x = mu * et * Math.pow(r, -1.5) * (s.j === 1 ? D : 1) * Math.pow(10, 0.4 * ph);
  return { m: mN - 2.5 * Math.log10(1 + x), mN, nucShare: 1 / (1 + x) };
}
function afrho(K, r, et) { return K * et * Math.pow(r, -1.5); }

// ---- Orbits (elements referred to the ecliptic and equinox J2000)
function jdOf(ms) { return ms / 86400000 + 2440587.5; }
function msOf(jd) { return (jd - 2440587.5) * 86400000; }
// heliocentric distance and true anomaly from q, e, Tp (JD) at JD J
function anomaly(el, J) {
  const q = el.q, e = el.e, dt = J - el.Tp;
  if (Math.abs(e - 1) < 1e-9) {                       // parabola: Barker's equation
    const W = 3 * KG * dt / Math.sqrt(2 * q ** 3), Y = Math.cbrt(W / 2 + Math.sqrt(W * W / 4 + 1)), D = Y - 1 / Y;
    return { r: q * (1 + D * D), nu: 2 * Math.atan(D) };
  }
  if (e < 1) {                                        // ellipse: E - e sin E = M
    const a = q / (1 - e), n = KG / a ** 1.5;
    let M = n * dt; M = M - 2 * Math.PI * Math.round(M / (2 * Math.PI));
    let lo = -Math.PI, hi = Math.PI;                  // f(E) increases monotonically
    for (let k = 0; k < 200; k++) { const E = 0.5 * (lo + hi); if (E - e * Math.sin(E) - M > 0) hi = E; else lo = E; if (hi - lo < 1e-15) break; }
    let E = 0.5 * (lo + hi);
    for (let k = 0; k < 3; k++) E -= (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
    const x = a * (Math.cos(E) - e), y = a * Math.sqrt(1 - e * e) * Math.sin(E);
    return { r: Math.hypot(x, y), nu: Math.atan2(y, x) };
  }
  const a = q / (e - 1), n = KG / a ** 1.5, M = n * dt;   // hyperbola: e sinh H - H = M
  let lo = -Math.asinh(Math.abs(M) / (e - 1)) - 1, hi = -lo;
  for (let k = 0; k < 300; k++) { const H = 0.5 * (lo + hi); if (e * Math.sinh(H) - H - M > 0) hi = H; else lo = H; if (hi - lo < 1e-15 * Math.max(1, Math.abs(H))) break; }
  let H = 0.5 * (lo + hi);
  for (let k = 0; k < 3; k++) H -= (e * Math.sinh(H) - H - M) / (e * Math.cosh(H) - 1);
  const x = a * (e - Math.cosh(H)), y = a * Math.sqrt(e * e - 1) * Math.sinh(H);
  return { r: Math.hypot(x, y), nu: Math.atan2(y, x) };
}
function earthHelio(J) {                              // low-precision Sun (Meeus), ~0.01 deg
  const T = (J - 2451545.0) / 36525, L0 = (280.46646 + 36000.76983 * T) % 360, M = (357.52911 + 35999.05029 * T) % 360;
  const Mr = M * DEG, C = (1.914602 - 0.004817 * T) * Math.sin(Mr) + 0.019993 * Math.sin(2 * Mr) + 0.000289 * Math.sin(3 * Mr);
  const lam = (L0 + C) * DEG, ecc = 0.016708634 - 0.000042037 * T;
  const R = 1.000001018 * (1 - ecc * ecc) / (1 + ecc * Math.cos(Mr + C * DEG));
  return [-R * Math.cos(lam), -R * Math.sin(lam), 0];
}
// r, Delta, phase angle [deg], elongation [deg] at time ms; el = {q, e, i, Om, w, Tp, full}
function ephem(el, ms) {
  const J = jdOf(ms), { r, nu } = anomaly(el, J);
  if (!el.full) return { r, D: NaN, alpha: NaN, elong: NaN };
  const i = el.i * DEG, O = el.Om * DEG, u = nu + el.w * DEG;
  const C = [r * (Math.cos(O) * Math.cos(u) - Math.sin(O) * Math.sin(u) * Math.cos(i)),
             r * (Math.sin(O) * Math.cos(u) + Math.cos(O) * Math.sin(u) * Math.cos(i)),
             r * Math.sin(u) * Math.sin(i)];
  const E = earthHelio(J), G = [C[0] - E[0], C[1] - E[1], C[2] - E[2]];
  const D = Math.hypot(G[0], G[1], G[2]), R = Math.hypot(E[0], E[1], E[2]);
  const cl = (v) => Math.min(1, Math.max(-1, v));
  return { r, D, alpha: Math.acos(cl((r * r + D * D - R * R) / (2 * r * D))) / DEG,
           elong: Math.acos(cl((R * R + D * D - r * r) / (2 * R * D))) / DEG };
}

// ---- Energy budget (criterion ii): largest water production of a projected disk of radius RN [km]
function qMax(RN, r, x, c) {
  const F = (1 - A_B) * S0 / (r * r), m = ICES[x][2] * AMU;
  return Math.PI * (RN * 1e3) ** 2 * F * eta(r, x, c) / latent(x) / m;
}
// SWAN envelope: mean exponent of constant-area production between q and max(2, q + 1) au
function envelope(q, c, h = Infinity) {
  const r2 = Math.max(2, q + 1), lq = (r) => Math.log(eta(r, "H2O", c, h)) - 2 * Math.log(r);
  return (lq(r2) - lq(q)) / (Math.log(r2) - Math.log(q));
}
// range of the exponent when the ice may lie under a mantle of any h >= 0.03 W m^-2 K^-1 (note on the mantle)
const H_MANTLE = [Infinity, 100, 30, 10, 3, 1, 0.3, 0.1, 0.03];
function envelopeMantle(q) {
  let lo = Infinity, hi = -Infinity;
  for (const h of H_MANTLE) for (const c of [0.25, 1]) { const p = envelope(q, c, h); lo = Math.min(lo, p); hi = Math.max(hi, p); }
  return [lo, hi];
}

if (typeof module !== "undefined") module.exports = { latent, balance, eta, kSlope, etaMix, kSlopeMix, muFrom, magnitude, afrho, ephem, anomaly, jdOf, msOf, qMax, envelope, envelopeMantle, ICES, ARCSEC, AU_CM };
