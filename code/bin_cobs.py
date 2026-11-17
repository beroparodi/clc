"""Bin raw COBS records into light curves with one fixed rule (note on the forecast test).

Input: cobs/raw/<key>.txt, one record per line: 'YYYY-MM-DDTHH:MM magnitude method' (COBS Observation
List API; our records were retrieved 2026 November 17; CC BY-NC-SA 4.0).  Rule, fixed before any forecast was made:
  * methods S, M, B, I (visual) and Z, V, C (CCD), magnitudes as reported (no aperture correction);
  * observations within the apparition window (WINDOWS in forecast_test.py: all records for these comets);
  * only records between the 1st and 99th percentile of their dates (this drops stray records far from
    the apparition, which otherwise stretch the bins: NEOWISE has records up to 2023);
  * bin width w = max(3 d, span/50), span = first to last kept record, bins counted from the first;
  * per bin: median magnitude, dated at the median time, at least 2 observations,
    error sqrt(0.15^2 + (1.4826 MAD)^2 / N); bins with an error above 0.6 mag are dropped.
Usage: python bin_cobs.py --fetch     (downloads the records of all comets but C/2014 UN271; needs internet)
       python bin_cobs.py [key ...]   (prints the light curves in the format of cobs_lightcurves.txt)
COBS keeps receiving and revising records, so a later download may differ slightly from ours."""
import os, sys
import numpy as np
from datetime import datetime, timedelta

HERE = os.path.dirname(os.path.abspath(__file__))
METHODS = set("SMBIZVC")

def read_raw(key, start=None, end=None, path=None):
    out = []
    for line in open(path or os.path.join(HERE, "cobs", "raw", key + ".txt")):
        p = line.split()
        if len(p) < 3 or p[2] not in METHODS: continue
        t = datetime.fromisoformat(p[0])
        if (start and t < start) or (end and t > end): continue
        try: out.append((t, float(p[1])))
        except ValueError: continue
    return sorted(out)

def fetch(key, des):
    """Download all COBS records of one comet (Observation List API) to cobs/raw/<key>.txt."""
    import json
    from urllib.parse import quote
    from urllib.request import urlopen
    lines, page = [], 1
    while True:
        url = f"https://cobs.si/api/obs_list.api?format=json&des={quote(des)}&page={page}"
        j = json.load(urlopen(url, timeout=120))
        for o in j.get("objects", []):
            method = (o.get("obs_method") or {}).get("key", "-")
            lines.append(f"{o['obs_date'][:16].replace(' ', 'T')} {o['magnitude']} {method}")
        if not j.get("info") or page >= j["info"]["pages"]: break
        page += 1
    os.makedirs(os.path.join(HERE, "cobs", "raw"), exist_ok=True)
    open(os.path.join(HERE, "cobs", "raw", key + ".txt"), "w").write("\n".join(lines) + "\n")
    return len(lines)

def bin_records(rec):
    N = len(rec)
    lo, hi = rec[int(np.floor(0.01 * (N - 1) + 0.5))][0], rec[int(np.floor(0.99 * (N - 1) + 0.5))][0]
    rec = [x for x in rec if lo <= x[0] <= hi]
    t0, t1 = rec[0][0], rec[-1][0]
    w = max(3, int(np.floor((t1 - t0).total_seconds() / 86400 / 50 + 0.5)))
    bins = {}
    for t, m in rec:
        bins.setdefault(int((t - t0).total_seconds() // (86400 * w)), []).append((t, m))
    rows = []
    for k in sorted(bins):
        b = bins[k]
        if len(b) < 2: continue
        ms = np.array([m for t, m in b]); ts = sorted(t for t, m in b)
        med = float(np.median(ms)); mad = float(np.median(np.abs(ms - med)))
        err = float(np.sqrt(0.15**2 + (1.4826 * mad) ** 2 / len(ms)))
        if err > 0.6: continue
        tm = ts[len(ts) // 2] if len(ts) % 2 else ts[len(ts) // 2 - 1] + (ts[len(ts) // 2] - ts[len(ts) // 2 - 1]) / 2
        rows.append((tm, med, err))
    return w, rows

if __name__ == "__main__":
    from forecast_test import COMETS, WINDOWS, COBS_DES
    if sys.argv[1:2] == ["--fetch"]:
        for key, des in COBS_DES.items(): print(key, des, fetch(key, des), "records")
        sys.exit()
    keys = sys.argv[1:] or [k for k in WINDOWS]
    for key in keys:
        start, end = WINDOWS[key]
        if not os.path.exists(os.path.join(HERE, "cobs", "raw", key + ".txt")):
            sys.exit(f"no records for {key}: run 'python bin_cobs.py --fetch' first (internet needed)")
        rec = read_raw(key, start, end)
        w, rows = bin_records(rec)
        print(f"#{key}")
        print(f"% {COMETS[key]['name']}: {len(rec)} records, {rec[0][0].date()} to {rec[-1][0].date()}, "
              f"bins of {w} d, {len(rows)} bins")
        for tm, med, err in rows:
            print(f"{tm:%Y-%m-%d} {med:.2f} {err:.2f}")
