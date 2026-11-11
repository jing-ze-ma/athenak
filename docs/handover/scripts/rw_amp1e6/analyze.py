#!/usr/bin/env python3
"""rw_amp1e6-1111: IDORT linear-wave (P, tau) grid at amp 1e-6, PROD vs FINAL; orders raw and E mean-removed.
   python3 analyze.py  -> markdown to stdout.  Reference at amp 1e-5: PROD19 Oct (idort1020/rw/runs)."""
import glob, json, os, sys, math
import numpy as np
sys.path.insert(0, os.environ.get("RW_HARNESS", os.path.dirname(os.path.abspath(__file__))))
import run_radwave as rw  # noqa
B = os.environ.get("RWA_ROOT", "/viper/ptmp2/jinma/rw_amp1e6_1111")
ARMS = {"PROD": B + "/runs_PROD", "FINAL": B + "/runs_FINAL", "P19_1e-5": os.environ.get("RWA_REF", "/raven/ptmp/jinma/idort1020/rw/runs"), "PROD_tt": B + "/runs_ttPROD", "FINAL_tt": B + "/runs_ttFINAL"}
G = {}
for arm, D in ARMS.items():
    for p in glob.glob(D + "/*/results.json"):
        r = json.load(open(p)); c = r["case"]
        if c["dir"] != "x1" or c["nt"] != 2 * c["N"] or c["extra"].get("arm") not in ("PROD19", "TIME"):
            continue
        if arm == "P19_1e-5" and (c["extra"].get("arm") != "PROD19" or "amp" in c["extra"] or len(c["extra"]) != 4):
            continue
        bg, refs = rw.references(c); w, vec = refs[r["ref_key"]]
        x = np.array(r["x1"]); e = np.array(r["e_end"]); rho = np.array(r["rho_end"]); t = r["t_end"]
        A = float(c["extra"].get("amp", 1e-5)); k = bg.k
        eex = bg.e0 * (1 + A * np.real(vec[3] / bg.e0 * np.exp(1j * k * x - 1j * w * t)))
        rex = bg.rho0 * (1 + A * np.real(np.exp(1j * k * x - 1j * w * t)))
        nE = A * abs(vec[3])
        G.setdefault((arm, round(c["prat"], 4), c["tau"]), {})[c["N"]] = dict(
            r=np.mean(np.abs(rho - rex)) / (A * bg.rho0), E=np.mean(np.abs(e - eex)) / nE,
            Enm=np.mean(np.abs((e - e.mean()) - (eex - eex.mean()))) / nE, mo=(e.mean() - bg.e0) / nE)


def p(d, v, a=128, b=256):
    return math.log(d[a][v] / d[b][v]) / math.log(b / a) if a in d and b in d else float("nan")


f = lambda x: "%.2f" % x if np.isfinite(x) else "-"  # noqa
print("| P | tau | arm | p rho | p E raw | p E mean-rm | p E raw 64->128 | L1 rho 256 | L1 E 256 | L1 E-mrm 256 | mean off 256 (A|dE|) | done N |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
cells = sorted({(k[1], k[2]) for k in G})
for P, tau in cells:
    for arm in ("PROD", "FINAL", "PROD_tt", "FINAL_tt", "P19_1e-5"):
        d = G.get((arm, P, tau), {})
        if not d:
            continue
        g = lambda v: "%.3e" % d[256][v] if 256 in d else "-"  # noqa
        print("| %.4g | %g | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            P, tau, arm, f(p(d, "r")), f(p(d, "E")), f(p(d, "Enm")), f(p(d, "E", 64, 128)),
            g("r"), g("E"), g("Enm"), g("mo"), ",".join(map(str, sorted(d)))))
# FINAL / PROD ratio at N256
print("\nFINAL/PROD L1 at N256 (rho, E):")
for P, tau in cells:
    a, b = G.get(("PROD", P, tau), {}), G.get(("FINAL", P, tau), {})
    if 256 in a and 256 in b:
        print("  P %.4g tau %g: %.3f %.3f" % (P, tau, b[256]["r"] / a[256]["r"], b[256]["E"] / a[256]["E"]))
