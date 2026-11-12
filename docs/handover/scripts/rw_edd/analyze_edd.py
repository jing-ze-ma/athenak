#!/usr/bin/env python3
"""rw_edd_1112: Eddington-closure (D = I/3) linear wave vs the analytic lab-frame Eddington root (radwave_disp 'edd').
   python3 analyze.py -> markdown table: L1 rho / E (raw, mean-removed) orders, omega error vs edd (and eddc)."""
import glob, json, math, sys
import numpy as np
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "rw_amp1e6"))
import run_radwave as rw  # noqa
B = os.environ.get("RWE_ROOT", os.getcwd())
G = {}
for arm in ("PROD", "FINAL"):
    for p in glob.glob(B + "/runs_%s/*/results.json" % arm):
        r = json.load(open(p)); c = r["case"]
        bg, refs = rw.references(c); w, vec = refs["edd"]
        x = np.array(r["x1"]); e = np.array(r["e_end"]); rho = np.array(r["rho_end"]); t = r["t_end"]
        A = float(c["extra"]["amp"]); k = bg.k
        eex = bg.e0 * (1 + A * np.real(vec[3] / bg.e0 * np.exp(1j * (k * x - w * t))))
        rex = bg.rho0 * (1 + A * np.real(np.exp(1j * (k * x - w * t))))
        nE = A * abs(vec[3])
        G.setdefault((arm, round(c["prat"], 4), c["tau"]), {})[c["N"]] = dict(
            r=np.mean(np.abs(rho - rex)) / (A * bg.rho0), E=np.mean(np.abs(e - eex)) / nE,
            Enm=np.mean(np.abs((e - e.mean()) - (eex - eex.mean()))) / nE,
            dw=complex(r["w_re"] - w.real, r["w_im"] - w.imag) / abs(w),
            dwc=complex(r["w_re"] - refs["eddc"][0].real, r["w_im"] - refs["eddc"][0].imag) / abs(w))


def p(d, v, a, b):
    return math.log(d[a][v] / d[b][v]) / math.log(b / a) if a in d and b in d else float("nan")


f = lambda x: "%.2f" % x if np.isfinite(x) else "-"  # noqa
print("| P | tau | arm | p rho 16-32-64-128-256 | p E 16..256 | p E mean-rm 128->256 | L1 rho 256 | L1 E 256 | |dw|/|w| 256 vs edd | vs eddc | done N |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for P, tau in sorted({(k[1], k[2]) for k in G}):
    for arm in ("PROD", "FINAL"):
        d = G.get((arm, P, tau), {})
        if not d:
            continue
        Ns = sorted(d); nm = Ns[-1]
        pr = " ".join(f(p(d, "r", a, 2 * a)) for a in (16, 32, 64, 128))
        pe = " ".join(f(p(d, "E", a, 2 * a)) for a in (16, 32, 64, 128))
        print("| %g | %g | %s | %s | %s | %s | %.3e | %.3e | %.2e | %.2e | %s |" % (
            P, tau, arm, pr, pe, f(p(d, "Enm", 128, 256)), d[nm]["r"], d[nm]["E"],
            abs(d[nm]["dw"]), abs(d[nm]["dwc"]), ",".join(map(str, Ns))))
