"""wbsgate-1112: arm vs arm+wbsolve per case.  dev = max|rho_W - rho_B| / (amp rho0) and max|E_W - E_B|/(amp |E1|)
at t = P (rho_end/e_end of results.json), dw = |w_W - w_B|/|w_ref| (1-mode); prediction g^4, g = (4/3) E0/(rho c^2).
usage: RW_C=10 python3 cmp_wbs.py RUNS_BASE LIST_BASE RUNS_W LIST_W"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "rw_amp1e6") if os.path.isdir(os.path.join(HERE, "rw_amp1e6")) else os.path.join(HERE, "..", "rw_amp1e6"))
import run_radwave as rw  # noqa
rb, lb, rwd, lw = sys.argv[1:5]
C = rw.C_LIGHT
rows, worst = [], {}
for b, w in zip(open(lb).read().split("\n"), open(lw).read().split("\n")):
    if not b.strip():
        continue
    cb, cw = rw.parse_case(b), rw.parse_case(w)
    pb = os.path.join(rb, rw.case_tag(cb), "results.json"); pw = os.path.join(rwd, rw.case_tag(cw), "results.json")
    key = (cb["clos"], cb["prat"], cb["tau"], cb["N"])
    g = (4.0 / 3.0) * 3.0 * cb["prat"] / C ** 2
    if not (os.path.exists(pb) and os.path.exists(pw)):
        rows.append((key, g, None, None, None, "missing B=%d W=%d" % (os.path.exists(pb), os.path.exists(pw)))); continue
    B, Wr = json.load(open(pb)), json.load(open(pw))
    bg, refs = rw.references(cb); wref, vec = refs[rw.ref_key(cb)]
    amp = float(cb["extra"].get("amp", 1e-5))
    dr = np.max(np.abs(np.array(Wr["rho_end"]) - np.array(B["rho_end"]))) / (amp * bg.rho0)
    de = np.max(np.abs(np.array(Wr["e_end"]) - np.array(B["e_end"]))) / (amp * abs(vec[3]))
    dw = abs(complex(Wr["w_re"], Wr["w_im"]) - complex(B["w_re"], B["w_im"])) / abs(wref)
    rows.append((key, g, dr, de, dw, "L1rho B %.3e W %.3e" % (B["l1_rho"], Wr["l1_rho"])))
print("| clos | P | tau | N | g^4 | dev rho | dev E | |dw|/|w| | max dev / g^4 | note |")
print("|---|---|---|---|---|---|---|---|---|---|")
for key, g, dr, de, dw, note in rows:
    if dr is None:
        print("| %s | %.4g | %g | %d | %.2e | - | - | - | - | %s |" % (key + (g ** 4, note))); continue
    m = max(dr, de)
    print("| %s | %.4g | %g | %d | %.2e | %.2e | %.2e | %.2e | %.2g | %s |" % (key + (g ** 4, dr, de, dw, m / g ** 4, note)))
    k2 = (key[0], key[1])
    worst[k2] = max(worst.get(k2, (0, 0)), (m / g ** 4, m))
print("\nper (closure, P): max over tau, N of (dev/g^4, dev)")
for k in sorted(worst):
    print(k, "dev/g^4 %.3g  dev %.3e  g^4 %.3e" % (worst[k][0], worst[k][1], ((4.0 / 3.0) * 3 * k[1] / C ** 2) ** 4))
