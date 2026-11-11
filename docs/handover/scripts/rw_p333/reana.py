#!/usr/bin/env python3
"""rw_p333-1112: analyse finished runs of CASELIST in RUNSDIR (tabs only, results.json not needed).
usage: reana.py CASELIST RUNSDIR  -> N nt cdt/dx w_re w_im l1rho l1E extra"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "rw_amp1e6"))
import run_radwave as rw  # noqa: E402
for s in open(sys.argv[1]):
    if not s.strip() or s.startswith("#"):
        continue
    c = rw.parse_case(s.strip())
    rdir = os.path.join(sys.argv[2], rw.case_tag(c))
    bg, refs = rw.references(c)
    w, vec = refs[rw.ref_key(c)]
    try:
        o = rw.analyse(rdir, c, bg, w, vec, 4, 1)
    except Exception as ex:  # noqa: BLE001
        print("%4d %5d MISSING %s" % (c["N"], c["nt"], ex)); continue
    ex = {k: v for k, v in c["extra"].items() if k not in ("arm", "ts", "vimp", "ref") and
          not k.startswith("rm_implicit_tol") and not k.startswith("rm_implicit_lin")}
    print("%4d %5d %5.3f %.8f %.8f %.4e %.4e t=%.6f %s" % (
        c["N"], c["nt"], 10.0 * 2 * 3.141592653589793 / w.real / c["nt"] * c["N"], o["w_re"],
        o["w_im"], o["l1_rho"], o["l1_E"], o["t_end"], " ".join("%s=%s" % kv for kv in sorted(ex.items()))))
