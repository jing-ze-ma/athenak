"""analyse existing radwave run dirs (no rerun).  usage: rw_ana.py LIST RUNSDIR"""
import os, sys, json
os.environ["RW_RUNS"] = sys.argv[2]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_radwave as R
from multiprocessing import Pool
def one(line):
    c = R.parse_case(line); tag = R.case_tag(c); rdir = os.path.join(sys.argv[2], tag)
    rp = os.path.join(rdir, "results.json")
    if os.path.exists(rp): return tag + " cached"
    bg, refs = R.references(c); rk = R.ref_key(c); w, vec = refs[rk]
    text, period, dt, icfl, nx2, nx3 = R.make_input(c, bg, w, vec)
    try:
        out = R.analyse(rdir, c, bg, w, vec, nx2, nx3)
    except Exception as ex:
        return tag + " FAIL %s" % ex
    out["case"] = c; out["ref_key"] = rk
    out["refs"] = {k: [refs[k][0].real, refs[k][0].imag] for k in refs}
    tail = open(os.path.join(rdir, "log.txt")).read().splitlines()
    ncyc = [int(s.split("cycle=")[1].split()[0]) for s in tail if s.startswith("time=") and "cycle=" in s]
    out["ncycle"] = ncyc[-1] if ncyc else None
    json.dump(out, open(rp, "w")); return tag + " ok"
lines = [l for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
with Pool(16) as p:
    for r in p.imap_unordered(one, lines): pass
print("ANA DONE")
