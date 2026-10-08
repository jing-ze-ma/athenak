"""Compare a Delta AG Car smoke run with the viper reference (and Raven numbers for A).
usage: python compare_smoke.py <A|B> <run dir>"""
import glob
import re
import sys

C, D = sys.argv[1], sys.argv[2]
REF = "/work/nvme/bivj/jma20/delta_1008/agcar_files/smoke_ref_viper"
RAVEN = {"A": {1: 1.2525591814773266e+04, 3: 6.1192879615024124e+32}}


def last_row(f):
    rows = [ln.split() for ln in open(f) if not ln.startswith("#") and ln.strip()]
    return [float(x) for x in rows[-1]]


def log_facts(f):
    s = open(f, errors="replace").read()
    out = {}
    for key, pat in [("Tcol", r"IC column T\(rho,eint\)/T_col - 1: max (\S+)"),
                     ("bal", r"he_ic_balance cells: max \|T/T_col - 1\| = (\S+)"),
                     ("pic_mean", r"Picard iterations mean=(\S+)"),
                     ("pic_max", r"Picard iterations mean=\S+ max=(\S+)"),
                     ("nonconv", r"NON-CONVERGED=(\S+)"),
                     ("cpu", r"cpu time used\s+=\s+(\S+)"),
                     ("cycle", r"^time=\S+ cycle=(\d+)")]:
        m = re.search(pat, s, re.M)
        out[key] = m.group(1) if m else None
    out["FATAL"] = s.count("FATAL")
    return out


mine = last_row(glob.glob(f"{D}/*.hydro.hst")[0] if glob.glob(f"{D}/*.hydro.hst")
                else glob.glob(f"{D}/*.hst")[0])
ref = last_row(f"{REF}/smoke{C}.hydro.hst")
print(f"case {C}  run {D}")
for col, name in [(1, "t"), (3, "mass"), (7, "tot-E")]:
    a, b = mine[col - 1], ref[col - 1]
    line = f"  {name:6s} delta {a:.16e}  viper {b:.16e}  rel {abs(a / b - 1):.2e}"
    if C in RAVEN and col in RAVEN[C]:
        line += f"  | raven rel {abs(a / RAVEN[C][col] - 1):.2e}"
    print(line)
fm, fr = log_facts(f"{D}/run.log"), log_facts(f"{REF}/smoke{C}.run.log")
for k in fm:
    print(f"  {k:9s} delta {fm[k]!s:14s} viper {fr[k]!s}")
if fm["cpu"] and fm["cycle"]:
    print(f"  s/cycle  delta {float(fm['cpu']) / int(fm['cycle']):.3f}"
          f"  viper {float(fr['cpu']) / int(fr['cycle']):.3f}")
