#!/usr/bin/env python3
"""Standalone driver for the 9 FOFC test files, against an arbitrary AthenaK binary.

The point is to run the SAME cases with the SAME pass criteria as
tst/test_suite/{hydro,mhd}/test_*fofc*_{cpu,mpicpu}.py, but against a binary that
was built elsewhere (e.g. a HIP build) and, on a cluster, through srun.  The test
modules are imported and their test functions called directly, so the criteria are
not copied -- they are the originals.  Only the launcher is replaced:

  testutils.run(input, flags)       -> <launch> <binary> -i <input> <flags>
  testutils.mpi_run(input, flags,n) -> <mpi-launch n> <binary> -i <input> <flags>

Every case prints one line:  PASS/FAIL  <case>  [t=..s]  metrics: ...
The metrics are harvested by wrapping the analysis helpers each test module defines
(max_fofc_count, read_log, check_mass, l1_rel, ...) plus every np.max() the module
evaluates -- which is where the tests' own "measured number" lives.  On a FAIL the
assertion / pytest.fail message is printed verbatim; it carries the measured number
AND the threshold.

Usage:
  python3 run_cases.py --binary /path/to/athena --outdir out_cpu \
      [--launch "srun -n 1"] [--mpi-launch "srun -n {n}"] [--nranks 4] \
      [--only SUBSTR] [--skip-mpi] [--timeout 900] [--results results.txt]
"""

import argparse
import importlib
import os
import shutil
import subprocess
import sys
import time
import traceback

REPO = os.environ.get("ATHENAK_REPO", "/viper/u2/jinma/ATHENAK/athenak")
TSTDIR = os.path.join(REPO, "tst")
VISPY = os.path.join(REPO, "vis", "python")

# ----------------------------------------------------------------- case table
# (module dotted name, test function, kwargs, needs MPI, label suffix)
CASES = []


def _add(mod, fn, kwargs=None, mpi=False):
    kwargs = kwargs or {}
    tag = ",".join(f"{k}={v}" for k, v in kwargs.items())
    label = f"{mod.split('.')[-1]}::{fn}" + (f"[{tag}]" if tag else "")
    CASES.append((mod, fn, kwargs, mpi, label))


H = "test_suite.hydro."
M = "test_suite.mhd."

# 1. hydro blast, 2 shock dirs x 2 EOS
for _sd in ["1", "2"]:
    for _ev in ["ideal", "general"]:
        _add(H + "test_hydro_fofc_blast_cpu", "test_fofc_blast_runs_and_fires",
             {"sd": _sd, "ev": _ev})
# 2. hydro Sod, plm x 3 Riemann solvers
for _fv in ["llf", "hlle", "hllc"]:
    _add(H + "test_hydro_fofc_sod_cpu", "test_fofc_sod_matches_no_fofc",
         {"rv": "plm", "fv": _fv})
# 3/4. spherical-polar and cubed-sphere hydro
_add(H + "test_hydro_fofc_sp_cpu", "test_fofc_sp_cpu")
_add(H + "test_hydro_fofc_cs_cpu", "test_fofc_cs_cpu")
# 5. well-balanced, 4 tests x 2 schemes
for _t in ["test_fofc_wb_fires_and_stays_finite",
           "test_fofc_wb_balanced_state_untouched",
           "test_fofc_wb_conserves_mass",
           "test_fofc_wb_fallback_sees_full_state"]:
    for _wb in ["static", "dynamic"]:
        _add(H + "test_hydro_fofc_wb_cpu", _t, {"wb": _wb})
# 6. hydro blast over MPI + SMR
for _ref in ["none", "static"]:
    _add(H + "test_hydro_fofc_blast_mpicpu", "test_fofc_blast_mpicpu_matches_serial",
         {"ref": _ref}, mpi=True)
# 7. MHD blast: 4 configs x 2 fallback solvers, x2 tests; + mass x4 configs; + vceil
_MHD_CONFIGS = ["cart_x1", "cart_x2", "cart_x3_3d", "spherical_polar"]
for _t in ["test_fofc_fires_and_stays_finite",
           "test_no_flagged_cell_is_bitwise_identical"]:
    for _cfg in _MHD_CONFIGS:
        for _rs in ["llf", "hlle"]:
            _add(M + "test_mhd_fofc_blast_cpu", _t, {"config": _cfg, "rsolver": _rs})
for _cfg in _MHD_CONFIGS:
    _add(M + "test_mhd_fofc_blast_cpu", "test_mass_conservation_matches_fofc_off",
         {"config": _cfg})
_add(M + "test_mhd_fofc_blast_cpu", "test_vceil_clips_and_flags")
# 8. MHD cubed sphere
for _rs in ["llf", "hlle"]:
    _add(M + "test_mhd_fofc_cs_cpu", "test_fofc_cs_fires_and_stays_finite_cpu",
         {"rsolver": _rs})
for _rs in ["llf", "hlle"]:
    _add(M + "test_mhd_fofc_cs_cpu",
         "test_fofc_cs_no_flagged_cell_is_bitwise_identical_cpu", {"rsolver": _rs})
_add(M + "test_mhd_fofc_cs_cpu", "test_fofc_cs_removes_floor_events_cpu")
_add(M + "test_mhd_fofc_cs_cpu", "test_fofc_cs_vceil_clips_and_flags_cpu")
# 9. MHD blast over MPI + SMR
for _ref in ["none", "static"]:
    _add(M + "test_mhd_fofc_blast_mpicpu", "test_mhd_fofc_blast_mpicpu_matches_serial",
         {"ref": _ref}, mpi=True)

# helpers whose return value is a "measured number" worth reporting
_METRIC_HELPERS = ["max_fofc_count", "max_vceil_count", "read_log", "check_mass",
                   "l1_rel"]

METRICS = []


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    if isinstance(v, tuple):
        return "(" + ",".join(_fmt(x) for x in v) + ")"
    return str(v)


class _NpProxy:
    """numpy, but np.max() records the scalars the test measures."""

    def __init__(self, np):
        self._np = np

    def __getattr__(self, name):
        return getattr(self._np, name)

    def max(self, *a, **k):
        out = self._np.max(*a, **k)
        try:
            if self._np.ndim(out) == 0:
                METRICS.append(("max", float(out)))
        except Exception:
            pass
        return out


def instrument(mod):
    """Wrap the module's analysis helpers and its numpy so measurements are logged."""
    import numpy as np
    for name in _METRIC_HELPERS:
        fn = getattr(mod, name, None)
        if fn is None or getattr(fn, "_instrumented", False):
            continue

        def make(name=name, fn=fn):
            def wrapper(*a, **k):
                out = fn(*a, **k)
                METRICS.append((name, out))
                return out
            wrapper._instrumented = True
            return wrapper
        setattr(mod, name, make())
    if hasattr(mod, "np") and not isinstance(mod.np, _NpProxy):
        mod.np = _NpProxy(np)


# ---------------------------------------------------------------- the launcher
def make_launchers(binary, launch, mpi_launch, nranks, timeout):
    import test_suite.testutils as testutils

    def _exec(cmd):
        t0 = time.time()
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, timeout=timeout)
        with open("run_log.txt", "a") as fp:
            fp.write(f"\n$ {' '.join(cmd)}   [{time.time()-t0:.1f}s rc={p.returncode}]\n")
            fp.write(p.stdout[-20000:])
        if p.returncode != 0:
            tail = "\n".join(p.stdout.strip().splitlines()[-25:])
            raise RuntimeError(f"binary failed (rc={p.returncode}): {' '.join(cmd)}\n{tail}")
        return True

    def run(inputfile, flags=[], **kw):
        return _exec(launch.split() + [binary, "-i", inputfile] + list(flags))

    def mpi_run(inputfile, flags=[], threads=None, **kw):
        n = nranks if (threads is None or threads > 1) else 1
        cmd = mpi_launch.format(n=n).split()
        return _exec(cmd + [binary, "-i", inputfile] + list(flags))

    testutils.run = run
    testutils.mpi_run = mpi_run
    return testutils


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--launch", default="", help='serial prefix, e.g. "srun -n 1"')
    ap.add_argument("--mpi-launch", default="mpirun -np {n}",
                    help='MPI prefix with {n}, e.g. "srun -n {n}"')
    ap.add_argument("--nranks", type=int, default=4,
                    help="ranks for the _mpicpu cases (the tests use 4)")
    ap.add_argument("--only", default="", help="substring filter on the case label")
    ap.add_argument("--skip-mpi", action="store_true")
    ap.add_argument("--timeout", type=float, default=1800.)
    ap.add_argument("--results", default="results.txt")
    args = ap.parse_args()

    binary = os.path.abspath(args.binary)
    if not os.path.exists(binary):
        sys.exit(f"no such binary: {binary}")
    outdir = os.path.abspath(args.outdir)
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    # testutils computes its log path relative to ../tst at import time
    os.makedirs(os.path.join(os.path.dirname(outdir), "tst"), exist_ok=True)
    # the tests open "inputs/fofc_*.athinput" relative to the cwd
    os.symlink(os.path.join(TSTDIR, "inputs"), os.path.join(outdir, "inputs"))
    os.chdir(outdir)

    sys.path.insert(0, VISPY)
    sys.path.insert(0, TSTDIR)
    testutils = make_launchers(binary, args.launch, args.mpi_launch,
                               args.nranks, args.timeout)
    testutils.ATHENAK_PATH = REPO

    results_path = os.path.join(outdir, args.results)
    rf = open(results_path, "w")

    def emit(line):
        print(line, flush=True)
        rf.write(line + "\n")
        rf.flush()

    emit(f"# binary   : {binary}")
    emit(f"# outdir   : {outdir}")
    emit(f"# launch   : '{args.launch}'   mpi: '{args.mpi_launch}' n={args.nranks}")
    emit(f"# started  : {time.strftime('%Y-%m-%d %H:%M:%S')}")

    npass = nfail = nskip = 0
    for modname, fnname, kwargs, is_mpi, label in CASES:
        if args.only and args.only not in label:
            continue
        if is_mpi and args.skip_mpi:
            emit(f"SKIP {label:<72} (MPI cases disabled)")
            nskip += 1
            continue
        mod = importlib.import_module(modname)
        if is_mpi and args.nranks != 4:
            mod._nranks = args.nranks
        instrument(mod)
        del METRICS[:]
        t0 = time.time()
        try:
            getattr(mod, fnname)(**kwargs)
            ok, msg = True, ""
        except Exception as exc:                                   # noqa: BLE001
            ok = False
            msg = f"{type(exc).__name__}: {exc}"
            with open("failures.txt", "a") as fp:
                fp.write(f"\n===== {label}\n{traceback.format_exc()}\n")
        dt = time.time() - t0
        seen, metrics = set(), []
        for name, val in METRICS:
            s = f"{name}={_fmt(val)}"
            if s not in seen:
                seen.add(s)
                metrics.append(s)
        mtxt = " ".join(metrics[:12]) if metrics else "(no scalar recorded)"
        if ok:
            npass += 1
            emit(f"PASS {label:<72} t={dt:6.1f}s  {mtxt}")
        else:
            nfail += 1
            emit(f"FAIL {label:<72} t={dt:6.1f}s  {mtxt}")
            for ln in msg.strip().splitlines():
                emit(f"       | {ln}")
        # keep the working directory clean between cases
        testutils.cleanup()

    emit(f"# SUMMARY  pass={npass} fail={nfail} skip={nskip} "
         f"of {npass+nfail+nskip} cases")
    rf.close()
    print(f"\nresults written to {results_path}")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
