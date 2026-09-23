#!/usr/bin/env python3
"""Runner + analysis of the linear radiation-modified acoustic wave (runs_3r_radwave;
copied to runs_3s_space2 with the extra key enth=upwind|central|plm ->
<rad_m1>/implicit_enthalpy; copied to runs_3v_vimplicit with vimp=true ->
<rad_m1>/implicit_vimp and nper=<n> -> run n periods; copied to runs_3x_hesdirk2 with
ts=be|hesdirk2 -> <rad_m1>/time_scheme and dbgfail=<cycle> -> time2_dbg_fail).

One CASE = (P_rad/P_gas, tau_lambda, direction, N, closure, transport, nt):

  dir     x1 (1-D mesh, or nx2=4 when the closure needs multi-D), xy (N x N), xyz (N^3)
  N       cells per box side = cells per wavelength projected on each axis
  clos    edd (closure = eddington) | vetf (vet_sc, vet_tensor = full) |
          vetu (vet_sc, uniaxial)
  trans   impl (transport = implicit, the full 7-point operator in multi-D, the x1
          column solve on a 1-D mesh) | expl (transport = explicit, PD-ARS/HLL)
  nt      time steps per period of the reference root (implicit: dt set through
          <rad_m1>/implicit_cfl = c dt/dx; explicit: ignored, dt = cfl_rad dx/c)

Every run starts from the EXACT eigenmode of its reference model
(<problem>/radwave_eig, radwave_disp.py) and runs ONE period P = 2 pi/Re(omega_ref).
Measured: omega from the complex Fourier amplitude of rho at the box fundamental over
17 dumps (phase slope, log-amplitude slope), and the L1 error of rho and E against the
exact damped eigenmode at t = P, normalised by the initial amplitude.

  python3 run_radwave.py run  CASELIST  --np 16     (CASELIST: one case per line)
  python3 run_radwave.py table                      (collect results.json of every case)
"""

import argparse
import glob
import json
import math
import os
import subprocess
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import radwave_disp as rd  # noqa: E402

ROOT = "/viper/ptmp2/jinma/h2_3x"
EXE = os.environ.get("RW_EXE", ROOT + "/new/build_mpi/src/athena")
RUNS = os.environ.get("RW_RUNS", ROOT + "/runs")
BINCONV = "/viper/u2/jinma/ATHENAK/athenak/vis/python"
C_LIGHT = 1.0e3
AMP = 1.0e-5
NDUMP = 16

KDIRS = {"x1": (1, 0, 0), "x2": (0, 1, 0), "xy": (1, 1, 0), "xyz": (1, 1, 1)}


def parse_case(line):
    """'prat tau dir N clos trans nt [key=val ...]' -> dict."""
    f = line.split()
    c = {"prat": float(f[0]), "tau": float(f[1]), "dir": f[2], "N": int(f[3]),
         "clos": f[4], "trans": f[5], "nt": int(f[6]), "extra": {}}
    for kv in f[7:]:
        k, v = kv.split("=", 1)
        c["extra"][k] = v
    return c


def case_tag(c):
    t = "%s_%s_%s_P%g_t%g_N%d_nt%d" % (c["clos"], c["trans"], c["dir"], c["prat"],
                                       c["tau"], c["N"], c["nt"])
    for k in sorted(c["extra"]):
        t += "_%s%s" % (k.split("/")[-1], c["extra"][k])
    return t


def ref_key(c):
    if c["clos"] == "edd":
        return "eddc" if c["trans"] == "impl" else "edd"
    return "vetqsc_d"


def background(c):
    kd = np.array(KDIRS[c["dir"]], dtype=float)
    lam = 1.0 / np.linalg.norm(kd)          # unit box, one wavelength per side
    return rd.Bg(c["prat"], c["tau"], lam=lam, c=C_LIGHT), kd


def drift(c, bg):
    """extra key v0: a uniform drift along the wave vector; 'cmix' = the equilibrium
    mixture sound speed.  0 when not named."""
    v = c["extra"].get("v0")
    if v is None:
        return 0.0
    return bg.mixture_cs() if v == "cmix" else float(v)


def references(c):
    bg, kd = background(c)
    r = rd.all_roots(bg, kdir=kd)
    v0 = drift(c, bg)
    if v0 != 0.0:
        # the roots are those of the gas frame; the lab frequency is Doppler shifted.
        # (The M1 model of the code is Galilean invariant only to O(v0/c); the drift
        # gate therefore also reports the N-extrapolated self-convergence.)
        r = {k: (w + bg.k * v0, v) for k, (w, v) in r.items()}
    return bg, r


def make_input(c, bg, w, vec):
    n = c["N"]
    d = c["dir"]
    vet = c["clos"] in ("vetf", "vetu")
    nx1 = n
    nx2 = n if d in ("xy", "xyz", "x2") else (4 if vet else 1)
    nx3 = n if d == "xyz" else 1
    if d == "x2":
        nx1 = 4
    if d == "x1" and "ny" in c["extra"]:
        nx2 = int(c["extra"]["ny"])       # the x1 wave on an N x ny multi-D mesh
    period = 2.0 * math.pi / w.real
    dxmin = min(1.0 / nx1, 1.0 / nx2 if nx2 > 1 else 1.0, 1.0 / nx3 if nx3 > 1 else 1.0)
    dt = period / c["nt"]
    icfl = C_LIGHT * dt / dxmin
    kap = bg.kappa
    ep = rd.eig_params(bg, vec)
    tr = "implicit" if c["trans"] == "impl" else "explicit"
    clos = "eddington" if c["clos"] == "edd" else "vet_sc"
    lines = []
    A = lines.append
    A("<comment>\nproblem = runs_3r_radwave %s\n" % case_tag(c))
    A("<job>\nbasename = rw\n")
    A("<mesh>\nnghost = 2")
    for ax, nn in ((1, nx1), (2, nx2), (3, nx3)):
        A("nx%d = %d\nx%dmin = 0.0\nx%dmax = 1.0\nix%d_bc = periodic\nox%d_bc = periodic"
          % (ax, nn, ax, ax, ax, ax))
    mb = int(c["extra"].get("mb", "0"))     # MeshBlock size along x2/x3 (0 = whole)
    A("\n<meshblock>\nnx1 = %d\nnx2 = %d\nnx3 = %d\n"
      % (nx1, mb if (mb and nx2 > 1) else nx2, mb if (mb and nx3 > 1) else nx3))
    A("<time>\nevolution = dynamic\nintegrator = rk2\ncfl_number = %s\nnlim = -1"
      % c["extra"].get("cfl", "0.8"))
    nper = int(c["extra"].get("nper", "1"))
    A("tlim = %.17g\nndiag = 1000000\n" % (nper * period))
    A("<hydro>\neos = ideal\nreconstruct = plm\nrsolver = hllc\n"
      "gamma = 1.666666666666667\n")
    A("<rad_m1>")
    A("c_light = %.17g\nchat_over_c = 1.0\ncfl_rad = 0.4\nthick_flux = none"
      % C_LIGHT)
    A("closure = %s\ne_floor = 1.0e-30\nsubcycle = true\nopacity = const" % clos)
    A("kappa_p = %.17g\nkappa_e = %.17g\nkappa_f = %.17g\nkappa_s = 0.0"
      % (kap, kap, kap))
    A("arad = %.17g\ncoupling = true\ngas_feedback = true" % bg.arad)
    A("transport = %s" % tr)
    A("implicit_cfl = %.17g" % (icfl if tr == "implicit" else -1.0))
    A("implicit_tol = %s" % c["extra"].get("tol", "1.0e-11"))
    A("implicit_maxit = 200\nimplicit_lin_tol = %s\nimplicit_lin_maxit = 400"
      % c["extra"].get("ltol", "1.0e-13"))
    A("implicit_solver = bicgstab\nimplicit_offdiag = operator")
    A("implicit_closure_lag = step\nimplicit_flux = central")
    A("implicit_recon = %s" % c["extra"].get("recon", "dc"))
    if "enth" in c["extra"]:
        A("implicit_enthalpy = %s" % c["extra"]["enth"])
    if c["extra"].get("vimp", "false") == "true":
        A("implicit_vimp = true")
    if "ts" in c["extra"]:
        A("time_scheme = %s" % c["extra"]["ts"])
    if "dbgfail" in c["extra"]:
        A("time2_dbg_fail = %s" % c["extra"]["dbgfail"])
    if vet:
        A("vet_nmu = %s\nvet_nphi = %s" % (c["extra"].get("nmu", "4"),
                                           c["extra"].get("nphi", "8")))
        A("vet_tensor = %s" % ("full" if c["clos"] == "vetf" else "uniaxial"))
        A("vet_x1_periodic = true")
        # periodic passes: the slowest ray (mu1 = the largest node) must see tau >= 14
        tau_x1 = kap * 1.0
        npass = int(math.ceil(14.0 * 0.95 / tau_x1)) + 1
        A("vet_x1_npass = %d" % int(c["extra"].get("npass", npass)))
    A("\n<problem>\npgen_name = rad_m1_beam\nm1_test = radwave")
    A("radwave_dir = %s\nradwave_amp = %.17g\nradwave_rho = 1.0\nradwave_t = 1.0"
      % (d, float(c["extra"].get("amp", AMP))))
    A("radwave_eig = true")
    if "v0" in c["extra"]:
        A("radwave_v0 = %.17g" % drift(c, bg))
    for k in sorted(ep):
        A("%s = %.17g" % (k, ep[k]))
    # the bin writer is single precision: double-precision tab dumps of ONE x1 line
    # (x2, x3 = first cell centres).  For the x1, xy and xyz waves the discrete solution
    # depends on i, i+j, i+j+k only, so one x1 line carries every phase.
    sl = ""
    if nx2 > 1:
        sl += "slice_x2 = %.17g\n" % (0.5 / nx2)
    if nx3 > 1:
        sl += "slice_x3 = %.17g\n" % (0.5 / nx3)
    for no, var in ((1, "hydro_w"), (2, "m1")):
        A("\n<output%d>\nfile_type = tab\nvariable = %s\ndata_format = %%.16e\n"
          "dt = %.17g\n%s" % (no, var, nper * period / NDUMP, sl))
    return "\n".join(lines) + "\n", period, dt, icfl, nx2, nx3


def fourier(arr, x1, x2, x3, kvec):
    g3, g2, g1 = np.meshgrid(x3, x2, x1, indexing="ij")
    ph = kvec[0] * g1 + kvec[1] * g2 + kvec[2] * g3
    return complex(((arr - arr.mean()) * np.exp(-1j * ph)).mean()), ph


def line(path, name, nx2, nx3):
    """(time, x1, x2s, x3s, values along the x1 line) of a tab dump."""
    import common  # noqa: E402
    dp = common.load_tab(path)
    x2s = 0.5 / nx2 if nx2 > 1 else 0.0
    x3s = 0.5 / nx3 if nx3 > 1 else 0.0
    return dp.time, np.asarray(dp.x1), x2s, x3s, np.asarray(dp.var(name)).ravel()


def analyse(rdir, c, bg, w, vec, nx2, nx3):
    kd = np.array(KDIRS[c["dir"]], dtype=float)
    kvec = bg.k * kd / np.linalg.norm(kd)
    hp = sorted(glob.glob(os.path.join(rdir, "tab", "rw.hydro_w.*.tab")))
    mp = sorted(glob.glob(os.path.join(rdir, "tab", "rw.m1.*.tab")))
    ts, zs = [], []
    for p in hp:
        t, x1, x2s, x3s, a = line(p, "dens", nx2, nx3)
        ph = kvec[0] * x1 + kvec[1] * x2s + kvec[2] * x3s
        zs.append(complex(((a - a.mean()) * np.exp(-1j * ph)).mean()))
        ts.append(t)
    ts = np.array(ts)
    zs = np.array(zs)
    ph = np.unwrap(np.angle(zs))
    wr = -np.polyfit(ts, ph, 1)[0]
    wi = np.polyfit(ts, np.log(np.abs(zs)), 1)[0]
    # L1 at the last dump against the exact damped eigenmode
    t, x1, x2s, x3s, rho = line(hp[-1], "dens", nx2, nx3)
    te, _, _, _, ee = line(mp[-1], "m1_e", nx2, nx3)
    phase = kvec[0] * x1 + kvec[1] * x2s + kvec[2] * x3s
    ex = np.exp(1j * phase - 1j * w * t)
    AMP = float(c["extra"].get("amp", 1.0e-5))   # noqa: N806
    rho_ex = bg.rho0 * (1.0 + AMP * np.real(ex))
    exe = np.exp(1j * phase - 1j * w * te)
    e_ex = bg.e0 * (1.0 + AMP * np.real(vec[3] / bg.e0 * exe))
    l1r = float(np.mean(np.abs(rho - rho_ex)) / (AMP * bg.rho0))
    l1e = float(np.mean(np.abs(ee - e_ex)) / (AMP * abs(vec[3])))
    return {"t": ts.tolist(), "absZ": np.abs(zs).tolist(), "argZ": ph.tolist(),
            "w_re": float(wr), "w_im": float(wi), "t_end": float(t), "l1_rho": l1r,
            "l1_E": l1e, "x1": x1.tolist(), "rho_end": rho.tolist(),
            "e_end": ee.tolist()}


def run_one(line):
    c = parse_case(line)
    tag = case_tag(c)
    rdir = os.path.join(RUNS, tag)
    res_path = os.path.join(rdir, "results.json")
    if os.path.exists(res_path):
        return tag, "cached"
    os.makedirs(rdir, exist_ok=True)
    bg, refs = references(c)
    rk = ref_key(c)
    w, vec = refs[rk]
    text, period, dt, icfl, nx2, nx3 = make_input(c, bg, w, vec)
    with open(os.path.join(rdir, "rw.athinput"), "w") as fp:
        fp.write(text)
    with open(os.path.join(rdir, "log.txt"), "w") as fp:
        cmd = [EXE, "-i", "rw.athinput", "-d", "."]
        if "np" in c["extra"]:                # MPI ranks (RW_EXE must be an MPI build)
            cmd = ["mpirun", "-np", c["extra"]["np"], "--oversubscribe"] + cmd
        rc = subprocess.call(cmd, cwd=rdir,
                             stdout=fp, stderr=subprocess.STDOUT)
    if rc != 0:
        return tag, "FAILED rc=%d" % rc
    try:
        out = analyse(rdir, c, bg, w, vec, nx2, nx3)
    except Exception as ex:  # noqa: BLE001
        return tag, "ANALYSIS FAILED %s" % ex
    out["case"] = c
    out["ref_key"] = rk
    out["refs"] = {k: [refs[k][0].real, refs[k][0].imag] for k in refs}
    out["period"] = period
    out["dt_target"] = dt
    out["implicit_cfl"] = icfl
    # the Picard/ncycle summary from the log
    with open(os.path.join(rdir, "log.txt")) as fp:
        tail = fp.read().splitlines()
    ncyc = [int(s.split("cycle=")[1].split()[0]) for s in tail
            if s.startswith("time=") and "cycle=" in s]
    out["ncycle"] = ncyc[-1] if ncyc else None
    out["log_tail"] = [s for s in tail if ("cycle" in s.lower() or "picard" in s.lower()
                                           or "nonconv" in s.lower()
                                           or "not converged" in s.lower()
                                           or "residual" in s.lower()
                                           or "bicgstab" in s.lower()
                                           or "vet_sc" in s.lower())][-12:]
    with open(res_path, "w") as fp:
        json.dump(out, fp)
    return tag, "ok w=(%.8g, %.6g) ref=(%.8g, %.6g) l1=%.3e" % (
        out["w_re"], out["w_im"], w.real, w.imag, out["l1_rho"])


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=("run", "input"))
    p.add_argument("caselist")
    p.add_argument("--np", type=int, default=8)
    args = p.parse_args()
    with open(args.caselist) as fp:
        lines = [s.strip() for s in fp if s.strip() and not s.startswith("#")]
    if args.mode == "input":
        for ln in lines:
            c = parse_case(ln)
            bg, refs = references(c)
            w, vec = refs[ref_key(c)]
            print(make_input(c, bg, w, vec)[0])
        return
    with Pool(args.np) as pool:
        for tag, msg in pool.imap_unordered(run_one, lines):
            print("%-60s %s" % (tag, msg), flush=True)


if __name__ == "__main__":
    main()
