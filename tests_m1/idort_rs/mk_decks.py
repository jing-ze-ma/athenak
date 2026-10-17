#!/usr/bin/env python3
"""mk_decks.py -- IDORT radiative-shock family rsx_m<M> (idort-rs-1017).

Builds Cartesian (vet_sc) and spherical-wedge (vet_col + vet_gd) decks for the
Lowrie & Edwards (2008) shocks from
  * the IDORT rs template deck (cart/rs_m2/FIN3_n512/in.athinput of the suite tree) and
    its args.txt = the production Cartesian method keys (FIN3/FIN5.1 set);
  * the sp-only production method keys (vet_col + vet_gd set of the IDORT sp decks,
    e.g. sp/rsw_r1000), physics/problem keys excluded;
  * ROL2U (production since 10-15): implicit_tol_rel=0.1 implicit_tol_rel_max=1e-4
    implicit_opac_lag=pred implicit_opac_lag_passes=2 implicit_gas_newton_u4=true;
  * the reference tables / far-field states written by le08_shock.py --athinput
    (refs/gen_m<M>.log, refs/ref_m<M>.txt).
Every deck starts FROM the semi-analytic profile and holds the shock in its own frame
(Dirichlet far-field states on both x1 faces, marshak radiation baths at the far-field
a T^4); tlim = 3 upstream flow-through times of the domain.  Nothing here is a method
key specific to this test: the method keys are the production sets verbatim.

usage: python mk_decks.py <template_dir> <refs_dir> <out_runs_dir> [P|L]
  case table (Mach, resolutions, cfl) is CASES below.
"""
import os
import re
import sys

ROL2U = ("rad_m1/implicit_tol_rel=0.1 rad_m1/implicit_tol_rel_max=1e-4 "
         "rad_m1/implicit_opac_lag=pred rad_m1/implicit_opac_lag_passes=2 "
         "rad_m1/implicit_gas_newton_u4=true")
# sp-only METHOD keys of the IDORT sp decks (sp/rsw_r1000/*/in.athinput), production set;
# vet_gd_rotate_dt -> auto (the sp decks name a test-time value, 0.3125 code units);
# implicit_precond not named (production default mg, as in the Cartesian deck)
SP_KEYS = """closure=vet_col vet_gd=true vet_gd_replace=true vet_col_lat_taucut=30.0
vet_gd_gl_nmu=10 vet_gd_gl_nphi=12 vet_gd_gl_stagger=true vet_gd_rotate_every=10
vet_col_lat_offdiag=operator implicit_op_stencil=true vet_gd_tangential=false
vet_gd_wall_interp=true vet_gd_twin=true vet_gd_twin_fuse=true vet_gd_twin_det=true
vet_source_noesrc=true vet_scatter=true vet_scatter_form=ma vet_scatter_j=relaxed
implicit_face_opac_n=true implicit_pos_floor=true implicit_pos_floor_solve=true
vet_gd_allow_clamp=true implicit_blend_xthin_mode=all implicit_blend_xthin_wmin=0.0
implicit_blend_xthin_r0=1.5 implicit_hr_recon=dc implicit_hr_recon_fresh=2
implicit_hr_damp=0.0 implicit_hr_pos=kill vet_gd_rotate_align=true vet_gd_rotate_t0=0.0
vet_gd_rotate_stretch=0.02 vet_gd_rotate_slope=correct vet_gd_rotate_dt=auto
vet_gd_twin_state=trans""".split()

# arm P (production): the IDORT rs template names three keys that differ from the
# production decks (AG Car E1 / He giant, b16b66db effective parameters): it sets
# implicit_bc_advect = false (production: true, the code default since m1-bcadv; with
# false the Marshak end faces cannot pass the lab-frame enthalpy flux A E and LE08
# shocks drift / carry an N-independent T error -- rad_m1_implicit.cpp:2744-2752) and
# the explicit-transport keys thick_flux = ap_hll / advect_split = true; the production
# implicit guards/opacity keys below are unset there.  arm L (legacy) = template + the old
# implicit_bc_advect=false (the suite rs decks dropped that line on 10-17 night, USER).
PROD_EXTRA = ("rad_m1/implicit_bc_advect=true rad_m1/thick_flux=none rad_m1/advect_split=false "
              "rad_m1/implicit_opac_update=true rad_m1/implicit_opac_newton=true "
              "rad_m1/implicit_lin_cnorm=1.0e-10 rad_m1/implicit_resid_fatal=1.0e-2").split()
# (Mach, (n_lo, n_hi), cfl list)
CASES = [(1.05, (512, 1024)), (1.2, (512, 1024)), (2, (512, 1024)), (3, (512, 1024)),
         (5, (1024, 2048)), (27, (1024, 2048)), (50, (1024, 2048))]
CFLS = (0.4, 0.2)
# spherical variant: one resolution (n_lo), cfl 0.4, wedge at R0 = RSP x Lx
SP_MACH = (2, 5, 27)
RSP = 1.0e4


def mstr(m):
    return ("%g" % m)


def read_gen(refs, m):
    kv = {}
    for line in open(os.path.join(refs, "gen_m%s.log" % mstr(m))):
        mm = re.match(r"^(\w+) = (\S+)$", line.strip())
        if mm:
            kv[mm.group(1)] = float(mm.group(2))
    return kv


def setval(txt, block, key, val):
    """replace `key = ...` inside <block> (first occurrence)."""
    pat = re.compile(r"(<%s>.*?\n)(\s*%s\s*=\s*)([^\n#]*)" % (re.escape(block), re.escape(key)),
                     re.S)
    out, n = pat.subn(lambda mo: mo.group(1) + mo.group(2) + str(val), txt, count=1)
    if n != 1:
        raise KeyError("%s/%s not in template" % (block, key))
    return out


def addval(txt, block, key, val):
    return txt.replace("<%s>\n" % block, "<%s>\n%-24s = %s\n" % (block, key, val), 1)


def setkeys(txt, kvs):
    """SET each blk/key=val in the deck text (replace if present in the block, else add
    at the top of the block), as ridort's setkeys.py: a FRESH start refuses an absent key
    given on the command line, so the method keys go into the deck itself."""
    add = {}
    for a in kvs:
        k, v = a.split("=", 1)
        b, key = k.split("/", 1)
        add.setdefault(b, {})[key] = v
    out, blk = [], None
    for ln in txt.split("\n"):
        s = ln.strip()
        if s.startswith("<") and s.endswith(">"):
            blk = s[1:-1]
            out.append(ln)
            for key, v in add.get(blk, {}).items():
                out.append("%-24s = %s" % (key, v))
            continue
        if (blk in add and "=" in s and not s.startswith("#")
                and s.split("=")[0].strip() in add[blk]):
            continue
        out.append(ln)
    return "\n".join(out)


def deck(tpl, kv, m, n, cfl, refpath, sp=False):
    t = tpl
    lx = kv["x1max"]
    xs = kv["m1_shock_xs"]
    v0 = kv["m1_shock_v_l"]
    tlim = 3.0 * lx / v0
    t = re.sub(r"problem = grey M1 radiative shock \(Lowrie & Edwards\), M0 = \S+",
               "problem = grey M1 radiative shock (Lowrie & Edwards 2008), M0 = %s, idort-rs-1017 "
               "rsx family (le08_shock.py reference, deck-implied P0/sigma_a)" % mstr(m), t)
    x0 = 0.0
    if sp:
        x0 = RSP * lx
        t = setval(t, "mesh", "x1min", "%.10e" % x0)
        t = setval(t, "mesh", "x2min", "1.1707963267948966")
        t = setval(t, "mesh", "x2max", "1.9707963267948967")
        t = setval(t, "mesh", "ix2_bc", "reflect")
        t = setval(t, "mesh", "ox2_bc", "reflect")
        t = setval(t, "mesh", "nx3", 4)
        t = setval(t, "mesh", "x3min", "0.0")
        t = setval(t, "mesh", "x3max", "1.0471975511965976")
        t = setval(t, "meshblock", "nx3", 4)
        t = addval(t, "mesh", "use_spherical_polar", "true")
        t = re.sub(r"slice_x2\s*=\s*0\.0", "slice_x2    = 1.4507963267948967\nslice_x3    = "
                   "0.10471975511965977", t)
    t = setval(t, "mesh", "nx1", n)
    t = setval(t, "mesh", "x1max", "%.10e" % (x0 + lx))
    t = setval(t, "meshblock", "nx1", n)
    t = setval(t, "time", "cfl_number", cfl)
    t = setval(t, "time", "tlim", "%.6e" % tlim)
    t = setval(t, "rad_m1", "implicit_ebath_x1min", "%.10e" % kv["ebath_l"])
    t = setval(t, "rad_m1", "implicit_ebath_x1max", "%.10e" % kv["ebath_r"])
    for k in ("rho_l", "v_l", "t_l", "rho_r", "v_r", "t_r"):
        t = setval(t, "problem", "m1_shock_" + k, "%.10e" % kv["m1_shock_" + k])
    t = setval(t, "problem", "m1_shock_xs", "%.10e" % (x0 + xs))
    t = setval(t, "problem", "m1_shock_ref", refpath)
    for o in ("output1", "output2"):
        t = setval(t, o, "dt", "%.6e" % (tlim / 3.0 * (1.0 - 1e-9)))
    return t, tlim


def main():
    tdir, refs, out = sys.argv[1:4]
    arm = sys.argv[4] if len(sys.argv) > 4 else "P"
    tpl = open(os.path.join(tdir, "in.athinput")).read()
    cargs = open(os.path.join(tdir, "args.txt")).read().split()
    lst = []
    for m, ns in CASES:
        kv = read_gen(refs, m)
        refpath = os.path.join(refs, "ref_m%s.txt" % mstr(m))
        for sp in (False, True):
            if sp and m not in SP_MACH:
                continue
            for n in ns:
                for cfl in CFLS:
                    if sp and (n != ns[0] or cfl != CFLS[0]):
                        continue
                    geo = "sp" if sp else "cart"
                    d = os.path.join(out, geo, "rsx_m%s" % mstr(m), "%s_n%d_c%g" % (arm, n, cfl))
                    os.makedirs(d, exist_ok=True)
                    txt, tlim = deck(tpl, kv, m, n, cfl, refpath, sp)
                    a = list(cargs)
                    if sp:
                        a = [x for x in a if not x.startswith("rad_m1/closure=")
                             and not x.startswith("rad_m1/vet_tensor=")]
                        a += ["rad_m1/" + k for k in SP_KEYS]
                    a += ROL2U.split()
                    if arm == "P":
                        a += PROD_EXTRA
                    else:   # legacy: the template's old explicit value (removed 10-17 night)
                        a += ["rad_m1/implicit_bc_advect=false"]
                    if sp:
                        txt = txt.replace("vet_tensor               = full\n", "")
                    txt = setkeys(txt, a)
                    open(os.path.join(d, "in.athinput"), "w").write(txt)
                    # the keys are IN the deck; args.txt stays empty (harness compatibility),
                    # keys.txt records what was set on top of the template
                    open(os.path.join(d, "args.txt"), "w").write("\n")
                    open(os.path.join(d, "keys.txt"), "w").write(" ".join(a) + "\n")
                    # cost estimate in zone-cycles: steps = tlim/dt, dt = cfl dx/(v0+a0)
                    lx = kv["x1max"]
                    v0 = kv["m1_shock_v_l"]
                    a0 = v0 / m
                    steps = tlim / (cfl * lx / n / (v0 + a0))
                    nz = n * 4 * (4 if sp else 1)
                    lst.append("%s %d %.3e" % (os.path.relpath(d, out), steps, steps * nz))
    open(os.path.join(out, "list_rsx_%s.txt" % arm), "w").write("\n".join(lst) + "\n")
    print("\n".join(lst))


if __name__ == "__main__":
    main()
