---
name: red-giant-envelope-project
description: "The red giant envelope model (src/pgen/red_giant.cpp) started 2026-09-07 on orion at the user's request: what exists, how the stellar opacities were sourced and merged, the conduction-module changes it needed (table_rho axis, flux-limited dt, radial blend flag), the first column results, and what is next"
metadata: 
  node_type: memory
  type: project
  originSessionId: 07f35476-06ea-4aaa-9c53-c127ba307299
  modified: 2026-09-07T18:49:22.769Z
---

**Ask (2026-09-07 evening):** after a discussion of which objects suit the ck-RT +
radiative-conduction + cubed-sphere machinery (brown dwarfs best; Sun/Saturn/M dwarfs
no), the user said "proceed with red giant". Plan agreed: reuse cubed sphere, general
EOS (its 1e6 K ceiling covers an adiabatic giant envelope down to r/R ~ 0.1; radiation
pressure ~1 % at 1e5 K, ~50 % by 3e5 K so eos_radiation MUST be on), well-balanced
scheme, radiative conduction; REPLACE the opacity source and the irradiated two-stream;
ADD point-mass gravity g(r) and a luminosity-injecting inner wall.

## Data: data/stellar_opac (untracked, PROVENANCE.md + .gitignore written)
The LLNL OPAL server RESETS every download from orion (HEAD works, GET dies) -- do not
retry it. Used instead: OPLIB atomic type-1 tables (Farag+2024, Zenodo 15277019, 193 MB
tarball in /orion/ptmp/jinma/Athenak/stellar_opac_raw/) and AESOPUS 2.0 low-T with
grains (Marigo+2022/23, Zenodo 8221362). GS98, X=0.70, Z=0.02. Merged by
`tools/stellar_opac/merge_rosseland.py` (flake8 clean) onto (log T 2.6-8.0 / 0.025,
log rho -14..0 / 0.05): blend 4.0-4.2 in log T; the two agree to median 0.07 dex there,
max 0.12 (OPLIB is ATOMIC -- below log T 3.9 it lacks molecules and is not used).

## Code (COMMITTED + PUSHED 2026-09-07: 35ea57ba conduction, 1d99e22e the extraction, 361be207 red_giant)
- `src/diffusion/conduction.{hpp,cpp}`: `rad_kappa_src = table_rho` (second table axis
  is log rho, all three lookup sites); `rad_blend_radial` (default true; false = the
  tau blend weights only the ANGULAR faces, radial flux-limited diffusion stays on
  everywhere -- the surface treatment for a self-luminous object with no two-stream);
  NewTimeStep now uses the FLUX-LIMITED diffusivity dF/dg = kappa (1+s^2)^-3/2 per
  direction from central-difference gradients. Reason: with the raw kappa the optically
  thin top gave dt = 6e-7 s (hydro dt ~ 200 s). MEASURE the effect (pending).
- `src/pgen/red_giant.cpp` (new, builds): mstar, lstar, teff, ptop, opac_table, mu
  (ideal only), vpert, column_dump, rin (Cartesian only). Point mass; phi = GM(1/rin -
  1/r); IC = hydrostatic column integrated INWARD from the outer wall with dlnT/dlnp =
  min(grad_rad, grad_ad), grad_rad = 3 kappa p L/(16 pi a c G M T^4), T_top^4 =
  Teff^4/2, in a single-thread DEVICE kernel (EOS + table live on device); user BCs
  continue the column into both ghost walls with v_r mirrored; gravity source with
  the WB dynamic option (face-sum form on curvilinear, plain on the column); inner
  flux L/(4 pi rin^2) handed to the conduction module when rad_flux_inner < 0. Works on
  cs / sp (x1 = r) and a 1D Cartesian column (x1 = r - rin).
- inputs/hydro/red_giant_column{,_eos}.athinput: 1.5 Msun, 500 Lsun, Teff 4000 K,
  rin 0.5 R, 0.6 R tall, 512 cells.

## Results so far (/orion/ptmp/jinma/Athenak/red_giant/)
- Column integrates for both EOS. Tabulated: grad_ad dips to 0.12 in the H ionization
  zone, T(0.5 R) = 1.7e5 K, p = 1e11; ideal (mu 1.3 fixed) is 2.3x hotter at the base.
  tau = 2/3 lands at 1.07 R (R nominal from L = 4 pi R^2 sigma Teff^4); convective from
  just below the photosphere down, as a giant should be.
- TRAP: `ptop` must exceed a T_top^4/3 (0.32 dyn/cm2 at 3364 K) or the general EOS
  with radiation has NO gas state and SolveDensity returns the floor -- a column of
  floor density that looks like a run. Now a FATAL in the pgen; ptop = 3 in the inputs.
- Report bug fixed: the RCB diagnostic tested grad_rad >= grad (always true).
- 300-cycle smoke: clean, but dt 6e-7 s -> the flux-limited dt change above.

## Next
1. Rebuild, rerun the smoke tests, confirm dt is hydro-limited and the column stays
   hydrostatic (v_r small vs c_s) with the WB scheme on; energy flux through the column
   in the radiative layers = L/(4 pi r^2).
2. 1D validation against a MESA profile is the honest gate (none on disk; user may have
   MESA -- run/mesa2.mplstyle exists).
3. 3D on the cubed sphere with vpert, then MHD. Commit in stages with tests.

## 21:30 findings: the thin-layer relaxation, and WHY a 1-D column cannot be steady
- dt was 6e-7 s even with the flux-limited diffusivity, because the limiter does
  nothing where T is FLAT (the isothermal top) while kappa_rad diverges there. So the
  tau blend is now REQUIRED by the pgen (rad_tau_lo/hi = 1/10) and, with weight 1-w,
  each thin cell relaxes toward T_eq^4 = 3/4 Teff^4 (tau + 2/3) on t_rad =
  rho c_v/(4 kappa_P rho sigma T^3) by the exact exponential (unconditionally stable).
  dt is now hydro-limited (~180 s); deep envelope hydrostatic to Mach 1e-5; runs end
  cleanly (RedGiantFinal must release EVERY file-scope View, incl. the opacity table).
- The IC follows grad_rad wherever tau < rad_tau_hi (the relaxation's domain).
- THE REAL RESULT: the photospheric layers still cooled 7-40 % and lost mass. Offline
  diagnostic (recomputing the operator's face flux exactly) showed the IC's discrete
  radiative flux is only 8 % of L at tau 20 and < 1 % deeper -- and it is NOT
  resolution (identical at 4x finer). It is grad_ad/grad_rad: the column is
  CONVECTIVE from tau ~ 10 down, radiation carries only that fraction of L, and 1-D
  cannot convect. So the 1-D column can only validate a RADIATIVE envelope:
  `problem/kappa_fac` (scales kappa in the IC AND the operator, ~1e-3) added for that.
- The general-EOS column behaved better at the top (Mach 0.08) once rcv > 0 guarded.

## 21:50 THE 1-D VALIDATION PASSES: constant-opacity grey column (`problem/kappa_const`)
kappa_fac CANNOT make the column radiative: at the photosphere kappa p ~ g tau whatever
the scale, and the convection is driven by kappa RISING inward (kappa p / (g tau) >> 1) --
the physical surface convection zone of every cool star. A CONSTANT kappa has kappa p =
g tau exactly, grad_rad = tau/(4 tau + 8/3) < 0.25 < 0.4, so the ideal-gas column is
radiative throughout and the equilibrium is closed-form. Run (ideal, kappa 0.01, 512
cells, 6000 cycles = 4e5 s ~ 15 atmospheric dynamical times, /orion/ptmp/.../red_giant/grey):
RCB reported "-1" (none); mass/energy drift 3e-5; deep Mach 1e-4..1e-7; top Mach 0.06;
max |T/T_ic - 1| 5 %, |rho/rho_ic - 1| 14 % (thin layers only). Against the SPHERICAL
Eddington profile T^4 = 3/4 (L/4 pi r^2 sigma)(tau + 2/3): max err 4.4 % (tau < 10),
2.3 % (10 < tau < 1e3), 6.4 % deep (mean -3.5 %, the Eddington closure vs diffusion).
=> gravity + WB source, inner flux, table plumbing, flux-limited dt, blend + relaxation,
and both walls are VERIFIED in 1-D. Remaining few-% offsets sit in the thin layers.

## STAGE 3 (2026-09-07 late): red_giant now RUNS THE BAND SOLVER
After [[correlated-k-shared-module]] made the two-stream reusable, `red_giant.cpp` was
wired to it: `problem/rt_ck = true` replaces the grey Eddington relaxation with the
correlated-k two-stream in the optically thin layers, with the tau blend handing over to
radiative diffusion underneath exactly as on the hot Jupiter.

Setup in UserProblem: reads the ck tables, runs both self-tests, and sets
`two_stream_rt::rt_int_at_cut = false` (the inner WALL carries L, not the ck cut),
`rt_tint_override = teff` (NEW knob in the module: a self-luminous object cannot use the
Thorngren T_int(T_eq) fit, which is for irradiated planets and floors at 100 K),
`rt_star_teff = 0`, and `hot_jupiter_param` with **Teq = 0**, which zeroes the stellar
sweep. Default is OFF, so nothing else changes.

**TRAP that cost a segfault:** `UserProblem` runs inside the ProblemGenerator
CONSTRUCTOR, and `main.cpp` assigns `pmesh->pgen` only when that constructor RETURNS. So
`pmy_mesh_->pgen->anything` is a null dereference in UserProblem. Fill `hot_jupiter_param`
as a plain member; only the source/BC functions (which run later) may go through
`pm->pgen->`.

**First comparison, cs 6x16x16x128, 200 cycles from the same IC:**

| | max abs rho/rho0 - 1 | max abs e/e0 - 1 | dt range |
| --- | --- | --- | --- |
| grey relaxation | 74 | 58 | 72-765 s |
| correlated-k | 0.73 | 1.00 | 4.5-765 s |

So the BAND solver is enormously better behaved in the thin layers -- the grey relaxation
was letting the top run away by ~2 orders of magnitude -- at the cost of a smaller
minimum dt (4.5 s) and ~30 % more wall time per cycle. Neither is settled yet: both are
still transients at 3e4 s, and no convergence or flux-conservation gate has been run.
NEXT: gate the emergent flux against sigma Teff^4 and the T(tau) structure, then run long.

### The FLUX GATE on the band solver (cs, first RT call, problem/ck_dump_file)
red_giant now forwards ck_dump_file/m/j/k to the solver, so the dhj column dump works here
too. On the initial state, one column (mu0 = -1, i.e. the anti-substellar direction, which
is meaningless without irradiation):

    top face  tau_R = 0     F_lw = 0.784 x L/(4 pi r^2)   (0.664 x sigma T_int^4)
    tau ~ 5   w_diff = 0.79 -> the blend hands over to diffusion, as designed
    tau > 200 w_diff = 1    -> pure diffusion, F_lw = 0

So the band solver carries **78 % of the required flux at t = 0**, and the missing 22 % is
NOT obviously a bug: the initial column was integrated with grad_rad from the ROSSELAND
mean, and it lands on the grey Eddington value at tau = 2/3 (T = 4001 K against
Teff = 4000). Handing that grey-consistent profile to a BAND solver and getting 22 % less
flux is the size of the line-blanketing effect that a grey mean cannot represent -- which
is the whole argument for using the band solver here. The gate to run next is whether the
ratio relaxes toward 1 as the atmosphere adjusts; the dump is one-shot (rt_dump_done), so
that needs either a re-dump hook or a global energy-balance check against L.

### THE OPACITIES ARE NOT CONSISTENT (user asked; measured 2026-09-07 late)
red_giant currently uses TWO unrelated opacity datasets either side of the tau blend: the
merged STELLAR table (OPLIB+AESOPUS) for the diffusion, and the Exo-FMS CORRELATED-K
tables for the two-stream. For the hot Jupiter this cannot happen -- `ck_build_rosseland_
table` derives the diffusion's mean FROM the k-table, which is the whole point of that
routine. Here I bypassed it.

New diagnostic `problem/opac_compare` writes the ck-derived Rosseland mean on its own
(T,p) grid, built into the Conduction object and then undone. Comparing against the
stellar table at the GIANT's own conditions (T 2500-6100 K, logR in [-8,1], rho from the
run's own EOS table), over 421 nodes:

    median kappa_ck / kappa_stellar = 0.70   (10-90 %: 0.46 .. 0.86, max 1.1)

So they disagree by ~30 % typically and up to ~2x. The blend therefore hands over between
two operators that disagree about the opacity by that much -- and, worse, about the
OPTICAL DEPTH SCALE, so tau_lo/tau_hi sit at a different physical depth than the RT thinks.

**Also found: the merge script CLAMPS out of range instead of refusing.** Both sources are
tabulated in logR = log rho - 3 log T + 18 over [-8, +1]. At hot-Jupiter conditions
(1000 K, 1-21 bar) logR is 4-6, far outside, and the clamp returns a value CONSTANT in
density over three decades -- silently wrong. In the giant column 2.1 % of points are
outside, all in the isothermal top above the photosphere (logR -9.4 .. -8.0), which is
exactly where the two-stream is meant to take over.

**Recommended fix (NOT done):** build both means and splice them -- ck-derived below
~5000 K where the k-table is valid and the molecular/atomic opacity is what the RT sees,
stellar above (the ck table stops at 6100 K while the envelope reaches 1e5-1e6 K), with a
blend in log T, exactly as merge_rosseland.py already does for OPLIB/AESOPUS. And make
that script REFUSE or flag out-of-range logR rather than clamp.

### WHY the two opacities differ -- and the splice REJECTED (2026-09-07, user asked)
Resolving the ck/stellar Rosseland ratio in TEMPERATURE separates two causes:

    T < 1500 K      ratio 0.00-0.02   GRAINS. The stellar low-T table is AESOPUS 2.0
                                      WITH SOLID GRAINS; the ck tables condense species
                                      out of the gas and never put the dust opacity back.
    1500-2500 K     0.02 -> 0.36      the condensation transition
    T > 2700 K      flat 0.6-0.74     NOT grains. Mostly METALLICITY: my table was
                                      Z = 0.020, the ck tables are 1x solar ~ 0.014.
                                      Rebuilding at Z = 0.014 moves the plateau to ~0.84.

The user's guess (no dust in ck) is CONFIRMED and is worth 1-3 ORDERS OF MAGNITUDE below
1500 K. It does not touch the red giant (its coldest cell is 3364 K) but it is decisive
for brown dwarfs (photosphere 1000-2000 K) and Saturn -- a number for the "clouds are the
deciding piece of work" claim in the earlier survey.

**The proposed fix FAILED its gate and is retracted.** Splicing the ck-derived mean into
the diffusion table (tools/stellar_opac/merge_ck_stellar.py, kept as a diagnostic) gives,
along the giant's own column in the shared layers, spliced/stellar median 1.00 but RANGE
0.19..4.90 -- factors of 3-5 either way across the photosphere. That manufactures
discontinuities rather than consistency, so the table was deleted.

**New reading:** at a giant photosphere (rho ~ 1e-9, T 3400-5000 K) the exoplanet ck
tables are outside their design regime, where stellar opacity is H- and metal lines. For a
STAR trust the stellar tables; a band solver for a red giant needs band opacities built
from stellar data. That also reopens the earlier "the band solver carries 78 % of the
flux" result -- the deficit may be the ck opacities being wrong here, not the grey IC.

### RETRACTION (same evening, after the user pushed back): the tables AGREE at the photosphere
The user asked "can you confirm... I thought the ck table also takes into account the
hydrogen ion". They were right and my factor-3-5 claim was WRONG -- three bugs on my side:
(1) `column_dump` writes p in dyn/cm^2, I read it as bar (1e6); (2) the ck T grid is
NON-UNIFORM (0.30 dex spacing at 100 K, 0.015 at 6100 K) and I indexed it as uniform,
wrong by up to 40x at the hot end -- the code's `RosselandTable()` bisects; (3) points
above the ck ceiling of 6100 K clamp, and I read the clamped values as physics.

With both lookups VALIDATED (ck reproduces the eight start-up `kappa_R/Freedman` values to
1.0000; stellar reproduces the pgen's own kappa column to 1.000), along the giant column
where the ck table is valid (T 3364-5995 K, tau 1e-4..6):

    kappa_ck / kappa_stellar(Z=0.014):  median 0.95, 10-90 % 0.84..1.01

Both carry H- bf/ff, the dominant continuum there, so agreement is what SHOULD happen.
**So there is no opacity inconsistency to fix at the handover** (once the metallicity is
matched: use the Z=0.014 table with the ck tables, not Z=0.020). The splice is
unnecessary; merge_ck_stellar.py stays as a diagnostic only.

The GRAIN result stands -- it came from a node-based comparison with no interpolation, so
bug 2 does not touch it: ck/stellar is 0.00-0.02 below 1500 K. Irrelevant to the giant,
decisive for brown dwarfs.

Consequence: the earlier "the band solver carries only 78 % of the required flux" is back
to being what I first said -- a grey-built IC handed to a band solver -- NOT an opacity
error. LESSON: validate a table lookup against numbers the code itself prints BEFORE
drawing any conclusion from it.

### READINESS FOR A 3-D PRODUCTION RUN: NOT YET. The ck source is radiatively STIFF.
Ran the full combination (cs 6x16x16x128, hydro + general EOS + etotgrav + WB polytropic
+ radiative conduction with the tau blend + correlated-k two-stream) for 1100 cycles:

    dt:  766 s at cycle 0  ->  12.5 (c100)  ->  2.70 (c200)  ->  1.85 (c1100), STILL falling
    t reached 2.1e4 s after 1100 cycles; a convective turnover is R/v_conv ~ 2.9e6 s
    the RT source limiter FIRED (rt_de_max clipping: the radiative time is below dt)

Same setup with the GREY relaxation instead of ck holds **dt 72-765 s**. The difference is
not the physics but the DISCRETISATION: my grey relaxation is an exact exponential toward
the equilibrium temperature and so is unconditionally stable, while the ck two-stream
source is operator-split and EXPLICIT, and merely clipped when it goes out of range.

So at ~1.9 s/step a single turnover needs ~1.1e6 steps ~ 73 h on 8 threads, and the trend
is still downward. **A production 3-D run is not viable in this configuration.**

Fixes, cheapest first:
 1. move the stiff layers out of the domain (lower x1max, or raise the floors);
 2. use the ck solver only below some tau and the IMPLICIT relaxation above -- the code is
    already structured for this (the blend weight w is available per face);
 3. make the ck source implicit: relax toward the band equilibrium with the same
    exponential. That is the real fix and is moderate work in two_stream_rt.hpp.
Also still open from before: the emergent flux is 0.78 of L at t=0 and has not been shown
to relax to 1.
