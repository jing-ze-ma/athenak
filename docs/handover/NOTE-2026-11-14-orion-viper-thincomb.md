# NOTE 2026-11-14 Orion -> viper: thincomb (base / min / aniso), posted per part: (2) G1 and (3) pulse ladder below; (1) IDORT set follows

**Build:** fork/rsgposlin-1113 = **48558c0dd8ef40718dae9c4ef04b1c07687c8256** (git archive; CMake = the wt_p18 build_none options: mpicxx gcc/13
openmpi/4.1.8, MPI ON, Release, -ffp-contract=off); athena md5 **3a34b876cfefc5b59fd1bffd2c2a6294**; `strings | grep -c implicit_thin_combine` = 2.
No OMP build (the OMP rows are optional). Arms: base (no key) | min (rad_m1/implicit_thin_combine=min) | ani (=aniso), each on top of the deck + EXTRA_HG
(keys absent from a deck are inserted with splitkeys, as for idortp1118). 4 p.exclusive nodes from ~05:30 PST 11-14 (orion034-039).

## (2) Grey atmosphere G1: **PASS as posted (all 3 arms 1.353e-3 <= 2e-2), but the combine key cannot act on this deck**
Deck task-thincomb-1114/grey/atm.athinput run as g1.sh: rad_m1/transport=implicit_x1, time/cfl_number=1e6, implicit_cfl 1 / 1e2 / 1e4, tlim 2000,
output dt 2000; t9_atmosphere.py (from the 48558c0d tree) on the last m1 tab dump. EXTRA_HG added (inserted).
- **asis (deck + EXTRA_HG):** all 9 runs rc=0, FATAL 0, NC 0; max|dE/E| = 1.353e-03, max|df/f| = 1.355e-03, L1(E) 5.008e-04 for base, min and ani
  at every cfl (= the 10-09 central value). The deck is closure m1, implicit_flux central, flux_beam closure, no thin switch: the combine key
  has nothing to act on, so base = min = ani is expected.
- **"Add the thin-switch keys of the IDORT decks": not possible on this deck.** Adding the D6 radiation/thin keys (closure vet_sc, vet_tensor full,
  vet 8x32, implicit_flux blend, flux_beam halfrange, blend idort_f, thin_switch fse h0 0.3 h1 0.6, hcface tau/x1, hcsf j, hr_recon dc, wbsolve,
  face_mmat, face_fct, blend_dtc 1e-6, k1_redist 0) FATALs on interlocks:
  1. with transport = implicit_x1 (1-D): "<rad_m1>/closure = vet_sc needs <rad_m1>/transport = implicit on a multi-D mesh";
  2. made multi-D (mesh/meshblock nx2 = 4, transport = implicit): "<rad_m1>/transport = implicit supports implicit_flux = central only (berthon |
     blend on the x1 faces of the spherical-polar wedge and the cubed sphere)".
  So a Cartesian grey slab cannot carry the fse thin switch. If you want G1 to test the combine rule, it needs an sp-wedge version of the atmosphere
  deck; tell me and I will run it.

```
G1 asis ani cfl 1e2 rc=0 fatal=0 Picard iterations mean=1.369531e+00 max=2.900000e+01 NON-CONVERGED=0.000000e+00 | PASS T9 atmosphere: q0=0.92963248 max|dE/E|=1.353e-03 max|df/f|=1.355e-03  top5: dE=1.353e-03 df=1.355e-03  L1(E)=5.008e-04 (tol 0.01) 
G1 asis ani cfl 1e4 rc=0 fatal=0 Picard iterations mean=3.807692e+00 max=2.500000e+01 NON-CONVERGED=0.000000e+00 | PASS T9 atmosphere: q0=0.92963248 max|dE/E|=1.353e-03 max|df/f|=1.355e-03  top5: dE=1.353e-03 df=1.355e-03  L1(E)=5.008e-04 (tol 0.01) 
G1 asis ani cfl 1 rc=0 fatal=0 Picard iterations mean=1.082988e+00 max=1.500000e+01 NON-CONVERGED=0.000000e+00 | PASS T9 atmosphere: q0=0.92963248 max|dE/E|=1.353e-03 max|df/f|=1.355e-03  top5: dE=1.353e-03 df=1.355e-03  L1(E)=5.008e-04 (tol 0.01) 
G1 asis base cfl 1e2 rc=0 fatal=0 Picard iterations mean=1.369531e+00 max=2.900000e+01 NON-CONVERGED=0.000000e+00 | PASS T9 atmosphere: q0=0.92963248 max|dE/E|=1.353e-03 max|df/f|=1.355e-03  top5: dE=1.353e-03 df=1.355e-03  L1(E)=5.008e-04 (tol 0.01) 
G1 asis base cfl 1e4 rc=0 fatal=0 Picard iterations mean=3.807692e+00 max=2.500000e+01 NON-CONVERGED=0.000000e+00 | PASS T9 atmosphere: q0=0.92963248 max|dE/E|=1.353e-03 max|df/f|=1.355e-03  top5: dE=1.353e-03 df=1.355e-03  L1(E)=5.008e-04 (tol 0.01) 
G1 asis base cfl 1 rc=0 fatal=0 Picard iterations mean=1.082988e+00 max=1.500000e+01 NON-CONVERGED=0.000000e+00 | PASS T9 atmosphere: q0=0.92963248 max|dE/E|=1.353e-03 max|df/f|=1.355e-03  top5: dE=1.353e-03 df=1.355e-03  L1(E)=5.008e-04 (tol 0.01) 
G1 asis min cfl 1e2 rc=0 fatal=0 Picard iterations mean=1.369531e+00 max=2.900000e+01 NON-CONVERGED=0.000000e+00 | PASS T9 atmosphere: q0=0.92963248 max|dE/E|=1.353e-03 max|df/f|=1.355e-03  top5: dE=1.353e-03 df=1.355e-03  L1(E)=5.008e-04 (tol 0.01) 
G1 asis min cfl 1e4 rc=0 fatal=0 Picard iterations mean=3.807692e+00 max=2.500000e+01 NON-CONVERGED=0.000000e+00 | PASS T9 atmosphere: q0=0.92963248 max|dE/E|=1.353e-03 max|df/f|=1.355e-03  top5: dE=1.353e-03 df=1.355e-03  L1(E)=5.008e-04 (tol 0.01) 
G1 asis min cfl 1 rc=0 fatal=0 Picard iterations mean=1.082988e+00 max=1.500000e+01 NON-CONVERGED=0.000000e+00 | PASS T9 atmosphere: q0=0.92963248 max|dE/E|=1.353e-03 max|df/f|=1.355e-03  top5: dE=1.353e-03 df=1.355e-03  L1(E)=5.008e-04 (tol 0.01) 
G1 d6 ani cfl 1e2 rc=1 fatal=1  | FAIL T9 atmosphere: no dump given 
G1 d6 ani cfl 1e4 rc=1 fatal=1  | FAIL T9 atmosphere: no dump given 
G1 d6 ani cfl 1 rc=1 fatal=1  | FAIL T9 atmosphere: no dump given 
G1 d6 base cfl 1e2 rc=1 fatal=1  | FAIL T9 atmosphere: no dump given 
G1 d6 base cfl 1e4 rc=1 fatal=1  | FAIL T9 atmosphere: no dump given 
G1 d6 base cfl 1 rc=1 fatal=1  | FAIL T9 atmosphere: no dump given 
G1 d6 min cfl 1e2 rc=1 fatal=1  | FAIL T9 atmosphere: no dump given 
G1 d6 min cfl 1e4 rc=1 fatal=1  | FAIL T9 atmosphere: no dump given 
G1 d6 min cfl 1 rc=1 fatal=1  | FAIL T9 atmosphere: no dump given 
```

## (3) Stiff pulse_cs time ladder dt = T/1..T/16: **base reproduces viper exactly; min keeps the finest-pair order; aniso loses it on m1_e**
Decks task-thincomb-1114/pulse/pulse_T_{1,2,4,8,16} (1 MeshBlock, serial); arms add only the combine key (min / ani), base adds nothing.
All 15 runs rc=0, FATAL 0, NON-CONVERGED 0, Picard max <= 4. order.py (bin_convert path set to the local tree) on the final so.m1.00001.bin:
```
== base
m1_e   dL1(T1-2,2-4,4-8,8-16) 2.566e-05 1.607e-05 4.503e-06 1.504e-06  orders 0.68 1.84 1.58
m1_f1  dL1(T1-2,2-4,4-8,8-16) 8.947e-06 1.280e-05 9.623e-07 3.278e-07  orders -0.52 3.73 1.55
m1_f2  dL1(T1-2,2-4,4-8,8-16) 1.740e-07 7.530e-08 1.259e-08 3.894e-09  orders 1.21 2.58 1.69
m1_f3  dL1(T1-2,2-4,4-8,8-16) 1.743e-07 7.712e-08 1.223e-08 3.794e-09  orders 1.18 2.66 1.69
== min
m1_e   dL1(T1-2,2-4,4-8,8-16) 2.975e-05 1.502e-05 4.724e-06 1.587e-06  orders 0.99 1.67 1.57
m1_f1  dL1(T1-2,2-4,4-8,8-16) 7.021e-06 5.175e-06 1.005e-06 3.426e-07  orders 0.44 2.36 1.55
m1_f2  dL1(T1-2,2-4,4-8,8-16) 1.046e-07 4.985e-08 1.262e-08 3.898e-09  orders 1.07 1.98 1.69
m1_f3  dL1(T1-2,2-4,4-8,8-16) 1.022e-07 4.845e-08 1.225e-08 3.798e-09  orders 1.08 1.98 1.69
== ani
m1_e   dL1(T1-2,2-4,4-8,8-16) 2.391e-05 1.680e-05 4.417e-06 1.899e-06  orders 0.51 1.93 1.22
m1_f1  dL1(T1-2,2-4,4-8,8-16) 1.811e-05 7.571e-06 1.132e-06 3.575e-07  orders 1.26 2.74 1.66
m1_f2  dL1(T1-2,2-4,4-8,8-16) 2.805e-07 8.835e-08 1.936e-08 5.491e-09  orders 1.67 2.19 1.82
m1_f3  dL1(T1-2,2-4,4-8,8-16) 2.793e-07 8.559e-08 1.882e-08 5.383e-09  orders 1.71 2.18 1.81
```
- **base = viper's reference to all digits** (E dL1 2.566e-5 1.607e-5 4.503e-6 1.504e-6, orders 0.68 / 1.84 / 1.58).
- **FLAGs (order < base - 0.1):**
  - min: m1_e 2->4->8 order 1.67 vs 1.84; m1_f1 2.36 vs 3.73; m1_f2/f3 1.98 vs 2.58/2.66. The finest pair (8->16) equals base on every field
    (1.57, 1.55, 1.69).
  - aniso: m1_e T1-2-4 0.51 vs 0.68, and **8->16 1.22 vs 1.58**; m1_f1 2.74 vs 3.73; m1_f2/f3 2.19/2.18 vs 2.58/2.66.
    At the finest pair aniso is better on the lateral fluxes (1.82 vs 1.69) but worse on E.
- **FLAGs (finest successive difference > 1.05 x base):** min m1_e dL1(8-16) 1.587e-6 (+5.5 %); m1_f1 3.426e-7 (+4.5 %, ok);
  **aniso m1_e 1.899e-6 (+26 %), m1_f1 3.575e-7 (+9 %), m1_f2 5.49e-9 (+41 %)**. min is the closer of the two to base on this ladder.

---
**Correction to the header above:** the 4 nodes started at 07:41 CET = **22:41 PST on 11-13** (sacct), not "~05:30 PST 11-14".

## (2b) sp-wedge G1 (TASK UPDATE, your request): **PASS for base and min at all cfl, aniso PASS at cfl 1 / 1e2 but NON-STEADY at cfl 1e4**
Built by an Orion sub-agent and checked by me against the run outputs. Deck: /orion/ptmp/jinma/thincomb_1114/g1sp/g1sp.athinput (generator
mkdeck.py; arm decks deck_{base,min,aniso}.athinput, with implicit_thin_combine inserted into the deck).
**Deck changes to make it legal on sp** (radiation block = the IDORT d6 sp production block of decks_d6/sp/rsx_m2/P_n512_c0.4; mesh and problem
keys after runs_5e_vetcol2/inp/sph_atm_vc.athinput):
- Mesh:
  - use_spherical_polar, r = 1000..1001 with 128 cells (z = r - r_in; curvature across the shell 1e-3, r^2 F varies by 2e-3);
  - theta = pi/2 +- 0.004 with 4 cells, reflect walls;
  - phi = 0..pi/6 with 4 cells, periodic (vet_gd needs 2 pi/nphi to divide the phi period: nphi 12 -> pi/6);
  - 1 MeshBlock = 1 rank (vet_gd_upwind needs one MeshBlock along x1, and MeshBlocks need >= 4 cells per dimension).
- Problem: m1_test=sph_atm, atm_rho0 = 0.128 e^(1/0.113) = 892.32687 with atm_scale_h 0.113, i.e. the Cartesian rho(z), tau_bot 100.8;
  e_out 1e-10 (vacuum start); atm_temp 1e-12 (cold held gas, see below).
- rad_m1 physics set to the G1 values: c 1, kappa_s 1, kappa_p = kappa_e = kappa_f = 0; implicit_bc_x1min flux 1.0; x1max marshak q 0.5.
- **gas_feedback false and implicit_vimp false** (radiation only, as the Cartesian G1).
- Kept: closure vet_col + vet_gd (GL 10x12) + vet_col_surface_q; implicit_flux blend (faces all), flux_beam halfrange, blend idort_f;
  thin_switch fse (h0 0.3, h1 0.6, hcface tau / x1, hcsf j); hr_recon dc; wbsolve; face_mmat; face_fct; blend_dtc 1e-6; k1_redist 0;
  hesdirk2; implicit_marshak_face linear.
- EXTRA_HG in the deck (time2_stage_local_be false, changed from the production true).
- Two tab radial slices at (theta, phi) cells (1, 1) and (3, 3) to check lateral uniformity.
**FATALs met and fixed:**
- tab output needs 1-D slices -> slice_x2 / slice_x3;
- vet_gd nphi / phi-period rule -> phi width pi/6;
- implicit_vimp needs gas_feedback -> with feedback on, the held pure-scattering gas went eint <= 0 / NaN at cycles 1-4 -> vimp and feedback off;
- held-gas explicit hydro NaN at huge CFL (cycle 535, cfl 1e2; also seen in a minimal deck) -> atm_temp 1e-12 (radiation result unchanged).
Not legal: a 2-D wedge (implicit_lzx needs 3-D) and a 1-D sp mesh (implicit transport needs multi-D).
**Scorer** g1sp/t9_sp.py (copy of t9_atmosphere.py): z = r - r_in; F_loc = F_in (r_in/r)^2; reference = Hopf/Milne with the exact-E_n VET q(tau)
(q(0) 0.577350, q(inf) 0.710446, from runs_5e_vetcol2/milne_v.py), the correct target for vet_col / vet_gd. The Levermore-M1 form stays available
as --ref m1. --hist prints the score of every dump.
Job 229566, p.exclusive orion039, 9 runs at once, 5 h 16 min. tlim 2000 at cfl 1e2 / 1e4 and 500 at cfl 1 (time-accurate; the score is final to
< 1e-7 by t 450, and dE/E between dumps is 1.3e-8 at t 500).

| arm | cfl | rc / FATAL / NC | Picard mean / max | max\|dE/E\| | max\|df/f\| |
|---|---|---|---|---|---|
| base | 1 / 1e2 / 1e4 | 0 / 0 / 0 | 1.015 / 2, 1.101 / 3, 1.902 / 4 | 6.475e-4, 6.475e-4, 6.473e-4 | 6.48e-4 |
| min | 1 / 1e2 / 1e4 | 0 / 0 / 0 | 1.015 / 2, 1.101 / 3, 1.902 / 4 | 5.650e-4, 5.649e-4, 5.647e-4 | 5.65e-4 |
| aniso | 1 / 1e2 | 0 / 0 / 0 | 1.015 / 2, 1.101 / 3 | 2.277e-3, 2.277e-3 | 2.28e-3 |
| aniso | 1e4 | 0 / 0 / 0 | 2.059 / 3 | **8.83e-3 at t 2000 = snapshot, not steady** | 1.11e-2 |

- **The combine key is active on the sp wedge** (the arms are not bitwise). min vs base: max|dE/E| 1.05e-4 at every cfl (F equal to 1e-11).
  aniso vs base: 1.70e-3 at cfl 1 and 1e2; 8.8e-3 at cfl 1e4.
- **Accuracy: min 5.65e-4 < base 6.48e-4 < aniso 2.28e-3** (3.5x base). All converged runs pass 2e-2.
- **aniso at cfl 1e4 never reaches steady state.** Score history: t 547 7.5e-2, 781 5.0e-2, 1016 3.3e-2, 1250 2.8e-2, 1562 3.1e-2, 1797 2.7e-2,
  2000 8.8e-3. E still changes 3e-2 between the last two dumps, flux conservation is off by 2.5e-3 (<= 2e-9 in every other run), and it is the only
  cfl-dependent arm. Count it as FAIL (it fails 2e-2 at most late dumps).
- Lateral spread between the two slices <= 2.5e-10. In every converged run the worst cell is in the top layer (tau 2e-3 to 5e-3).
- **Geometry reference:**
  - the minimal plain sph_atm_vc deck (vet_col + central + surface_q, G1 physics) scores 2.287e-3 with r^2 F conserved to 8e-11 (Cartesian G1
    M1-central: 1.353e-3), so curvature is under control at the 1e-3 level;
  - **possible missing interlock:** the production block with implicit_flux central / flux_beam closure / thin_switch none (plus the two interlock
    fixes faces x1, blend tau_f) converges to a **wrong steady state with no FATAL**: max|dE/E| 3.35 (E ~ +8 too high for tau <~ 1, f_top 0.49).
    A 16-key flip bisection isolates **implicit_face_mmat = true**: flipping only it gives 8.4e-4 (cref runs: 8.384e-4 at cfl 1e2 and 1e4).
    face_mmat with central flux probably should FATAL.

Verdict (2b): **base PASS, min PASS (and best), aniso PASS at cfl 1 / 1e2 but FAIL at cfl 1e4 (non-steady).** Files: /orion/ptmp/jinma/thincomb_1114/g1sp/.
