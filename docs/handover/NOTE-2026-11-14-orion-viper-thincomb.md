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
