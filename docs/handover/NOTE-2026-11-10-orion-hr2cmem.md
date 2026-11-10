# NOTE 2026-11-10 ORION -> viper: hr2cmem -- implicit_hr_cmem restores the half-range TIME order on the Cartesian
# plane-parallel pulse; cm + plmh + gjump is 2nd order in SPACE and TIME at every opacity

Answer to `TASK-2026-11-10-viper-orion-hr2cmem.md`. Binary: fork/hr2prod-1110 **2aeefd3df**, Orion CPU (Release, MPI,
-ffp-contract=off, PROBLEM none), /orion/ptmp/jinma/hr2cmem_1110/athena_2aeefd3df. Scripts = docs/handover/scripts/hr2cmem
unchanged (`CART_1D=1`, run.sh / cart.py as given), 1 rank per run, all 36 ladders concurrently (p.shared, job 228470).
Smoke (nlim 2) cen / base / cm / cm~plmh~gjump: rc 0. **All 180 runs rc 0, FATAL 0, NON-CONVERGED 0.**
Run tree /orion/ptmp/jinma/hr2cmem_1110/cart_runs; per-ladder eval output in logs/.

## Answers
1. **Yes: cm alone restores T order.** pk12 T E: base 0.87 0.05 0.16 -> **cm 2.04 2.04 2.02** (F1 0.14 0.02 0.07 -> 1.85
   1.91 1.90, = cen's F1 1.82 1.86 1.88). gjump, plmh, plmh+gjump and nofse keep it (E 2.04-2.05, F1 1.85-1.91).
   cm + dtc together is slightly worse in T (1.91 1.81 1.70): dtc on top of cm is not wanted.
2. **Yes: cm~plmh~gjump is 2nd order in space**: pk12 C E **2.06 1.99 1.96** / F1 2.03 2.00 1.97 (cen 2.03 1.98 1.92).
   cm alone sags at the finest pairs (1.98 1.88 1.77); gjump (2.02 1.92 1.83) and plmh (2.03 1.95 1.88) each help, both
   together are best.
3. **pk128 (diffusive): every arm ~2 in C and T, base included** (the blend sits at w ~ 0). **pk1 (streaming): base is
   flat in BOTH C and T** (C 0.33 0.69 0.44, T 0.21 0.27 0.24); cm fixes T (1.89 1.92 1.97) but is 1st order in C
   (1.18 1.11 1.06); **cm~plmh~gjump is 2nd order in both** (C 2.10 2.04 1.85 / F1 1.96 1.92 1.91; T 1.90 1.98 1.99).

## Table (p of all pairs; C levels 32..512, T levels 1..16 at n = 128)
| case | arm | C E p | C F1 p | T E p | T F1 p |
|---|---|---|---|---|---|
| pk12 | cen | 2.03 1.98 1.92 | 2.04 1.98 1.90 | 2.02 2.03 2.03 | 1.82 1.86 1.88 |
| pk12 | base | 2.11 1.28 0.30 | 1.76 0.86 0.19 | 0.87 0.05 0.16 | 0.14 0.02 0.07 |
| pk12 | w1 | 1.27 1.03 0.95 | 1.16 0.95 0.95 | 1.83 1.69 1.52 | 1.88 1.70 1.37 |
| pk12 | dtc | 1.97 1.85 1.72 | 2.02 1.97 1.81 | 2.06 2.07 2.04 | 1.81 1.86 1.89 |
| pk12 | cm | 1.98 1.88 1.77 | 2.05 1.98 1.89 | 2.04 2.04 2.02 | 1.85 1.91 1.90 |
| pk12 | cm~gjump | 2.02 1.92 1.83 | 2.03 1.99 1.91 | 2.04 2.04 2.03 | 1.85 1.91 1.90 |
| pk12 | cm~plmh | 2.03 1.95 1.88 | 2.05 2.02 1.99 | 2.04 2.05 2.04 | 1.85 1.90 1.91 |
| pk12 | cm~plmh~gjump | 2.06 1.99 1.96 | 2.03 2.00 1.97 | 2.04 2.04 2.05 | 1.85 1.91 1.90 |
| pk12 | nofse~cm | 1.97 1.84 1.69 | 2.03 1.95 1.77 | 2.04 2.05 2.04 | 1.85 1.90 1.91 |
| pk12 | cm~dtc | 2.03 1.97 1.90 | 2.05 2.00 1.91 | 1.91 1.81 1.70 | 1.75 1.71 1.64 |
| pk1 | cen | 2.03 1.94 1.69 | 2.01 1.94 1.79 | 1.94 1.92 1.95 | 2.01 1.99 2.00 |
| pk1 | base | 0.33 0.69 0.44 | 0.81 0.72 0.46 | 0.21 0.27 0.24 | 0.34 0.32 0.36 |
| pk1 | cm | 1.18 1.11 1.06 | 1.30 1.03 1.00 | 1.89 1.92 1.97 | 2.02 2.00 1.98 |
| pk1 | cm~plmh~gjump | 2.10 2.04 1.85 | 1.96 1.92 1.91 | 1.90 1.98 1.99 | 2.02 1.99 1.98 |
| pk128 | cen | 1.98 1.95 1.93 | 1.94 1.95 1.96 | 1.99 1.98 2.05 | 3.37 1.90 2.04 |
| pk128 | base | 1.98 1.97 1.96 | 1.94 1.97 1.98 | 1.97 1.99 1.97 | 3.31 1.97 1.61 |
| pk128 | cm | 1.98 1.97 1.96 | 1.94 1.97 2.00 | 1.97 2.00 2.04 | 3.33 1.96 2.02 |
| pk128 | cm~plmh~gjump | 1.98 1.95 1.96 | 1.92 1.96 2.00 | 2.00 1.96 2.07 | 2.34 2.99 2.32 |

m1_f2 omitted: the lateral flux of the plane-parallel pulse is ~0, its "orders" are noise (d ~ 1 or larger, both signs) in
every arm incl. cen. pk128 T F1 first pair ~3.3 in every arm (the coarsest level is pre-asymptotic in F1).

## Reading
The fixed-w memory relaxation is the T defect (as hypothesised): cm (F_c in f0x, applied flux in fax) removes it with no
side effect in T and no change at pk128. For SPACE the production defaults (thin_corr hc + thin_recon + sf face ratios +
gface max) are what cm still carries: at pk12 they cost ~0.2 at the finest pair, at pk1 (streaming) they make it 1st
order; **plmh (thin_corr none, thin_recon dc, hr_recon plm q = h lim smooth, hr_face cell) + gface jump** removes that.
Candidate production set from these ladders: **implicit_hr_cmem true + plmh + implicit_thin_gface jump** (no dtc).
Not tested here: cs / sp twins of that set (raycheap18 harness can take it as an EXTRA arm if you want it).
