# NOTE 2026-11-24zzzt orion -> Caltech: t9l62 ray-set test (T9HOT4), P1; run when a node is free, after Compton stops

**Shortened by the coordinator (USER: no unnecessary tests):** 8 cycles per arm (time/nlim = rst cycle + 8) instead of 20, since only cycles 4-8 are read; ~50 GB of dumps per arm. Skip the test if a node isn't free within 1 h, or if the preflight FATALs (post the message). Not on production nodes; delete the full dumps after extraction.


## Why
- On arm B's sub-box the radiation force (median kappa F0/c 1.1e6 cm/s^2 in 23 cells with |v| > 1e8) is the only term large enough for |v| 1-3e8. The budget did not close. This tests one thing: whether that flux changes with the angular set.

## Arms (one test node, 4 GPUs, one after the other; not a production node)
Binary: T2 `athena_ia123rc4v_262affd6f` (md5 7d7c4a71...). Start: `rst_keep/t9l62.00015` (t ~4.75e5, dt 9-13).
Args: those of your Compton arm C0, unchanged (T2 t9l62 production args, every key named), plus the list below.

| arm | `rad_vet/vet_gd_gl_nmu` | `vet_gd_gl_nphi` | directions |
|---|---|---|---|
| G1 | 7 | 18 | 126 (production) |
| G2 | 13 | 36 | 468 |

- **Not 14x36 and not 10x24 (code reading, not run).** With `vet_gd_gl_stagger = true` odd rings are shifted half a step (VgdProduct). For an EVEN nmu a ring and its equatorial mirror have opposite shifts, so `vet_gd_theta_seam_sym = true` should FATAL ("not mirror-symmetric"). 13x36 is the nearest legal set.
- nphi 36 passes the phi-wall rule (2 pi/36 divides pi). 468 is above the default cap of 432: the binary must have been built with `Athena_VET_VGD_NMAX` >= 468 (your fast-build recipe has 512). If it FATALs "too many directions", or is out of memory, use **7x36** (252) and say so.
- Memory estimate: your B log has 17.5 GB of intensities on rank 0 with 126, so ~65 GB with 468 IF the lateral band stays 57 cells. It may grow, and the free memory is not known to us.

Same in both arms (name all of them):
```
rad_vet/vet_gd_gl_stagger=true  rad_vet/vet_gd_ls=0
rad_vet/vet_gd_rotate_dt=0  rad_vet/vet_rotate_avg=0  rad_vet/vet_rotate_avg_ema=0
rad_vet/vet_gd_wall_sym=true  rad_vet/vet_gd_theta_seam_sym=true  rad_vet/vet_gd_wall_interp=true
hydro/dt_report=1  time/nlim=<rst cycle + 8>  job/final_restart=false
output2/dcycle=1  output3/dcycle=1      # hydro_w ; vet (E, F)
output6/dcycle=1  output7/dcycle=1      # vet_vet (D) ; vet_face
output5/dt=1e30                          # no rst
```
Please check in the rst's echo that outputs 2 / 3 / 6 / 7 are hydro_w / vet / vet_vet / vet_face and that 5 is rst. A new output block cannot be added on the command line.

**Preflight:** G2 with `time/nlim` = rst cycle + 1 and no dumps. Stop and post if it FATALs.

## Is the set change legal on this restart? (code reading, not run)
- **Yes for the file layout.** The outer wall is marshak, so no ray intensities are in the rst (`cln_scmem` false): its channel count does not depend on the number of directions.
- **But old-set state is in it:** the extrap1 closure history (moments of the last stage 1), the radiation E and F, and the D used for the first dt. No key drops the history. Step 1's closures use it, and step 2 should be the first free of it. G1 and G2 therefore share the same start, and **cycles 1-3 are a transient: do not read them.**
- Rotation and the ring are off, so there is no xfade or ring state.

## Dumps
- `gid` is NOT safe for `bin` output under MPI here: the non-slice writer still writes every rank's blocks. So the dumps are full: ~23 single-precision variables x 6.7e7 cells = ~6.2 GB per cycle, ~125 GB per arm. If that is too much or slow, use `dcycle=2` on outputs 6 and 7 only.
- After each arm, extract the box and delete the full dumps: gid 0, local k 70-86, j 43-57, i 50-70 (the T9HOT2 box, 17 x 15 x 21) to .npz per cycle: dens, velx/y/z, eint, vet_e, vet_f1-3, the 6 vet_vet, the 5 vet_face, x1v / x2v / x3v, time, dt. Add 2 cells on each side of the box.
- **`vet_face` under staggered clean holds the LAB face momentum flux** (code comments), not F0. Orion forms F0 = F - (E v + v.P) from F, E, D, v. The face F at the thin/dense interfaces is therefore available.
- **Not available as outputs:** opacity (kappa_T), `vet_comov` (FATAL under clean), and per-cell FOFC flags. `vet_fs` is not expected to work with the production band storage (it needs a full intensity array). kappa comes only from DTREP's tauT_h / (rho h) at the binder; Orion uses 0.34 elsewhere.

## Post per arm
1. Per cycle: dt; max and p99 |v| over rho < 1e-11 (box, and whole mesh if cheap); kappa F0/c and F0/(cE), median and p90, over the cells with rho < 1e-11 and |v| > 5e7 in the box; the same for the fixed cell list that G1 has at cycle 4.
2. L_top / L_in, Krylov its (t2s2, t2s3 sums), passes, s/cyc, CW / FATAL / DT_RETRY.
3. Startup lines: "direction set", "lateral halo, band", "ragged per-shell band intensity", and "band-clamped lateral reads" at the end (must be 0 in both arms, or the arms are not comparable).

## Reading (fixed before the run; cycles 4-8, where the gas is nearly the same in both arms)
Compare F0 cell by cell in G1's cells with rho < 1e-11 and |v| > 5e7: the median of |F0(G2)| / |F0(G1)|, and of the cosine between them.
- **Ray effect:** ratio < 0.6, or cosine < 0.8. |v| then falls in G2 only over the ~20 cycles (~200 s, 1-2 cell crossings), so |v| is the weaker test.
- **Not this ray effect:** ratio 0.8-1.25 and cosine > 0.9. This does NOT show the flux is physical: the grid, the lateral interpolation and the closure are the same in both arms.
- **Between:** inconclusive. Say so.

Limits:
- One start, one box, few fast cells (23 over B's 7 dumps).
- No repeat arm, so the noise floor (solver tolerance 1e-5 and trajectory divergence) is unknown. Differences of a few % mean nothing.
- The 1/nphi amplitude scaling is NOT established: ZAMSVPHI measured 1.0e-2 (7x18 staggered) and 2.1e-2 (unstaggered) in F_phi/F_r of a symmetric free-streaming layer; I found no posted 7x36 result. It is a different regime (here F/(cE) is 0.05-0.11 and there is no symmetric reference), so no numerical amplitude is predicted for this test.
- Both nmu and nphi change, so the test cannot say which matters.

Cost (estimate, not measured): B's sweeps were ~0.75 of 3.0 s/cyc; 468 directions put a cycle at ~5-6 s, so compute is ~3 min for both arms, plus two restarts and the dumps. The dumps may exceed the 10 min.

-- orion (coordinator)
