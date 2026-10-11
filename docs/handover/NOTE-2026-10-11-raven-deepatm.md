# NOTE deepatm (deepatm-1011): second order in the thick bottom of the atmosphere

Branch `deepatm-1011` (local, NOT pushed) in the worktree /raven/u/jinma/ATHENAK/wt-deepatm, on vetacc-1011 42aa6bbc.
Work dir: /raven/ptmp/jinma/deepatm_1011.

**Commit**
- f0b6b8a4: opt-in `rad_m1/vet_bc_diff = face` (default `cell` = old behaviour, bitwise).

**Binaries** (bin/)
- athena_cpu_base42aa: 42aa6bbc code (built from `git archive 42aa6bbc` in src_base/), md5 06fb6d2b...
- athena_cpu_f0b6b8a4 (= athena_cpu_wip1): md5 1c351a85...

## Verdict

**The cause is the bottom boundary intensity of the vet_sc formal solution.**
- At the x1min face the upward rays start with the diffusion intensity eps = E + 3 F.Omega/c and are carried half a cell to the bottom cell centre.
- The carry already uses the FACE source, S_f = 1.5 S_0 - 0.5 S_1, linear in tau (exact for the linear S of the diffusion regime), and the Milne diagnostic already starts from S_f.
- But the diffusion branch took E of the bottom CELL, not of the face.
  - The isotropic part of the entering intensity is low by E_0 - E_f = (dE/dtau) dtau_cell/2 = 1.5 F dtau_cell. That is O(dx), and it is O(1) when the bottom cell is thick: dtau_cell is 14 at 32 cells and 0.44 at 2048.
  - The sweep carries the deficit about one mean free path upward. That gives a boundary layer in J_sc and H_sc, so f = K/J is wrong there, and E = P/f inherits the error.
- Bottom-cell H_sc/H (exact 1), tJe, levels 128..2048: base 0.958 0.895 0.846 0.854 0.893. With face: **0.984 0.993 0.997 0.999 0.9997**.
- Why it shows as order 1.5:
  - The layer is ~1 tau thick, i.e. ~1/dtau_cell cells, and its amplitude is O(dtau_cell). Its L1 sum relative to the band is O(1/N), so it is first order asymptotically.
  - At 512-2048 cells it is pre-asymptotic: the layer is just being resolved. That gives the ~1.5.
  - The layer shares the same physical depth at all N, so a 10x longer run cannot change it.

**Fix: `vet_bc_diff = face`.**
- E is moved to the face by the source's own extrapolation: E_f = E_0 + (S_f - S_0). This is exactly 1.5 E_0 - 0.5 E_1 for scattering (S = E), and the deep-LTE value for a thermal S.
- It is implemented in all four sweep paths: the single block (VetShortChar), the banded MB sweep (VetRayK), the per-cell MB kernel, and vet_mb_angles (VetSweepAng). The MB paths need no new data, because S_1 is already in the planes.

| atm self-convergence, tau > 10, orders 32..2048 | cell (base 42aa6bbc) | **face** |
|---|---|---|
| tJe (final thinfix arm) | 2.01 2.04 1.77 1.49 1.47 | **1.99 2.00 2.00 2.00 2.00** |
| tEs | 2.01 2.05 1.78 1.50 1.48 | **1.99 2.00 2.01 2.03 2.05** |
| tZb | 2.01 2.04 1.77 1.49 1.47 | **1.99 2.00 2.00 2.00 2.00** |
| cen | 2.06 2.13 1.91 1.57 1.38 | **2.04 2.08 2.20 2.36 2.18** |
| hr (production xthin, dt-dependent, broken at 2048 as before) | 1.26 0.17 -0.42 -0.42 -2.85 | same |
| Kn < 0.03 band, tJe | 2.02 2.05 1.70 1.37 1.39 | **1.99 2.00 2.00 2.01 2.00** |
| whole-domain E self-conv., tJe | 2.00 2.03 1.79 1.53 1.50 | **1.99 2.00 2.00 2.00 2.00** |

- The tau < 0.1, 0.1-1 and 1-10 bands are unchanged to every printed digit for all arms. tJe stays 1.99-2.09 / 2.00 / 2.00.
- **tJe + face is 2nd order in every band, at every level from 32 to 2048.**
- The Hopf L1 at 2048 goes from 3.455e-4 to 3.430e-4. It is still the angular-quadrature floor.

**I propose making `face` the default.** It is clearly a bug:
- The boundary intensity and the half-cell carry referred to different points.
- The Milne branch next to it already used the face.

The default is kept bitwise (`cell`) as asked.

## D. Diagnosis

### D1. Where the first-order part lives: error profile vs depth
- Data: cpu/RESULTS/bottom_profile.txt (ana/bot.py, ana/prof.py).
- Self-error e_n = E_n - R E_2n in the first bottom cells, tJe base:
  - 512: 8.8e-2 9.5e-3 2.9e-3 2.6e-3 ...
  - 1024: 6.6e-2 2.5e-2 8.1e-3 2.9e-3 1.3e-3 8.6e-4 ...
- So the bottom cell does not converge (0.088 -> 0.066), and the error decays upward by about e per dtau = 1. That is a transport boundary layer of ~1 mean free path.
- The interior deep cells converge with ratio 4.5-6.
- At 1024 the bottom cell alone holds 36 % of the whole-domain L1 self-error for tJe (20 % for cen).
- With face, the same cells give 4.2e-3 3.0e-3 2.7e-3 (512) and 1.0e-3 9.2e-4 7.7e-4 (1024): a clean ratio of 4, with no layer.
- The m1_vsc output confirms it at the bottom cell, base cen 512: J_sc = 302.21 against E = 302.46, and H_sc = 0.846 against F = 1.

### D2. Separating the closure from the transport discretisation
- cenE is closure = eddington (fixed 1/3, no formal solution) with the same moment solve and the same flux bottom BC.
  - tau > 10 orders: 2.01 2.03 2.07 2.14 2.31.
  - So the moment discretisation, the flux bottom BC, the face opacities and the 2:1 restriction of the measure are all 2nd order in the deep layer.
  - The first-order term was in the tensor, i.e. in the formal solution.
- cenC (vet_bc_diff = cell given explicitly) is bitwise identical to base cen.

### D3. The other candidates
- **Bottom flux/moment BC.** Ruled out by cenE, and the face fix leaves it untouched.
- **Ghost-cell fill.** Not used by the SC sweep (the boundary is the face intensity).
- **Cell-average vs point source/opacity, T^4.** The atm is pure scattering with a prescribed gas. cenE shows the opacity/face averages are 2nd order.
- **f -> 1/3 round-off / twin noise.** No: the interior deep cells converge at ratio 4.5-6 before the fix as well.
- **Region mask.** The band is fixed by the analytic tau at the coarse centres. The order fix appears in every level pair and in the Kn band (which is defined from the solution), so it is not a mask effect.
- **The order measure.**
  - It is the L1 sum over the band, normalised by sum |E|, and E is largest at the bottom. So yes, the boundary layer dominates it (D1).
  - But the effect is real: the bottom-cell error did not converge.

### D4. Same pattern elsewhere (not changed)
- vet_gd (rad_m1_vetgd.cpp ~2000, `ib = sa + g0*muf` with sa = S of the first shell's CELL, carried with the face source sbt) has the same cell-vs-face mismatch at the inner boundary.
- So does vet_lat (rad_m1_vetlat.cpp ~713).
- The analogous fix there would be sa -> sbt in ib. That affects the He-star / AG Car GPU runs, so it is left for a separate decision.
- vet_col's VC_CORE uses e0 + 3 f0 mu with sbot. It should be checked the same way.

## Validation (base = tJe with base42aa, new = tJeD = tJe + vet_bc_diff = face)

**Default bitwise.** verify/run_verify.sh gives **0 differing data files** in 7 tests:
- ba0 hr, ba0 tJe (bin)
- G1 atm2d hr, G1 atm2d tJe (tab + hst)
- atm tJe 128
- pulse_k0.128 tJe hesdirk2 X128
- rw_t10 tJe be X64

**Multi-block / MPI** (mbchk/, cpu/RESULTS/mbchk.txt; atm cen 256, 4 x1 blocks, implicit_partition = gather):
- 1 rank, 4 ranks, vet_mb_kernel = cell, vet_mb_angles, and vet_mb_tblock = 2 all agree with the single block to 1e-11. That is the implicit tolerance; face against cell differs by 4.4e-4.
- The bottom H_sc is the same on every path: 0.9929 (cell: 0.8946).

**1-D X/T/C ladders** (cpu/RESULTS/order_eval_tJe_base.txt, order_eval_tJeD.txt; pulses k0.128/k12.8/k1280, rw_t10, rw_t1000, be + hesdirk2, P/Pvs T; 122 lines each):
- pulse_k1280 and rw_t1000 are bitwise identical. They have x1-periodic SC with a thick first pass, and the boundary is used only in the first, non-wrapped pass.
- pulse_k0.128 and k12.8 differ by at most 1.2e-6 of the perturbation. Every order is equal to the printed 2 decimals.
- rw_t10 differs at the round-off-drift level of thinfix D2: max dE 4e-11 on E = 300.
  - Raw last orders: be X E 1.70 -> 1.73, hesdirk2 X E 1.69 -> 1.77.
  - Demeaned: be 2.02 -> 2.04, hesdirk2 1.94 -> 1.93.
- 15 pulse `tJePvs` / `tJeDPvs` runs fail the same way on both binaries. That is not a regression: implicit_vimp needs gas coupling, and thinfix ran Pvs only on rw.

**Gates** (cpu/RESULTS/gates.txt):

| | tJe base | tJeD |
|---|---|---|
| G1 cfl 1e2: L1 / NON-CONV | 7.080e-4 / 0 | **6.994e-4** / 0 |
| G1 cfl 1e4: L1 / NON-CONV | 7.081e-4 / 16 | 6.996e-4 / 16 |
| G3 c1 / c10 | 0.0139 / 0.0142 | 0.0139 / 0.0142 (PASS) |
| G5 k1280 / k12800 | 1.016652 / 0.999962 | same (PASS) |
| G5 k12.8 / k128 (fail for every arm) | 0.5079 / 0.5338 | 0.5070 / 0.5457 |

**Beams** (xb20, ba0, ba20, cyl, shd3b, 100 cycles):
- All 5 are stable, with 0 NON-CONVERGED.
- Picard means: 4.86 / 4.19 / 4.88 / 4.02 / 3.00.
- They are **bitwise identical to tJe**: the x1min is a marshak bath (vet_bc_bath), so the diffusion branch is not used.

## Files

**cpu/** is a copy of thinfix cpu/, with these additions:
- ana/od.py: arms `<a>D` (+ face) for cen, tJe, tEs, tZb, hr and the P/Pvs forms (tJeDP, tJeDPvs); cenE (eddington); cenC (cell given explicitly).
- battery_raven.sh: A[tJeD], A[cenD].
- atm_jobs.sh: the atm ladder 32..2048.
- ladders.sh: the 1-D ladders.

**Other**
- ana/bot.py, ana/prof.py: bottom error profile.
- verify/run_verify.sh.
- mbchk/run.sh.
- build_cpu.sh: `sbatch build_cpu.sh <src> <tag> <builddir>`. Note that sbatch on Raven does not export the caller's environment, so pass arguments.

**Run trees**
- cpu/run_base: the 42aa6bbc binary.
- cpu/run: f0b6b8a4.

**Results**
- cpu/RESULTS/{atm_self_regions, atm_hopf, bottom_profile, order_eval_*, gates, verify, mbchk}.txt

**Rules followed**
- CPU only, short jobs.
- The held agc jobs and the other agents' dirs and worktrees were not touched.
- The commit is local only, not pushed.
