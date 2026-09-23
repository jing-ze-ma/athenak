# runs_3v_vimplicit: gas velocity implicit in the M1 enthalpy flux (option 1 of runs_3u)

Phase A (below, first part) is the discrete model and the plan. Phase B (last section) is
the implementation, `<rad_m1>/implicit_vimp`. It is on branch `m1-vimp`, from `m1-space2`
89dc28d7, which has `implicit_enthalpy = plm`. Build and run tree: `/viper/ptmp2/jinma/vimp_3v`.

## Files

| file | content |
|---|---|
| `vimp_model.py` | discrete 1-D Fourier model (staggered E/F0 layout of `rad_m1_implicit.cpp`): `order`, `stab`, and `iters` (a 2-D operator for the BiCGStab count) |
| `vimp_upw_drift.py` | the upwind-cell variant, with a background drift v0 |
| `RESULTS_order.txt`, `RESULTS_stab.txt`, `RESULTS_iters.txt`, `RESULTS_upw_drift.txt` | outputs |

## The model

Staggered grid: E, T, rho and v sit at cell centres, F0 on the faces. Each term follows the code:
- **Face flux:** `F0'_f = theta (F0^n - c^2 a (E_{i+1}-E_i)/(3 dx))`, with Eddington closure.
- **Force:** `dm_i = (dt/c) * 1/2 (k_l F0_l + k_r F0_r)`, the write-back rule.
- **Temperature coupling:** linearised and eliminated.
- **Hydro rows:** centred differences plus Rusanov dissipation at c_s,gas, advanced with Heun.

There are two enthalpy discretisations:
- **`cen`:** `A_f = (4/3) E0 v_f`, where `v_f` is the face mean of the two cells.
- **`upw`:** `A_f = (4/3) E0 v_up`, the velocity of the upwind cell. This is what the code does
  today, with `a = v_cell (1 + w)` taken at the upwind cell.

The model is linearised about v0 = 0. There, the Newton form `a(v^k)E' + (4/3)E^k (v'-v^k)`
equals the fully implicit `(4/3) E0 v'`. The time schemes are:
- **`imp`:** H-ESDIRK2 with v' in the stage solve (option 1).
- **`rhs`:** H-ESDIRK2 with the code's lagged `M1_IW_V1`.
- **`code`:** the present scheme, Heun gas-only followed by backward Euler with v lagged.
- **`codeimp`:** the present scheme with v' in the backward-Euler solve. This is phase B's
  first gate.

### Order

The reference is the semi-discrete acoustic root of the same spatial operator, with 64 cells
per wavelength. The table gives p(1024->2048) over the 12 radwave cases (`RESULTS_order.txt`).

| scheme | cen | upw |
|---|---|---|
| imp | 1.97-2.00 in all 12 | 1.97-2.00 in all 12 |
| rhs | 1.00 at (1,1e3), (10,1e3), (100,10), (100,1e3); 1.21 at (0.1,1e3); 1.56 at (10,10) | same |
| code, codeimp | 0.99-1.05 (backward Euler) | same |

### Boundedness

The table gives max rho(G) - 1 over m = 1..32 (N = 64) and tau = 1..1e6 per unit length
(`RESULTS_stab.txt`), for enthalpy `cen`. A value of 0 means <= 1e-15. Negative values are
round-off or damping.

| cfl | P | imp (Lh=0 / full) | rhs | code | codeimp |
|---|---|---|---|---|---|
| 0.15 | 1 | 0 / <0 | 6.6e-4 / <0 | 0 / 0 | 0 / 0 |
| 0.15 | 10 | 0 / <0 | 2.3e-2 / <0 | 0 / 0 | 0 / 0 |
| 0.15 | 100 | 0 / <0 | 0.23 / 0.15 | 0 / 0 | 0 / 0 |
| 0.30 | 1 | 0 / <0 | 5.4e-3 / <0 | 0 / 0 | 0 / 0 |
| 0.30 | 10 | 0 / <0 | 9.4e-2 / 1.3e-3 | 0 / 0 | 0 / 0 |
| 0.30 | 100 | 0 / <0 | 0.76 / 0.52 | 3.9 / 2.2 | 0 / 0 |
| 0.80 | 100 | 0 / <0 | 9.5 / 2.8 | 47 / 15 | 0 / 0 |

**Finding: the upwind-cell velocity is anti-diffusive.** With `upw`, every scheme grows at
P >= 10, including `imp` and the present `code`. Examples:
- `imp`, full, cfl 0.15: 0.15 at P = 10 and 0.81 at P = 100.
- `code`, full, cfl 0.15: 0.14 at P = 10.

Adding a real drift v0 with its own upwind advection does not remove it
(`RESULTS_upw_drift.txt`). At v0 = c_s, P = 100, cfl 0.15, full, it is still 0.23 (`imp`)
and 0.25 (`code`). At v0 = 0.1 c_s, P = 10 it is 0.11 (`imp`).

The reason is that `(v_i - v_{i-1})/dx` in the E row equals `dv/dx - (dx/2) d2v/dx2`. For one
of the two radiation-acoustic directions this gives a growth rate of about
`(1/6) sqrt(E0) dx k^2`. The gas dissipation at c_s,gas cannot cover it once P_rad >> P_gas.
This behaviour of the present code has only been found in the model. A code gate is proposed below.

Consequence for the plan: the implicit velocity, and preferably the whole velocity in `a`,
must be the **face mean** `v_f`, which the code already computes for the upwind sign. E stays
upwinded. That is `A_f = [v_f + (v_f . D_up)] E'_up`.

### BiCGStab iterations

The model is 2-D, 64 x 64, with a random kappa spread of x10, c = 1e3 and a max-norm tolerance
of 1e-10. The preconditioner is a right x1-line preconditioner (`RESULTS_iters.txt`, `cen`).
Each entry is none / old / tri / penta:
- **none:** today's operator.
- **old:** the new operator with today's tridiagonal preconditioner.
- **tri:** the tridiagonal part of the new operator.
- **penta:** tri plus the x1 +-2 terms.

| a | cfl | P | tau_cell 0.01 | 1 | 100 | 1e4 |
|---|---|---|---|---|---|---|
| g dt | 0.15 | 100 | 118/118.5/121.5/120.5 | 55.5/52.5/54.5/55 | 8/8/8/7.5 | 2.5/5/5/4 |
| g dt | 0.30 | 10 | 126/128.5/127.5/140 | 64.5/66/70/66 | 10.5/10.5/10.5/10.5 | 2.5/3.5/3.5/3.5 |
| g dt | 0.30 | 100 | 124.5/144.5/150/158 | 71/67/67.5/68 | 10.5/11/11/10 | 2.5/8.5/8.5/6.5 |
| dt (BE) | 0.30 | 10 | 154.5/163.5/154.5/157.5 | 100.5/97.5/98/98 | 18/16.5/17/16.5 | 3.5/9.5/9.5/7 |
| dt (BE) | 0.30 | 100 | 165.5/162.5/156/157 | 106.5/91/102/102.5 | 18/33.5/33.5/24 | 3.5/32.5/32.5/22 |

- **Where the count changes:** only in thick cells at P >= 10. In the worst case (backward
  Euler, P = 100, cfl 0.3, thick) it rises by up to 9x, but the absolute count stays at or
  below the thin-cell counts.
- **Tri vs old:** the tri preconditioner is no better than the old one. With uniform
  coefficients the centred Newton term is a pure +-2 wide Laplacian, so its +-1 part is zero.
- **Penta:** it gains at most 1.5x, because the x2 (and x3) +-2 terms stay unpreconditioned.

Recommendation: keep the existing line preconditioner (rbgs/pcr) and do not build a
pentadiagonal PCR.

### Threshold variant

The optional variant applies the coupling only where P_rad/P_gas and the cell optical depth are
large. The homogeneous model cannot rank a spatial threshold. It does show the cost in cells
left lagged:
- **Order:** `rhs` still loses order at P = 0.1, tau = 1e3 (p = 1.21), with 3x the error of
  `imp` (1.89e-5 vs 6.24e-6 at nt = 1024).
- **Iterations:** the threshold saves stencil flops but no iterations. The extra iterations
  come from exactly the cells the threshold would keep (thick, P >= 10).

## Implementation plan (phase B; not started)

Switch: `<rad_m1>/implicit_vimp = false`. When it is off, no new branch touches any
arithmetic, so the result is bitwise. It requires `transport = implicit` and `have_hydro`
with `feedback`.

1. **Face-mean velocity in `a`.** This is in step (b) `m1_impl_lag` (rad_m1_implicit.cpp
   ~4719-4815) and the assembly `m1_impl_asm` (~5167-5260). With vimp on, the enthalpy
   coefficient is built per face from `v_f = 1/2 (v_i + v_ip)` and `(1 + D_dd)` of the
   upwind cell, instead of `M1_IW_ADV` of the upwind cell. This applies on x1 and, in
   `ImplicitTransverseTerms` (~1379-1650), on x2/x3 with `M1_IW_A2/A3`. Store `a_f` on the faces:
   one new iw slot per axis. The write-back `fp1` (~5795) keeps using the cell `a E`, since it
   only reconstructs the cell-centred flux.
2. **The v' Jacobian.** The face flux after the solve (~5470-5545: `fn = th*(f0n - ...)`,
   plus the HLL part and the boundary faces `fb`) is affine in E'. Write it as
   `F0'_f = cL_f E'_L + cR_f E'_R + r_f`. Store (cL, cR) per face and axis: 6 slots. Then:
   - `dv_i = (dt/(c rho_i)) (wl k_l F0'_l + wr k_r F0'_r)`, with the same wl/wr/bmhalf rule as
     the write-back (~5832-5848).
   - The E row of cell i gains
     `(dt/dx_d)(chat/c)[(1 + D_dd,up) E^k_up dv_f(E')]_{i+1/2} - [...]_{i-1/2}`, with
     `dv_f = 1/2 (dv_i + dv_ip)` and `E^k` the Picard iterate.
   - Its constant part (`r_f`, v*) goes to the RHS.
   - Any part of F0' that is not in the operator is evaluated at E^k and put in the RHS, so
     that the converged Picard state is exactly `a(v')E'`. This covers `od` under
     `od_now != operator`, the line-Jacobi transverse terms, the g0 term, and the D_de (e != d)
     cross-velocity terms of a non-Eddington closure.
   - Stencil: +-2 cells per axis, with the diagonal and +-1 parts folded into TA/TB/TC and
     CJM..CKP. With Eddington that is 13-point in 3-D; with vet_sc/tau and the off-diagonal
     operator it is 19 + 6 = 25 slots.
3. **Operator paths.** Add one kernel, `ImplicitVImpOp(xc, yc)`, with `y += (+-2 terms)`,
   called after each of the three ApplyOp branches (stencil, od cache, legacy) and in the
   residual. Later, the fast path can fold it into `ost` by growing it to 25 slots with a new
   `M1StIdx` range 19..24. The +-1 and diagonal parts go into the existing row, so the line
   preconditioner sees them unchanged in structure.
4. **Halo.** E' needs a 2-cell halo. `ImplicitKrylovHalo` already exchanges ng layers; check
   ng >= 2 in `ImplicitInit`. The new face coefficients (cL, cR, rho, k_t) must be in the
   lagged-quantity halo (`M1_NHALO_T/Q`, ImplicitX1Halo / TransverseHalo) so that a block
   boundary's cells i = is-1, is-2 can form dv. Physical boundary faces:
   - Marshak: `F0 = sgn c mq (E' - eb)`, so its E' coefficient goes into the Jacobian.
   - flux, reflect, efix: constant.
   - Periodic x1 (cyclic): no halo and a wrap index; support it only on the od_cache path, as
     for the stencil path now.
5. **Energy and work.** Keep the write-back as it is (dm from the solved F0', work from v and
   v'). Optional step 2, `implicit_vimp_work`: put the linearised work
   `-(chat/c) v_mid^k . dm(E')` (3-point per axis) into the row, and give the gas what the
   solve removed, as step (a) does with SRCR/SRCB, so that conservation holds to round-off.
   The linear model cannot judge this term, because it is second order about v0 = 0.
6. **Positivity.** The Newton term is a positive wide Laplacian with uniform coefficients: its
   +-2 entries are <= 0. Non-uniform E^k and theta create +-1 entries of either sign, so the
   M-matrix property can be lost. The existing positivity fallback (`od_now -> none`) must
   also drop vimp for that step, and count it.
7. **Restart.** No new persistent state. The predictor state (`M1PRED01`) is unchanged, so a
   restart should be bitwise; gate it.
8. **Threshold variant** (`implicit_vimp_pmin`, `implicit_vimp_taumin`, default 0 = all
   cells): zero the Jacobian of a face unless both cells pass. Report it as a variant only.

### Gates, in order

1. **Off-bitwise:** radwave plus the He slab, with the switch off.
2. **Radwave** (runs_3r_radwave patches), backward-Euler scheme, vimp on:
   - no regression at P <= 10;
   - bounded at P = 100, cfl 0.3 (the model gives `codeimp` = 0 against 2.2 for `code`);
   - phase error vs `tables.json`.
3. **Drift instability:** a radwave with a uniform background v0 = c_s at P = 100 and tau = 1e3,
   vimp off vs on. This checks the upwind-cell anti-diffusion found above in the code itself.
4. **He slab 2-D:** NON-CONVERGED = 0; Picard passes and inner iterations vs off.
5. **Restart:** bitwise.
6. **GPU cost:** 3-D He box, 1 GPU on apudev, same binary, interleaved arms with repeats;
   report ms/cycle, passes and inner iterations. The model predicts more inner iterations
   only in thick cells at P >= 10: 2.5 -> 8.5 at g dt, and 3.5 -> 32.5 at dt, cfl 0.3.

## Phase B: `<rad_m1>/implicit_vimp` (implemented, backward-Euler gates)

### What the code does

The code is in `rad_m1_implicit.cpp`: `ImplicitVimpBuild`, `M1VimpRow`, `M1EnthEf`, and the
constants in `rad_m1_implicit.hpp`.

- **The switch.** `implicit_vimp` defaults to false. It is read only when named, so the
  parameter dump does not change, and it must be written in the input file. It needs
  `transport = implicit` on a multi-D mesh with `implicit_solver = bicgstab`, nghost >= 2,
  and coupling, gas_feedback and dbg_gas_force on. Each of these is checked, with a fatal
  error if missing.
- **The Newton form.** Per face of axis d, the enthalpy flux gains
  `E_f^k [(1 + D_dd) dv_d^k + sum_{e!=d} D_de dv_e^k] + E_f^k (1 + D_dd) [P dv_d(E') - P dv_d(E^k)]`.
  - `dv^k` is the velocity change of the iterate's face fluxes, by exactly the write-back
    rule. It uses f0x1/f0x2/f0x3, the bmhalf weights and the fref reference acceleration.
  - `P` is its Jacobian in E'. It is built from the face-normal eliminated flux as step (g)
    forms it: theta, D_dd, the AP-HLL blend and the Marshak faces. The od and g0 terms are
    lagged.
  - At the Picard fixed point the flux is `E_f a_f(v')`, whatever P is. This was measured:
    `implicit_vimp_jscale = 0.5 / 1 / 2` agree to 5e-11 on the P = 100, tau = 1e3 wave.
- **The face value (the plm consistency).** `E_f^k` is the face E that the
  `implicit_enthalpy` mode itself carries, from `M1EnthEf`:
  - `plm`: the van Leer value from the side upwind of the plm a_f;
  - `central`: the mean;
  - `upwind`: the donor cell.

  So the implicit increment rides on the same face E as the plm deferred correction, which
  stays as it was: `a_f(v^n) E_f` at the iterate, minus the donor matrix. The velocity
  increment is the face MEAN of the two cells. Phase A found that the upwind-cell velocity
  is anti-diffusive.
- **Operator.**
  - The diagonal and the x1 +-1 part go into TA/TB/TC, so the line preconditioner (Thomas
    or PCR) sees them and keeps its structure.
  - x1 +-2 and x2/x3 -2..+2 are applied by `M1VimpRow` in all three operator paths
    (stencil, od_cache, legacy), including the fused-reduction ones.
  - `J E^k` and the lagged part go to the right-hand side.
  - The line-Jacobi fallback after a BiCGStab breakdown lags M1VimpRow at E^k. Without this,
    the thin tau = 0.1 waves at N = 128 went into a Picard 2-cycle (4096 NON-CONVERGED); that
    was found and fixed in this gate.
- **Halo.** P and dv^k (12 components) are exchanged once per pass, with their own exchange
  object. `ImplicitHaloCopy/Direct` now offset the component for `nq > 1`; with `nq = 1`
  nothing changes.
- **Positivity.** A non-positive E from the solve drops vimp for the rest of the step. This
  is counted (`implicit_vimp positivity fallbacks`) and never fired in any gate.
- **Energy.** The write-back is unchanged: dm from the solved F0', and work from v and v'.

### Gates (backward Euler, the present time scheme)

| gate | result | source |
|---|---|---|
| switch off, bitwise vs 89dc28d7 (2-D He slab 200 s: default and plm; hst, bin, rst) | BITWISE | gate/cpu e_def_*, e_plm_* |
| radwave x1, 12 (P, tau), N = 64/128, nt = 4096, Nx4 bicgstab, plm | NON-CONVERGED 0 in all 48 runs; omega unchanged where the coupling does not act; where it acts (P >= 10, tau = 1e3) the damping rate grows by the BE error (P=100: -6.61 vs -6.52, P=10: -6.09 vs -6.05) | RESULTS_gates_rw.txt |
| time convergence P=100 tau=1e3 N=32, nt 256..4096 | vimp w_im -8.00, -7.25, -6.88, -6.69, -6.60 -> first order to the same limit as off (-6.51) | tconv |
| stability P=100, tau = 1e3 and 10, cfl_gas 0.3 and 0.15 (nt 24/48), with and without drift v0 = c_s, 10 periods | bounded in every arm, grid noise <= 6e-4 of the amplitude, NC 0. Off is bounded too (the model's grid-scale growth of `code` does not appear in the code) | stab |
| accuracy at cfl_gas 0.3, P=100 tau=1e3 | vimp keeps 0.22 of the wave after 10 periods where the exact decay is 3.6e-3 (off 2.9e-3). This is the BE error of the implicit acoustic pair; the discrete model gives 0.224 / 0.115 for nt 24 / 48, the code 0.224 / 0.115 | stab, vimp_model codeimp |
| drift v0 = c_mix and c_s, P=100 tau=1e3, N 64/128 | NC 0, noise equal to off; L1 up by the BE error (2.2e-3 vs 5.3e-4 at N=128, c_s) | drift |
| He slab 2-D 200 s (plm + vimp): 1 block, 2 blocks, 2 ranks, op_stencil off | NON-CONVERGED 0, Picard 2.048 (= off), inner 9.12 (off 9.12); vs plm: KE1 6.6e-5, KE2 2.2e-5, totE 9.7e-9; stencil off BITWISE; 2 blocks / 2 ranks 1.1e-5 / 1.4e-5 hst (plm-off control 1.2e-5 / 1.3e-5) | gate/cpu e_vim* |
| vet_sc full + vimp slab | NC 0, Picard 2.614 (= off), KE1 3.4e-5 vs plm | v_vim |
| 3-D box 60 cycles, 1 rank / 2 ranks | NC 0, Picard 2.233 (off 2.217), inner 14.0 (off 14.1) | e3_vim* |
| restart (slab vimp, rst at 100 s -> 200 s) | BITWISE bin and rst | e_vim_rst |

The first-order damping of the implicit acoustic pair under BE is expected: `codeimp` in
the model matches the code to three digits. It is the reason H-ESDIRK2 is needed. vimp under
BE is therefore NOT a production option by itself.

### GPU cost (job 11949611, apudev, 1 GPU, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1)

The 3-D He box is 84x104x104 in 4 blocks, run for 120 cycles; ms/cycle is measured over
cycles 20-120. One binary (md5 b5a5e716) was used, with the arms interleaved in three
repeats (a, b reversed, c). Log: `gate/gpu_time_11949611.log`.

| arm | ms/cycle a, b, c | Picard/step | inner/step | NC |
|---|---|---|---|---|
| Eddington plm | 48.6, 48.9, 47.1 | 2.125 | 38.10 | 0 |
| Eddington plm + vimp | 51.7, 52.0, 51.8 | 2.125 | 38.26 | 0 |
| vet_sc full plm | 53.6, 53.6, 54.9 | 2.733 | 39.56 | 0 |
| vet_sc full plm + vimp | 58.2, 58.2, 58.8 | 2.733 | 39.90 | 0 |

vimp costs +7.5 % (Eddington) and +8 % (vet_sc) per cycle. Picard passes are unchanged and
inner iterations rise by 0.4-0.9 %, so the cost is the per-pass build (two kernels plus a
12-component exchange) and the extra reads in every operator application, not the solver.
Folding the ±2 terms into the 19-point stencil (25 slots) is the obvious next saving.
