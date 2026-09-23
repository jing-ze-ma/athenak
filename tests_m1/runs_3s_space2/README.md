# runs_3s_space2: second-order enthalpy flux in the implicit M1 operator

Date 2026-09-23, viper, branch `m1-space2` (from rt-integration 51a9adb2). The runs_3r_radwave
patches (`pgen_radwave_eig`, `vet_x1_periodic`) are committed on this branch too.
Build and run tree: `/viper/ptmp2/jinma/space2_3s/` (`base` = 51a9adb2 + those patches,
`new` = this commit; `runs/`, `gate/`, `lists/`).

## Problem

runs_3r_radwave found the implicit operator second order in the complex frequency and
FIRST order in the eigenmode shape (L1 of rho) wherever the enthalpy flux matters: tau = 1e3
at every P_rad/P_gas and P = 100, tau = 10. The cause was the advective enthalpy flux
A_d = a_d E, a_d = v_d + (v.D)_d (Eddington: (4/3) v E). It was donor cell on every face:
velocity AND E both taken from the cell upwind of the face velocity.
- Linearised about v = 0, the error is (4/3) E0 (dx/2) sign(v) dv/dx. It does not
  shrink with the amplitude.
- It is an upwind (one-sided) difference of v, which is anti-diffusive on one
  radiation-acoustic branch (coordinator note, tests_m1/runs_3v_vimplicit).

Every other term of the operator was checked:
- The face flux F0 is centred: w E differences, face-mean kt/theta, v_f g0_f and the
  off-diagonal Eddington divergence.
- The post-solve cell flux F = mean(F0) + a E is cell centred.
- The explicit advective split (`advect_split`, `split_vel = recon`,
  rad_m1_fluxes.cpp) is not used by `transport = implicit` (rad_m1_tasks.cpp replaces
  the whole explicit chain), and it already reconstructs v.

## Change: `<rad_m1>/implicit_enthalpy = upwind | central | plm`

The default is `upwind`, bitwise the old scheme. The key is read only when it is named, so an
input that does not name it keeps its parameter dump. It must be written in the input file,
because a command-line override cannot add a key.

- **The MATRIX is unchanged** in every mode. It keeps the donor-cell a_up E'_up, so it is
  still an M-matrix, E' > 0 is guaranteed by construction, and the 7/19-point stencil,
  od_cache and preconditioners are untouched.
- The difference (high-order face flux - donor-cell face flux), evaluated at the lagged
  Picard iterate E^k, goes to the right-hand side. This is a deferred correction, like
  `implicit_recon = plm_dc`. It is re-evaluated every Picard pass, so at convergence the
  flux is the high-order one, fully implicit in E.
  - The velocity is frozen over the step, as before: the flux stays inside the implicit
    operator and nothing moves to the explicit side.
  - It uses only the solve's dt and the current iterate, never E^n, so it works unchanged
    for a stage solve with dt -> gamma dt and an arbitrary "old" vector.
- **x1 faces**: `rr -= nu cr corr` in the row assembly (`m1_impl_asm`).
- **x2/x3 faces**: `corr` is added to the face fluxes fp/fm/gp/gm of
  ImplicitTransverseTerms. These reach TRHS (the right-hand side) and never the row. With
  bicgstab the right-hand side (KB) is then the high-order flux at E^k minus the low-order
  operator at E^k.
- **Face values** (`M1EnthCorr`, rad_m1_implicit.cpp):
  - `central`: a_f = (a_L + a_R)/2, E_f = (E_L + E_R)/2.
  - `plm`: a_f is the mean of the two van Leer PLM face values of a. At an extremum of a
    this is the central mean; where a is smooth it is a 4-point interpolation. E_f is the
    van Leer PLM value of E from the side upwind of a_f, which is monotone and gives
    E_f <= 2 E_donor, so it cannot produce a negative face value. The face uses central
    where cell L-1 or R+1 is unavailable.
- **Halo**:
  - Multi-D uses the module's ordinary halo, which fills ng layers of EP and a.
  - The 1-D partitioned x1 solve uses part_nlay = min(ng, 2) layers.
  - Periodic x1 wraps inside the block.
  - PLM is used only when BOTH L-1 and R+1 exist. The two MeshBlocks of a shared face
    therefore always make the same choice, and the face flux is single valued.

## Gate (a): radwave x1, all 12 (P, tau), Eddington

- N = 32/64/128; nt = 4096/8192, Richardson in time.
- 216 runs, `lists/x1grid.txt`.
- The upwind arm reproduces runs_3r (L1_64 4.48e-3 ... 5.28e-2).
- `pL1` is the order of L1(rho) from N = 32 to 64 and from 64 to 128. `err_w` is
  |omega - omega_ref|/|omega_ref| at N = 64, against `eddc`.

| P, tau | upwind pL1 | plm pL1 | upwind L1_64 | plm L1_64 | upwind err_w64 | plm err_w64 |
|---|---|---|---|---|---|---|
| 0.1, 0.1 | 2.24 / 2.17 | 2.24 / 2.17 | 3.31e-3 | 3.31e-3 | 7.33e-4 | 7.33e-4 |
| 0.1, 10 | 1.89 / 2.05 | 1.89 / 2.05 | 8.40e-3 | 8.40e-3 | 7.13e-4 | 7.14e-4 |
| 0.1, 1e3 | **1.27 / 1.04** | 2.12 / 2.16 | 4.48e-3 | 2.04e-3 | 7.96e-4 | 9.32e-4 |
| 1, 0.1 | 2.02 / 2.28 | 2.02 / 2.28 | 7.18e-3 | 7.18e-3 | 7.26e-4 | 7.26e-4 |
| 1, 10 | 1.88 / 2.01 | 1.88 / 2.01 | 8.48e-3 | 8.47e-3 | 7.10e-4 | 7.14e-4 |
| 1, 1e3 | **0.98 / 1.01** | 2.23 / 2.08 | 4.78e-3 | 4.21e-4 | 9.51e-4 | 8.82e-4 |
| 10, 0.1 | 1.89 / 2.03 | 1.89 / 2.03 | 8.59e-3 | 8.59e-3 | 7.13e-4 | 7.14e-4 |
| 10, 10 | **1.69 / 1.38** | 1.95 / 2.10 | 7.99e-3 | 6.38e-3 | 7.23e-4 | 7.47e-4 |
| 10, 1e3 | **1.01 / 1.01** | 2.11 / 2.02 | 1.35e-2 | 1.84e-4 | 2.21e-3 | 2.76e-4 |
| 100, 0.1 | 1.88 / 2.01 | 1.88 / 2.01 | 8.56e-3 | 8.55e-3 | 7.10e-4 | 7.14e-4 |
| 100, 10 | **0.85 / 0.94** | 2.13 / 2.11 | 6.80e-3 | 5.72e-4 | 6.40e-4 | 8.16e-4 |
| 100, 1e3 | **1.13 / 1.04** | 2.15 / 2.04 | 5.28e-2 | 8.99e-4 | 2.04e-3 | 4.54e-4 |

- **plm: p_L1 >= 1.88 everywhere, and >= 2.01 from 64 to 128.** The frequency order p_w is
  1.93-2.20 in every case.
- At N = 64 L1(rho) drops by up to 73x (P = 10, tau = 1e3) and 59x (P = 100, tau = 1e3).
- **Frequency.**
  - err_w is 5-8x smaller where the enthalpy flux matters most: P >= 10 at tau = 1e3.
  - It is unchanged to 3 digits in the 7 cases the flux does not touch.
  - It is larger in three: P = 0.1, tau = 1e3 by +17 %; P = 100, tau = 10 by +27 %; and
    P = 10, tau = 10 by +3 %. These are second order with a different constant: at
    N = 128 they are 2.28e-4, 2.03e-4 and 1.91e-4.
- `central` (face-mean a and E) is also second order in L1, but its frequency is worse:
  err_w64 at P = 100, tau = 10 is 1.21e-3, and p_w is 1.58 (32 to 64) at P = 1,
  tau = 1e3. **`plm` is the recommended mode.**
- Picard passes (1-D column solve, tol 1e-11): mean over runs 1.48 plm and 1.63 upwind,
  max 2, NON-CONVERGED 0 in all 216 runs.
- The full tables are in `RESULTS_space.txt` (`python3 space2_table.py`).

## Drift gate: uniform background drift v0 = c_mix (coordinator request)

`<problem>/radwave_v0` (new, read only when named) adds a uniform v0 along the wave vector,
with the matching lab flux (4/3) v0 E.
- The reference is the gas-frame root, Doppler shifted: omega + k v0.
- v0 = the equilibrium mixture sound speed: 11.6 at P = 100 and 3.9 at P = 10. That makes
  v0/c = 1.2e-2 and 3.9e-3.
- N = 32/64/128, nt 4096/8192 with Richardson: 36 runs, `lists/drift.txt`.
- `err_g` is (Im w - Im w_ref)/|w_ref|: positive means too little damping, negative too
  much.

| P, tau | arm | err_c 32/64/128 | err_g 32/64/128 | p_L1 | L1_64 |
|---|---|---|---|---|---|
| 10, 1e3 | upwind | -9.6e-2 / -5.4e-2 / -2.8e-2 | **+2.5e-2 / +6.4e-3 / +1.1e-3** | 0.92 / 0.97 | 9.8e-2 |
| 10, 1e3 | plm | 5.6e-4 / 7.3e-5 / 4.0e-6 | -1.6e-3 / -2.3e-4 / -3.6e-5 | 2.42 / 2.37 | 8.4e-4 |
| 100, 10 | upwind | 2.8e-3 / 8.7e-4 / 2.4e-4 | -1.0e-3 / +3.1e-4 / +2.6e-4 | 1.98 / 1.69 | 4.4e-3 |
| 100, 10 | plm | 2.6e-3 / 7.8e-4 / 2.0e-4 | -2.2e-3 / -2.9e-4 / -3.7e-5 | 1.92 / 1.99 | 4.5e-3 |
| 100, 1e3 | upwind | -9.7e-2 / -5.0e-2 / -2.5e-2 | -5.5e-3 / -9.7e-3 / -6.7e-3 | 0.91 / 0.96 | 1.45e-1 |
| 100, 1e3 | plm | 6.9e-4 / 1.5e-4 / 1.8e-5 | -2.2e-3 / -3.9e-4 / -6.8e-5 | 2.36 / 2.40 | 1.7e-3 |

- **With a drift the donor-cell flux is first order even in the frequency.**
  - The phase speed is 5 % slow at N = 64 for tau = 1e3.
  - At P = 10, tau = 1e3 the damping is 20 % too SMALL at N = 32 (Im w = -4.85 against
    -6.05). No run grew.
  - The drift-free runs of gate (a) are the zero-drift control.
- **plm** converges at 1.9-2.8 in frequency and L1. Its damping error at N = 64 is <= 4e-4
  of |w| in all three cases.
- Picard: mean 1.39 plm and 1.33 upwind, NON-CONVERGED 0.
- The code's M1 model is Galilean invariant only to O(v0/c). The plm errors that keep
  converging at second order down to 4e-6 bound that model error here.

## Gate (b): oblique (xyz) wave, P = 100, tau = 1e3

- N^3 box, `lists/xyz.txt`: nt 1024/2048 with Richardson, MPI.
- N = 32 runs as 4 blocks of 32x16x16 on 4 ranks; N = 64 as 16 blocks of 64x16x16 on 16
  ranks.

| arm | err_w 32 | err_w 64 | p_w | L1 32 | L1 64 | p_L1 |
|---|---|---|---|---|---|---|
| upwind | 8.78e-3 | 2.04e-3 | 2.10 | 1.155e-1 | 5.28e-2 | **1.13** |
| plm | 1.90e-3 | 4.58e-4 | 2.05 | 4.05e-3 | 9.09e-4 | **2.16** |

- The N = 32 plm run on 1 block and 1 rank gives the same omega and L1 to every printed
  digit as the 4-block, 4-rank run (w = 125.57005 - 11.2864 i, L1 4.658e-3, nt = 1024).
  So the x2/x3 halo reads (2 ghost layers) are decomposition independent.

## Gate (c): 2-D He slab and 3-D He box (runs_3p_fastdefault setup), CPU, gcc/openmpi

- Inputs are `gate/inp/*` (slab2d_def + one line). Runs are in `gate/cpu/*`, from
  `cpu_gate.sh` and `cpu_gate2.sh`.
- 2-D slab: 200 s. 3-D box: 84x32x32, 4 blocks, 60 cycles.

**Neutrality, switch off (base = 51a9adb2 + radwave patches, new = this commit).** The
following are bitwise identical (`cmp`):
- 2-D Eddington: hst, bin and rst.
- 2-D vet_sc full: hst and bin.
- 3-D Eddington: hst and bin.
- The radwave `tab` dumps (2 cases).

**Convergence.** In every run: NON-CONVERGED 0, positivity fallbacks 0, and the same min E
as upwind.

| run | Picard mean / max | inner its total |
|---|---|---|
| 2-D Edd upwind | 2.268 / 3 | 25 611 |
| 2-D Edd plm | 2.048 / 3 | 23 185 |
| 2-D Edd central | 2.048 / 3 | 23 217 |
| 2-D Edd plm, 2 blocks, 1 rank / 2 ranks | 2.059 / 3 | 21 644 / 21 645 |
| 2-D vet_sc upwind | 2.675 / 3 | 25 621 |
| 2-D vet_sc plm | 2.614 / 3 | 25 451 |
| 3-D Edd upwind / plm / plm 2 ranks | 2.233 / 2.217 / 2.217, max 4 | 2018 / 1871 / 1867 |
| 3-D vet_sc upwind / plm / plm 2 ranks | 2.483 / 2.467 / 2.467, max 4 | 2051 / 1993 / 1991 |

**Decomposition** (`tools/stats.py`, `tools/cmp.py`):
- plm, 2 blocks or 2 ranks against 1 block: 1e-9 to 1e-12 relative in every column (KE1,
  KE2, F1top, totE, dt). That is the ±1-ulp control level of runs_3p.
- `implicit_op_stencil = false` gives a slab bitwise identical to the default (stencil).

**Physics change, plm against upwind:**
- 2-D slab: F1top/Fin mean 4e-7, KE1 3.6e-3, KE2 5e-3, totE 2.3e-6.
- vet_sc: KE1 2.1e-3, V1max 4e-2.
- 3-D, 9.68 s: KE1 2e-5, KE2 1.6e-3.

## GPU cost (job 11946784, apudev, 1 GPU)

- One binary: md5 dcc72327.
- `HSA_XNACK=1` and `HSA_NO_SCRATCH_RECLAIM=1`.
- Arms interleaved, repeat b in reverse order.
- 3-D He box 84x104x104, 4 blocks, 120 cycles, ms/cycle over cycles 20-120.
- Log: `gpu_time_11946784.log`.

| arm | ms/cycle (a, b) | Picard/step | inner/step | NON-CONV |
|---|---|---|---|---|
| Eddington upwind | 49.4, 49.0 | 2.233 | 41.49 | 0 |
| Eddington plm | **47.3, 47.2** | 2.125 | 38.10 | 0 |
| vet_sc full upwind | 53.8, 52.8 | 2.733 | 39.99 | 0 |
| vet_sc full plm | 52.9, 53.0 | 2.733 | 39.56 | 0 |

- plm is 4 % FASTER with Eddington, because it needs fewer Picard passes and inner
  iterations. It is neutral within noise with vet_sc.
- The extra work is 8 loads and one PLM pair per face in two kernels that already exist.

## Gate (d): vet_sc (vet_tensor = full), x1 wave on an N x 4 mesh, tau 10 / 1e3

- 96 runs, `lists/vetgrid.txt`. The reference is `vetqsc_d`. The vet_x1_periodic patch is
  used.
- `plm` gives p_L1 >= 1.88 (N = 32 to 64) and >= 1.90 (64 to 128), and p_w >= 1.88.
- `upwind` gives 0.85-1.27 at tau = 1e3 and at P = 100, tau = 10.

| P, tau | upwind pL1 | plm pL1 | upwind L1_64 | plm L1_64 |
|---|---|---|---|---|
| 0.1, 1e3 | 1.27 / 1.04 | 2.11 / 2.15 | 4.48e-3 | 2.04e-3 |
| 1, 1e3 | 0.98 / 1.02 | 2.22 / 2.05 | 4.78e-3 | 4.25e-4 |
| 10, 10 | 1.69 / 1.37 | 1.95 / 2.10 | 8.00e-3 | 6.38e-3 |
| 10, 1e3 | 1.01 / 1.01 | 2.02 / 1.90 | 1.35e-2 | 2.09e-4 |
| 100, 10 | 0.85 / 0.94 | 2.13 / 2.11 | 6.80e-3 | 5.72e-4 |
| 100, 1e3 | 1.13 / 1.04 | 2.09 / 1.96 | 5.29e-2 | 9.74e-4 |

The tau = 10 cases with P <= 1 are unchanged, as in Eddington. Picard mean 1.46 plm and
1.67 upwind, BiCGStab 1.02 and 1.06 inner its per step, NON-CONVERGED 0.

## Not done / caveats

- Time is still first order. This change does not touch the time integration.
- The thin-regime over-damping of runs_3r (tau = 0.1) is unchanged. The enthalpy flux does
  not act there, and those cases were already second order.
- No production run uses `plm` yet; it is default off.

## Files

- `run_radwave.py`: copy of runs_3r with new keys:
  - `enth`, which writes implicit_enthalpy into the input file;
  - `v0`, the drift (Doppler-shifted reference);
  - `mb`, the MeshBlock size in x2/x3;
  - `np`, MPI ranks.
- `space2_table.py`: the tables plus Picard and inner statistics, into `RESULTS_space.txt`.
- `tables.py` and `radwave_disp.py` are copies from runs_3r.
- `cpu_gate*.sh`, `cpu_run.sh`, `time.sh` and `lists/*`: exactly what was run.
