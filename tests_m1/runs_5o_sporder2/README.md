# runs_5o_sporder2: second order in space for implicit M1 + VET on the sp wedge

- **Date.** 2026-09-25, viper.
- **Branch.** `m1-sp-order2` from rt-integration 7aabf6fd, worktree `/viper/ptmp2/jinma/wt_sporder2`.
- **Work dir.** `/viper/ptmp2/jinma/sporder2_0925`:
  - `bin/`: ref = `git archive 7aabf6fd`; new = snapshots of this branch (`scripts/snap.sh`, `build.sh`).
  - `cpu/`, `gate/`, `gpu/`, `inp/`, `logs/`.
- **Binaries (md5).**
  - ref: none cpu d7b613bf, box cpu b31d126c, none gpu 2132bd47, box gpu 1fd86da3.
  - new:
    - `f3` (the commit): none cpu 3eee50af, box cpu 3df8f83e (gate_cpu_final.log: 39/39 BITWISE);
    - `f2`: none cpu 8a1552aa (comments differ);
    - `f1`: none cpu af8fbc84 (the same rad_m1 code; its sph_atm has no `lat_amp`);
    - box cpu c8068fe3, none gpu 26ff2ff3, box gpu 65e5da3e (f1 source).
  - `o4`: the vet_col tensor study. It has the same tensor code as the commit. It differs only in the outer
    Marshak face, which that study does not touch.
  - `pg`: ref plus the `sph_shell` `lat_amp` / `gas_vr` pgen options only.
- **Toolchain.** gcc 14 + openmpi 5, Release, MPI. GPU: rocm 6.3, gfx942, apudev, `HSA_XNACK=1`,
  `HSA_NO_SCRATCH_RECLAIM=1`.

## 1. What was first order on sp, and what was not

The brief suspected the face reconstruction. Measured on the ref binary:

**Already second order (unchanged):**
- **The transport operator** under `transport = implicit`.
  - The face-normal reduced flux F0 lives on the faces. Its equation is the centred two-point form over
    `dxface`.
  - E takes `dt A_f/V_i` per face.
  - M1SphDrr is the exact integrating factor.
  - M1SphCurv is the face mean of the two cell-centre values.
  - `implicit_enthalpy = plm` (the default for Eddington and vet_col) is the plm deferred correction on
    uniform and stretched r.
- **`implicit_recon` is not a lever here.** plm_dc acts only on the HLL part of `implicit_flux != central`.
  ImplicitInit makes `transport = implicit` fatal with anything but `central` on every mesh, so plm_dc has
  nothing to act on. The three refusals of SphericalS1Check stay: lifting them would validate nothing.

Self-convergence, L1_V(E_n - R E_2n)/L1(E):
- wedge r = 1..3, theta pi/2 +- 0.4, 4 x 4 angular cells;
- hesdirk2 with dt proportional to dx, so time is second order as well;
- `radial.py`, `lateral.py`.

| test | n (or 16 s) | L1 differences | orders |
| --- | --- | --- | --- |
| thick Gaussian shell, Eddington, rho kappa = 1000, uniform r | 32..512 | 2.18e-3 5.41e-4 1.35e-4 3.38e-5 | 2.01 2.00 2.00 |
| same, he4 stretched r | 32..256 | 3.96e-3 9.91e-4 2.47e-4 | 2.00 2.01 |
| thin shell, Eddington, rho kappa = 1e-3, uniform | 32..512 | 2.20e-2 5.65e-3 1.42e-3 3.53e-4 | 1.96 2.00 2.00 |
| same, stretched | 32..512 | 4.13e-2 1.12e-2 2.84e-3 7.12e-4 | 1.88 1.98 2.00 |
| intermediate, Eddington, rho kappa = 1 | 32..512 | 2.40e-2 5.98e-3 1.51e-3 3.76e-4 | 2.00 1.99 2.00 |
| thick shell, vet_col | 32..256 | 2.19e-3 5.47e-4 1.38e-4 | 2.00 1.99 |
| lateral theta mode, Eddington, thick (`pg`) | s = 1..4 | 2.45e-3 6.12e-4 | 2.00 |
| lateral theta mode, Eddington, thin (`pg`) | s = 1..8 | 1.67e-3 4.21e-4 1.06e-4 | 1.98 1.99 |
| lateral phi mode, Eddington, thick (`pg`) | s = 1..8 | 1.92e-2 4.84e-3 1.21e-3 | 1.99 2.00 |
| moving gas, v_r = 0.01 sin, thick, plm enthalpy, uniform (`pg`) | 32..128 | 1.13e-2 1.39e-3 | 3.0 |

- **The moving-gas row needs reflecting hydro walls.**
  - With `outflow` hydro walls the gas is first order (a zero-gradient ghost against v' != 0 at the wall).
  - That gave 1.04 / 1.02 on sp, and 1.07 on the Cartesian twin as well.
  - This was the test's hydro boundary, not the radiation.

**First order (fixed here):**
1. **The Marshak end face** (both closures). F = c q E of the end CELL is taken half a cell inside the face.
   - Steady Eddington atmosphere, analytic reference (`steady.py`):

     | r grid | L1 | order |
     | --- | --- | --- |
     | uniform, n = 32..512 | 6.68e-2 3.32e-2 1.65e-2 8.22e-3 4.10e-3 | 1.01 |
     | stretched, n = 32..512 | 3.41e-2 .. 1.90e-3 | 1.03-1.07 |

     The top cell error was 28 % at n = 32.
   - T-S2 free streaming (m1 closure, `ts2.py`): E/E_exact - 1 was 2.06e-2 1.04e-2 5.19e-3 2.60e-3
     (uniform) and 1.02e-2 .. 1.23e-3 (stretched).
2. **The vet_col formal solution**, i.e. the tensor handed to the faces (`tensor.py`).
   - Setup: f_K of the first build against a 2048-shell build.
   - It is order 1.1-1.4 in L1, and Linf 0.2-1.0.
   - Three causes, each found by switching it off:
     - the trapezoid in mu over p-ray nodes that crowd towards mu = 0 like sqrt(dr): O(dr^1.5);
     - S (and chi) taken linear in the PATH along every segment. Along a near-tangent segment r is
       quadratic in the path, which costs O(dr^1.5) in every non-core ray. This is the main cause;
     - the core rays started from the bottom CELL's E and F. With a reflecting inner boundary they
       started from an isotropic intensity, where the correct start is the mirror of the incoming one.

## 2. Changes (both new keys are read only when named; defaults bitwise)

- **`<rad_m1>/implicit_marshak_face = cell | linear`** (sp only; fatal on other meshes).
  - Files: `rad_m1_implicit.hpp` (`M1SphMarshakCoef`, `M1SphMarshakFaceE`), and the sp row and face
    branches of `m1_impl_asm` and `m1_impl_face`.
  - The face E is the linear extrapolation of G = r^2 E from the end cell and its neighbour. That is second
    order, and the free-streaming state r^2 E = const is exact.
  - It is IMPLICIT in the row: `ca` on the diagonal, `-cb` on the neighbour, ca > cb >= 0. The
    off-diagonal is non-positive and the row diagonally dominant, so the M-matrix sign pattern stays.
  - G_f is limited to [0, 2 G_end], so E_f >= 0. The limiter remainder is a deferred correction, zero on
    smooth profiles.
  - The implicit_bc_advect outflow uses the same face E.
  - Not used at the outer face under `vet_col_surface_q`. There q = H(face)/J(top cell) is built for the
    cell E. A face-J variant was tried: J(face) carries the grazing layer, and the atmosphere transient
    dropped to order 1.2-1.5.
  - A first version as a pure deferred correction stalled Picard: 118-180 NON-CONVERGED in the restart
    gate input. The implicit form gives 0.
- **`<rad_m1>/vet_col_order2 = true`** (`rad_m1_vetcol.cpp`, both kernels, team and column).
  - **Quadrature.** Base weights from the piecewise-quadratic interpolant in mu: the mean of the two
    overlapping triples, never across the core-edge node. Then the old 1, mu, mu^2 moment correction.
  - **S and chi linear in r along every segment.**
    - gb = mean of g(u) = (r(u) - r_a)/(r_b - r_a) per ray-segment is precomputed on the host
      (`vcol_gb`, `vcol_gb0`).
    - The weights are `VcolW2`: S = S_up + dS (u + a2 u(u-1)), a2 = +-(3 - 6 gb) clamped.
    - Both weights are >= 0.
  - **Core rays** start from E extrapolated to the inner face and the stored face flux f0x1. With
    `implicit_bc_x1min = reflect` they start from the mirrored incoming intensity (H = 0).
  - No change to the M1 solve and no change on the GPU layout.
- **pgen `rad_m1_tests2.cpp`** (test only, read only when named): `sph_shell` `lat_amp`/`lat_l2`/`lat_l3`,
  `gas_vr`; `sph_atm` `lat_amp`.
- `rad_m1_sph.cpp`: the header records the above.

## 3. Convergence with the new options (CPU)

| test | old | new |
| --- | --- | --- |
| steady Eddington atmosphere, uniform, n = 32..512 (tol 1e-11/1e-12) | L1 6.68e-2 .. 4.10e-3, order 1.0 | 9.96e-4 2.55e-4 6.45e-5 1.62e-5 4.06e-6, orders 1.96 1.98 1.99 2.00 |
| same, stretched | 3.41e-2 .. 1.90e-3, 1.03-1.07 | 1.85e-3 4.66e-4 1.17e-4 2.92e-5 7.32e-6, orders 1.99 2.00 2.00 2.00 |
| T-S2 free streaming, E/E_exact - 1 (uniform / stretched) | 2.1e-2 .. 2.6e-3 / 1.0e-2 .. 1.2e-3, order 1 | 1.6e-11 .. 2.6e-11 / 7e-12 .. 3e-11 (exact to the Picard floor) |
| vet_col tensor, extended atmosphere (T-S6 (ii) state), f_K L1, n = 32..512, ncore 8 | 1.87e-3 .. 4.54e-5, orders 1.28-1.40 | 5.58e-4 1.29e-4 3.03e-5 7.05e-6 1.49e-6, orders 2.11 2.09 2.10 2.25 |
| same, ncore 64, nsub 4 | 6.60e-4 .. 2.07e-5, orders 1.11-1.37 | 2.03e-4 .. 9.76e-7, orders 1.89-1.96 |
| vet_col tensor, Gaussian shell, reflect r, ncore 64, nsub 4 | L1 7.36e-4 .. 2.04e-5 (1.1-1.5), Linf order 0.3-0.9 | L1 6.00e-4 1.55e-4 3.94e-5 9.81e-6 2.34e-6 (1.95-2.06), Linf 1.93e-3 .. 1.16e-5 (1.95-2.00) |
| vet_col atmosphere TRANSIENT (T-S4 atmosphere, pure scattering, fixed dt, `atmt.py`) | 2.27e-3 6.26e-4 1.69e-4, orders 1.86 1.89 | 2.24e-3 5.82e-4 1.34e-4, orders 1.95 2.12 |
| lateral theta mode on that atmosphere, vet_col (`lateral.py th_atm_vc`) | 7.14e-3 1.89e-3 5.03e-4, orders 1.92 1.91 | 7.09e-3 1.81e-3 4.57e-4, orders 1.97 1.98 |
| same, Eddington | 1.10e-2 3.56e-3 1.29e-3, orders 1.63 1.47 | 7.15e-3 1.81e-3 4.55e-4, orders 1.98 1.99 |
| Gaussian pulse, rho kappa = 1, Marshak top, Eddington | 2.16e-2 5.42e-3 1.51e-3, orders 2.00 1.84 | 2.25e-2 5.62e-3 1.51e-3, orders 2.00 1.90 |

**T-S4** (the runs_5d/5e grey spherical atmosphere, vet_col + surface_q, be, steady; exact-transfer
reference `ts4_exact.npz` via `runs_5d_vetcol/ts4v.py`; ncore 8; `scripts/ts4.sh`):

| n | old L1 | old Linf | new L1 | new Linf |
| --- | --- | --- | --- | --- |
| 32 | 4.81e-4 | 2.44e-3 | 1.38e-4 | 6.81e-4 |
| 64 | 2.50e-4 | 1.85e-3 | 3.00e-5 | 2.28e-4 |
| 128 | 1.18e-4 | 1.05e-3 | 9.74e-6 | 1.71e-4 |
| 256 | 5.22e-5 | 4.89e-4 | 4.64e-6 | 1.74e-4 |

- The old column reproduces runs_5e to all digits.
- New is 3.5-12x better in L1. It flattens at ~5e-6 from n = 128 on (the reference or ncore 8), and Linf
  at 1.7e-4 (the top cell).
- L_out/L_in unchanged (7.8e-6 at 256).
- NON-CONVERGED 0; Picard mean 1.13.

**Not second order, open:**
- **A vet_col pulse whose K/J drops below 1/3** (a radiating shell).
  - The tensor clamp to [1/3, 1] puts a kink into f_K where K/J crosses 1/3. The dump at t = 1 shows the
    crossing at r = 2.68, next to the front.
  - The solve then converges at ~1 there: Gaussian pulse, rho kappa = 1, Marshak top.

    | arm | L1 differences | orders |
    | --- | --- | --- |
    | old | 2.66e-2 1.15e-2 6.92e-3 | 1.22 0.73 |
    | new | 2.57e-2 1.22e-2 7.24e-3 | 1.07 0.76 |

  - The error is localised at the outgoing front. With the tensor itself second order at t = 0, the clamp
    is the likely cause.
  - Relaxing the clamp is a closure change and was not done.
- **Reflecting OUTER x1 with vet_col.** The formal solution always takes a vacuum top, so the tensor is
  inconsistent with a reflecting M1 wall. The lateral mode on a reflecting-r shell does not converge, old
  and new alike (`th_mid_vc`).
- **Time.** The vet_col tensor stays lagged at U^n in both hesdirk2 stages (runs_5h: order 1.2-1.3 in
  time). The spatial studies use a fixed small dt, or dt proportional to dx for Eddington.
- One shell per build keeps the trapezoid weights: a quadratic weight turned negative, and the log says so.

## 4. Gates

- **CPU bitwise ref vs new** (`scripts/gate_cpu.sh`, `gate_cpu.log`): **39/39 BITWISE**.
  - A: Cartesian, `be` named, 10 cases.
  - B: Cartesian, key absent, 7 cases, including Milne vet_col.
  - C: sp S1/S2/S5 plus vet_col, `be` named, 15 cases.
  - F: sp key absent (the hesdirk2 default), 5 cases.
  - G: the new keys NAMED with their default values, 2 cases (payload after `<par_end>`).
  - Groups C, F and G rerun with the commit's none binary (f2): `gate_cpu_f2.log`.
- **Restart with the options on** (R): straight vs restart at cycle 250, 4 ranks.
  - `sph_sym_gas_rst` plus `implicit_marshak_face = linear` and `vet_col_order2 = true`, Eddington and
    vet_col.
  - The rst payloads are BITWISE (2/2 each).
  - NON-CONVERGED 0 in all four runs.
- **`tests_m1/gates/gates.py`** (new box binary): GATES: PASS, 15 pairs (`gates_eval.txt`).
- **`tst/test_suite/rad_m1`** (ATHENAK_M1_DATA = faces_0924/m1data): 3 passed (`pytest.log`).
- **cpplint** (`--filter=-build/include_subdir`) on the six changed files: clean.
- **flake8** (90 columns) on the scripts: clean.

## 5. GPU (apudev job 11971737, 2 x gfx942; `scripts/gate_gpu.sh`, `gpu/log.out.11971737`)

- **Bitwise ref vs new:** slab 1 rank be, box3d_nd 2 ranks be, box3d_nd key absent (hesdirk2 + vimp), and
  the He wedge vet_col default (96 x 128 x 128, 16 blocks, 20 steps): **4/4 BITWISE**.
- **Correctness vs CPU.** T-S4 vet_col wedge, n = 64, 100 steps, new options, 1 GPU vs 1 CPU rank:
  - max |dE|/max E 5.1e-11, F1 6.9e-11 (tolerance 1e-10);
  - identical Picard counts, NC 0.
- **Cost** on the He wedge grid (hesdirk2, vet_col, 40 steps, same binary, interleaved, 2 repeats):

  | arm | wall (s) | ms/cycle | build | inner its/solve |
  | --- | --- | --- | --- | --- |
  | default | 1.255 / 1.258 | 31.4 | 3.69 / 3.75 ms | 10.60 |
  | new (`implicit_marshak_face = linear`, `vet_col_order2 = true`) | 1.303 / 1.308 | 32.6 (+3.9 %) | 4.74 / 4.76 ms | 10.79 |

  NC 0 in all four.

## Files

- Python:
  - `solib.py` (input writer, readers, restriction);
  - `radial.py`, `lateral.py`, `steady.py`, `ts2.py`, `tensor.py`, `atmt.py` (run/eval drivers;
    env `SO_RAD`, `SO_PROB`, `SO_DT`, `SO_SUF`).
- `scripts/`: `snap.sh`, `build.sh`, `gate_cpu.sh`, `gate_gpu.sh`, `ts4.sh`, `atmt.sh`, `addkeys.py`,
  `gcmp.py`, `rstcmp.py`.
