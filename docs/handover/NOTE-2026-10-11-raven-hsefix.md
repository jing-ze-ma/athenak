# NOTE hsefix (hsefix-1011): A100 frA/frA0 check, and the sph1d HSE IC on AG Car B

Branch `hsefix-1011` (from prodfinal-1011 8279c0a6), worktree /raven/u/jinma/ATHENAK/wt-hsefix, MPCDF Raven.
Work dir /raven/ptmp/jinma/hsefix_1011. Commit **eef1f8f6** (code). Nothing pushed.

## 1. A100: frA vs frA0 not bitwise (j.31079211)

**Verdict: the cause is a deck difference. GPU nondeterminism is ruled out.**

The resolved PAR_DUMPs of frA and frA0 (rst 00000) differ in one key only: `rad_m1/implicit_det_reduce`. frA has 0, frA0
has 1. (he_ic_rad_hse prints as `true` vs `1`; that is the same value.)

How the 0 got into the deck:
- prod/agcar_{A,B}_fresh.athinput:259 name `implicit_det_reduce = 0`. mkdeck.py copied it from a **CPU** PAR_DUMP.
- The key's default depends on the build (src/rad_m1/rad_m1_implicit.cpp:1677-1683): on a device build it defaults to
  det_ok (true for the production solver); on a host build it defaults to false.
- The CPU pair was bitwise because the CPU default is 0 as well.
- No other GetOrAdd default depends on the build.

Where the difference shows:
- It is already present at cycle 0, after the ic_rad fs_relax implicit solves: m1_f2/f3 2.9e-10 on about 9000 cells.
- By cycle 20 it reaches velz 7.8e-9.

Repeat job 31080690 (gpudev, 4 A100, athena_he_a100_18ff91c8):

| pair | result |
|---|---|
| frA_1 vs frA_2 (det_reduce 0) | BITWISE; also bitwise vs old j.31079211/frA |
| frA0_1 vs frA0_2 (det_reduce 1) | BITWISE; also bitwise vs old frA0 |
| frAf_1 vs frAf_2 (he_ic_rad file, 2 cycles) | BITWISE |
| frA0_1 vs frA_1 | maxrel 7.755e-9 (velz), the original number |

Fix:
- In prod/agcar_{A,B}_fresh.athinput, set `implicit_det_reduce = 1` or drop the line.
- prod/restart_keys_B_runB.txt also has `rad_m1/implicit_det_reduce=0`.
- mkdeck.py and mkrestart.py should skip this key.

## 2. AG Car B: "sph1d: the VET iteration failed" with he_ic_rad_hse = true

### 2.1 Root cause

The failure was reproduced on CPU at init. It is a host solve, run with 1 rank and nlim = 0. The iteration was
instrumented with a debug build (the instrumentation is not committed); runs are in cpu/runs/j.31080789.

**The cause is not the Fe zone.** There, Gamma = kt F/(c g) is 0.83-0.85, below 1.

**The AG Car B file column is super-Eddington just below R_ph:**
- The layer is at T 20-25 kK, tau_R 0.5-3, r 6.99-7.05e12 cm.
- kt rises there from 0.47 to 0.62 and Gamma reaches 1.3.
- The column already has a density inversion there: 678 fine cells with rho increasing outward.

What the fixed-T HSE march does in that layer:
- With Gamma > 1 it grows rho outward, and kappa_t(rho) grows with it.
- Above the layer (Gamma_es about 0.5 at the RE temperature, 15.6 kK) the radiation-supported atmosphere is 1e3-1e5 x
  denser than the file's gravity-only isothermal atmosphere.
- The outer E-solve <-> HSE Picard loop therefore runs away:
  - max |dln rho| per outer is 6.0, 2.7, 2.1, 1.9, 1.2, 0.75, 1.1, ...;
  - the Gamma > 1 zone wanders (Gamma_max 1.2-1.9; 46 at outer 9);
  - at outer 9 the secant fails (residual 0.9998) and rho overflows;
  - at outer 10 tau is NaN, the solve grid collapses to 16 NaN points, and the HsSph1dPass check fires
    (he_star_m1.cpp:965 in 18ff91c8).

AG Car A has Gamma ≤ 0.80 everywhere and a monotone outer residual (0.54 → 5.7e-7 in 29 outers), so it never hit this.

### 2.2 Fix (eef1f8f6, src/pgen/he_star_m1.cpp, HSE part of the sph1d solve only)

1. **Gamma clamp** `problem/he_ic_rad_hse_gmax` (default 0.9, as he_bc_hse_gmax):
   - In the HSE march, g_eff = (1 - gmax) g wherever kt F/c > gmax g.
   - Unclamped cells use the original expression, so they are bitwise unchanged.
   - Clamped cells keep a net outward force (Gamma - gmax) g. The file IC has the same imbalance there.
2. **Damping after the first rise of the outer residual** `problem/he_ic_rad_hse_damp` (default 0.5):
   - From then on, ln rho ← (1 - w) ln rho + w ln rho_HSE.
   - Before the first rise w = 1, which is the old update.
3. **Separate outer limit** `problem/he_ic_rad_hse_maxit` (default 400):
   - It was he_ic_rad_maxit = 60, which is shared with the VET passes.
   - A WARNING is now printed when the loop does not converge. Before, an unconverged loop passed silently.
4. The log line reports the clamped-cell count, the outermost clamped r, and the outer at which damping started.

Both the clamp and the damping are needed (init-only scans, cpu/runs/j.31081014 and j.31081130):

| variant | outers | result |
|---|---|---|
| clamp only, gmax 0.9 / 0.95 | 60 | oscillates, dln rho about 2-3 |
| damp 0.5 or 0.25, no clamp (gmax 2) | 200 | oscillates, dln rho 0.3-5 |
| clamp 0.8 | 23 | converges, but clamps the whole Fe envelope (9222 cells) and R_ph moves to 0.93 |
| **clamp 0.9 + damp 0.3 / 0.5 / 0.7** | 154 / **146** / 267 | **converges** to 1e-6; the same column in all three (R_ph 7.12966e12) |
| clamp 0.95 + damp 0.3 / 0.5 | 400 | does not converge |
| (halving damping) | 200 | stalls at w = 0.05 |

The B column that results (gmax 0.9, damp 0.5):
- R_ph (tau_R 2/3) moves from 7.0472e12 to 7.1297e12 (+1.2 %).
- The HSE atmosphere reaches the 1e-16 floor at 7.44e12 (1.056 R_ph0). The file's gravity-only atmosphere reached it at
  7.21e12.
- 538 fine cells are clamped, all at r ≤ 7.0355e12 (the super-Eddington layer).
- Cost: 146 outers / 705 VET passes. That is 49-57 s on a CPU host (A: 25 s; on the A100 host A takes 12.6 s).

### 2.3 Validation (reduced CPU AG Car, 480x32x32, 4 ranks, FRESH prod decks, production dt, 60 cycles; cpu/runs/agc.31081366)

The deck copies are inp/agcar_{A,B}_v*.athinput: outputs every 10 cycles, plus the three new keys named.

**Bitwise checks against prodfinal 18ff91c8** (all 18 files each):
- A fresh (sph1d re_atm hse) eef1f8f6 vs 18ff91c8: **BITWISE**. The IC log line is identical (331 VET passes, 29
  outer). Clamp 0 cells, damping never triggers.
- B fresh with hse false (the current prod B deck) eef1f8f6 vs 18ff91c8: **BITWISE**.

**AG Car B photosphere** (window 0.97-1.03 R_ph0; rho drift = max over shells of |<rho>/<rho>_c0 - 1|):

| arm | Picard mean / max | max\|v\| c10 / c30 / c60 [km/s] | rho drift c10 / c30 / c60 | L_top c1 / c60 |
|---|---|---|---|---|
| **B hse true (fixed)** | **2.73 / 8** | **0.22 / 0.53 / 0.82** | **0.003 / 0.013 / 0.026** | 1.0001 / 1.0005 |
| B hse false (current prod deck) | 2.92 / 8 | 4.0 / 12.2 / 19.0 | 0.091 / 0.62 / 3.54 | 1.0000 / 1.0006 |
| B he_ic_rad = file | 2.84 / 7 | 4.3 / 10.9 / 18.9 | 0.086 / 0.59 / 3.19 | 1.0000 / 1.0006 |
| A hse (unchanged; for scale) | 3.02 / 6 | 0.07 / 0.18 / 1.06 | 0.001 / 0.017 / 0.14 (the 1.03 R_ph floor shock, NOTE-icgen) | |

All arms: rc 0, NON-CONVERGED 0.

What the B numbers show:
- With hse, the B photosphere is quasi-static: v is 23x lower and the rho drift 140x lower at cycle 60 than with
  hse false.
- Without hse it expands at about 19 km/s and rho rises 4.5x near 1.005 R_ph.
- The remaining hse motion is at 0.992 R_ph0: outward <v_r> 0.64 km/s, growing about linearly. That is the clamped
  Gamma > 0.9 layer.
- max |v| < 0.95 R_ph (convective seed) and > 1.05 R_ph (floor gas) is the same in all B arms.

### 2.4 A100 (gpudev job 31081367, athena_he_a100_eef1f8f6, full size)

Full size, 4 A100, deck copies with `implicit_det_reduce = 1`:
- **B hse true** (inp/agcar_B_vh_g, 30 cycles):
  - rc 0, NON-CONVERGED 0, Picard 2.53 / 8.
  - The sph1d solve is identical to CPU (146 outers, R_ph 7.12966e12, 538 clamped cells). It takes 46 s on the GPU
    node host.
  - Photosphere window: max|v| 0.22 / 0.52 km/s and rho drift 0.003 / 0.013 at c10 / c30, the same as on CPU.
  - L_top c1 / c30: 1.0001 / 1.0002.
- **A**, eef1f8f6 vs 18ff91c8 (10 cycles): **BITWISE** (8 files).

## 3. Open

- **gmax 0.9 is a choice.**
  - Above about 0.93 the B outer loop does not converge.
  - Below 0.85 it starts clamping the Fe-zone envelope of B (Gamma 0.83-0.85).
  - The clamped layer is genuinely super-Eddington in this column. No static solution at fixed T was found there.
- B is much slower to start: about 50 s CPU host vs 0.5 s with hse false.
- The prod B FRESH deck still names hse false. To switch, set `he_ic_rad_hse = true` and name the three new keys
  `he_ic_rad_hse_gmax = 0.9`, `he_ic_rad_hse_damp = 0.5` and `he_ic_rad_hse_maxit = 400`. The A deck should name them
  too, because command-line keys must exist in the deck.
