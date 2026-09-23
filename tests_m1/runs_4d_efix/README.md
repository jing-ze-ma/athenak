# runs_4d_efix: `implicit_bc = efix` broke down at the T7 outflow end (fixed)

Branch `m1-efix` (from rt-integration 22418492). CPU only (login node, serial Release
build, nice 10). Runs: `/viper/ptmp2/jinma/efix_0923` (`diag/`, `gate/`, `bitwise/`,
`expl/`). Runner: `run.sh` here. Table: `validate_0923/scripts/t7_table.py`.

## Symptom (runs_3z_validate, constraint 4)

- Lowrie-Edwards T7 with `<rad_m1> transport = implicit`, `implicit_bc_x1min|max = efix`.
- After 1-1.5e-10 s the x1max end cell drops to `e_floor` and F/cE goes to 0.999.
- Both `be` and `hesdirk2`. `marshak` ends do not show it.

## Diagnosis (M2, be, N512, cfl 0.4, t = 2e-10; runs in `diag/`)

- The breakdown is an **x2 mode** (nx2 = 4) at the **x1max end column**.
  - F2 in the end cell grew by x400 per 5e-11 s: 4.9e10, 2.0e13, 8.5e15, 3.6e18.
  - It doubles per cell toward x1max.
  - The pattern in j is (a, b, -a, -b), and gas vely grows with it.
  - The growth is faster at higher N: x40, x400 and x4000 per 5e-11 s at N256, N512
    and N1024.
- Control arms:

| arm | outcome |
| --- | --- |
| `implicit_vimp = false` | still grows |
| `implicit_offdiag = none` | bitwise identical to base (Eddington closure) |
| efix at x1min, marshak at x1max | stable |
| F2 of the efix cell set to 0 (no transverse force) | still breaks down |
| **E of the efix cell kept at the Dirichlet value** | **stable, end F2 ~ 1e3** |
| no gas coupling in the efix cell | stable |

## Cause

`src/rad_m1/rad_m1_implicit.cpp`, the write-back kernel `m1_impl_wb`, `ep -= (ch/cl)*work`
(line 6950 at 22418492).

- The efix row is replaced by E' = EN (`m1_impl_asm`, about line 6266). That row has no
  transverse coupling.
- After the solve, the write-back still subtracted the gas work term from E in that cell.
- So the "Dirichlet" cell was moved every step, by a j-dependent amount.
- Nothing in its row damps x2 structure, so the kicks accumulated.
- At the outflow end they fed back through the gas (radiation force, work) into a growing
  x2 mode until E hit the floor.

## Fix

In `m1_impl_wb`, the work term is not applied to an `M1_IBC_EFIX` end cell (`efc`).

- The cell keeps E' = EN: its initial value under be, and zero stage slope under
  hesdirk2.
- The gas in that cell still gets the momentum and heat exchange.
- Every other cell and every other BC runs the same arithmetic as before.

## Gates (t = 2e-10, efix at both ends; `gate/`)

Gate: L1 < 2 % in rho, T_gas and T_rad, and NON-CONVERGED = 0.

| run | L1 rho | L1 T_gas | L1 T_rad | NC | marshak (runs_3z) rho / Tg / Tr |
| --- | --- | --- | --- | --- | --- |
| m2 be N512 c0.4 | 4.70e-4 | 5.21e-4 | 1.08e-4 | 0 | 4.80e-4 / 5.21e-4 / 1.21e-4 |
| m2 be N1024 c0.4 | 4.48e-4 | 3.58e-4 | 7.28e-5 | 0 | 4.52e-4 / 3.88e-4 / 1.07e-4 |
| m2 h2 N512 c0.4 | 4.37e-4 | 4.80e-4 | 5.93e-5 | 0 | 4.48e-4 / 4.87e-4 / 7.69e-5 |
| m2 h2 N1024 c0.4 | 4.40e-4 | 3.48e-4 | 4.99e-5 | 0 | 4.46e-4 / 3.83e-4 / 9.20e-5 |
| m5 be N512 c0.4 | 1.19e-3 | 7.24e-4 | 3.52e-4 | 0 | 1.98e-3 / 5.11e-4 / 1.69e-4 |
| m5 be N1024 c0.4 | 5.32e-4 | 8.14e-4 | 4.18e-4 | 0 | 1.81e-3 / 4.48e-4 / 9.65e-5 |
| m5 h2 N512 c0.4 | 8.60e-4 | 5.12e-4 | 1.94e-4 | 0 | 1.72e-3 / 3.93e-4 / 9.11e-5 |
| m5 h2 N1024 c0.4 | 4.48e-4 | 7.74e-4 | 3.39e-4 | 0 | 1.80e-3 / 4.38e-4 / 6.56e-5 |
| m2 be N1024 c0.2 | 4.36e-4 | 3.60e-4 | 5.87e-5 | 0 | |
| m2 h2 N1024 c0.2 | 4.37e-4 | 3.53e-4 | 5.01e-5 | 0 | |
| m5 be N1024 c0.2 | 4.63e-4 | 7.97e-4 | 3.81e-4 | 0 | |
| m5 h2 N1024 c0.2 | 4.48e-4 | 7.74e-4 | 3.29e-4 | 0 | |

- The x1max end cell now holds its initial E exactly.
  - M2 N1024: 3.18285103825375781e12 at every output.
  - F/cE there is 4e-4 (M2) and 1.1e-3 (M5).
- hesdirk2 took 1 backward-Euler step in each run, as the marshak runs do, and had 0
  stage fallbacks.

## Bitwise

- M2 N512 with marshak ends (the input default) gives bitwise identical m1 and hydro_w
  tabs at t = 2e-10 between the base binary and the fixed one, under both be and
  hesdirk2 (`bitwise/`, `diag/H_marshak`).
- The fix touches only efix cells.

## Explicit transport with the fixed-E ghost (RadM1ShockBC), M2 N512 nx2 = 4, t = 1e-10

- L1 rho, T_gas, T_rad at t = 1e-10: 4.25e-4, 4.94e-4, 8.88e-5. This is the same as the
  1-D explicit run.
- End-cell F/cE was 5e-4 and E moved by -3.6e-5 relative.
- F2 is exactly 0 everywhere: the explicit scheme keeps x2 symmetry to the bit, so there
  was no seed.
- Structurally the explicit path has no analogue of the bug.
  - Its Dirichlet state lives in the GHOST cells, re-imposed each stage by the user BC.
  - The end active cell is an ordinary cell, with its own transverse fluxes.
