# TASK 2026-11-06 Caltech -> viper (USER): `his1d` on his real radial grid, unseeded

## Why
`his1d_30d` (384 uniform, dr 0.15625 at the wall, seeded) is our Athena++ 1-D reference, but two things differ from what a shell mean of his 3-D sees:
- **Radial grid:** near the wall his 3-D has dr 0.3125 Rsun (root 192 on 20-80 Rsun; level 1 only at 35-55). `his1d` has half that at the wall. The wall layers we study are only a few cells thick.
- **Seed:** his cell-by-cell seed averages out over a shell in 3-D. In a column it stays as radial noise (rms N^2 ~22e-10 at t = 0). That the shell-mean noise is small in 3-D is our reasoning, not a measurement.

Our AthenaK `Lg` (`vstarfix_1106/col/Lg`) has this grid: 48 x 0.3125 + 128 x 0.15625 + 80 x 0.3125 Rsun, unseeded, GM 19.78 Msun, dt ~4.4 s. So far its KE is 1.36e40 at 0.5 d against your 1.17e40.

## Arms (his code, units, tables, GM 260.795, cfl_rad 0.025, nzeta 5; everything else as `his1d_30d`)
1. **`his1d_g`:** his radial grid, seed off (`problem/seed_amp`). This is the main arm.
2. **`his1d_gs`:** his radial grid, seeded as coded, only if cheap. It separates the grid effect from the seed effect.

His radial grid means real SMR in 1-D if that works (root 192, level 1 at 35-55 after block rounding). If 1-D SMR is awkward, use an equivalent non-uniform grid and say how it differs at the two jumps in dr.

Run to at least 25 d, so it covers his 17-19 d inflow, as `his1d_30d` did.

## Please report
Same columns as `his1d_30d`, so that we can compare directly:
- `*_wall_profiles.txt` (first 24 cells, double precision, cgs), at 0.5, 1, 2, 3, 5, 7, 10, 14, 18, 22 d if possible. Dumps every 0.5 d would be better still.
- KE in the pi/3 sr wedge against time at hst cadence (a text file on this branch), so that we can compare the oscillations with `Lg`.
- N^2 on faces, T/IC, rho/IC, Lc0/L_*, dE, dEgrav, dU, R(1e-10 g/cc), Jacobi stats.
- What changes against `his1d_30d`: any wall layer at dr 0.3125, the size and period of the ~1 d oscillation, and any inflow.

Short NOTE on this branch, please; first results when they're easy, not only at the end. Don't delay the Orion 3-D for this.
