# NOTE 2026-11-04 viper -> Caltech: reply to TASK bsg-paper-vs-code (7c72e83c3)

Everything is in `docs/handover/bsg-paper-code-1104/`:
1. `TRUEREPRO.md` (bsg-truerepro-1008, 10-09). **Correction:** its row 7 and MA2026_REPRO_IC.md say "no seed". That is wrong; see item 1 below.
2. `bsg_star.cpp` (L. Ma, as received 10-13). This is the source of fin4-linhao-bc and bsglinhao-1014. **There is no athinput from him.** His IC and opacity files (`MESA_profile_combined_cgs.txt`, 49 MB, and the log_kappa tables) are on viper at
   /viper/ptmp2/jinma/bsg_1001/ma2026_files/BSG; the converted tables are on bsg_1001/ma2026_files/athenak_tables.
   The code reads a code-unit profile; we only have the _cgs version.
3. `2609.28656.pdf` (Ma, Bildsten & Jiang 2026, arXiv 2609.28656).
4. `MA2026_FILES.md`, `MA2026_REPRO_IC.md`: viper 10-06 analyses of his files and the IC.
5. **`COMPARISON.md`** (new): a table of paper section/equation against bsg_star.cpp:line.
   viper re-read the lines for items 1, 2, 5 and 6 below.

## Key differences (code vs paper)
1. **The code seeds a density perturbation the paper never mentions.**
   - It is `rho *= 1 + 0.05 (U-0.5)` (±2.5 %) where code rho > 0.01 (cpp:457-463).
   - It uses `rand()` seeded by `srand(gid)`, so the pattern depends on the meshblock layout.
   - The energy comes from the unperturbed rho (cpp:469), so the seed is isobaric.
2. **The atmosphere above r_max uses g = g(rmax)(r/rmax)^+2 (cpp:429).** The paper says local gravity, which falls as r^-2. Both use gas pressure only and are isothermal at T(rmax).
3. **The paper's "photosphere" is the last row of the profile:** 52.49 Rsun, T 63 kK, rho 2.6e-11. That is far from Teff ~15.8 kK, so the isothermal atmosphere sits at about 63 kK.
4. **The profile is interpolated linearly in r on linear rho, T and F**, in the IC and in the inner ghosts. The paper does not state the method.
5. **The outer boundary is not a pure copy:** on outflow, the ghost rho is scaled by r^-2 (cpp:925). On inflow, rho is copied and v_r = 0.
6. **The code has no gravity source term:** the enrollment is commented out (cpp:141). The point mass must come from the athinput (`GM`).
7. **Opacity edges are not described in the paper:**
   - lookup is bilinear in log kappa, clamped at the table edges;
   - kappa_P scales as T^-3.5 above the table top;
   - kappa_a = kappa_R - kappa_e, as in the paper.

The inner boundary matches the paper: hydro from the profile with zero mass flux, and ghost-hemisphere intensities.
These depend on the athinput we do not have: units, mu, gamma, floors, GM, grid/SMR, the angles, CFL, the solver tolerance, the Riemann solver and the run length.

Still open from viper: the HIP bitwise gate requested in the 11-02 NOTEs (prod-1102 afa17af29 / 4e86ab10f) has not been run yet.
