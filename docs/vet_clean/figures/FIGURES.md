# Figures of the tests

The key plots of the tests, exported from the executed notebooks by `export_figures.py` (run it
from the repository root after re-executing a notebook, and commit the PNGs with the notebook).
The documentation site embeds them: on the Tests page in this order, on the Properties page
after each property whose "Test:" line names one of the figure's tests, and the `star3d_*`
figures on the 3-D stars page.

Columns: `figure` is the PNG in this directory; `source` is the notebook and `#N` its N-th PNG
output (1-based, in cell order); `tests` are suite test names (comma-separated); `caption` says
what the plot shows and the pass criterion, in one line.

| figure | source | tests | caption |
|---|---|---|---|
| conservation.png | tests/vet_clean/notebooks/shell.ipynb#1 | conservation_sp, conservation_cs | Relative drift of total energy, mass and `Lz_tot` from cycle 0 on the closed rotating shell (left spherical polar, right cubed sphere; grey: gas + radiation L_z without the ledger). Pass: energy, mass and `Lz_tot` below 1e-13 (dashed). |
| restart.png | tests/vet_clean/notebooks/restart.ipynb#1 | restart_sp, restart_cs | Radiation energy per cycle of the straight run (dots) and of the run restarted half-way (crosses). Pass: every dump and history row identical (difference 0). |
| order.png | tests/vet_clean/notebooks/order.ipynb#1 | order_time, order_space | Relative L1 difference of successive levels against dt (left, deck D) and against N (right, linear wave), with a slope-2 line. Pass: order of the finest pair >= 1.7. |
| radwave.png | tests/vet_clean/notebooks/radwave.ipynb#1 | radwave | Left: decay of radiation-acoustic waves from the free-streaming (x = 0.01) to the trapped limit (x = 100); points code, lines the analytic root. Right: phase-speed and damping errors against x. Pass: phase speed within 3e-3 (dashed), damping within 5e-3 per radian (dotted). |
| wb.png | tests/vet_clean/notebooks/wb.ipynb#1 | wb_sp, wb_cs | Largest \|v\|/c_s of an envelope in hydrostatic and radiative equilibrium against time (spherical polar to t = 10, cubed sphere to t = 2.5). Pass: below 1e-8 (dashed). |
| wb_2.png | tests/vet_clean/notebooks/wb.ipynb#2 | wb_sp | The spherical-polar envelope at t = 0 and t = 10: density, internal and radiation energy (left; the curves overlap) and the radial velocity, of order 1e-10 (right). |
| wedge_sym.png | tests/vet_clean/notebooks/wedge_sym.ipynb#1 | wedge_sym | Left: F_phi per phi column on the 20, 3.7 and 7.3 deg wedges. Right: the wedge-mean \|ΣF_phi\| / Σ\|F_r\| of each wedge. Pass: below 1e-8 (dashed); the 3.7 deg wedge without `vet_gd_wall_sym` stops with the wedge FATAL. |
| nlte.png | tests/vet_clean/notebooks/nlte.ipynb#1 | nlte_eps1e-1, nlte_eps1e-4 | Mean intensity J against optical depth in a two-level-atom atmosphere (epsilon = 0.1 and 1e-4) with the exact transport solution (left) and the relative difference (right). Pass: L1 within 5 % of the reference value. |
| thermal_eq.png | tests/vet_clean/notebooks/thermal_eq.ipynb#1 | thermal_eq | Gas temperature and radiation energy relaxing to equilibrium against the exact ODE at two time steps (left), their relative errors (middle) and the total-energy drift (right). Pass: T error below 2e-4 (dashed), order >= 1.8, energy drift below 1e-13 (dashed). |
| floor_pay.png | tests/vet_clean/notebooks/floor_pay.ipynb#1 | floor_pay | A temperature floor that acts in most cells every cycle, paid from the radiation: total-energy drift from cycle 1 on (left) and the floor energy paid per cycle (right). Pass: drift below 1e-13, energy paid > 0, unpaid 0, no `CLEAN WARNING`. |
| momentum_box.png | tests/vet_clean/notebooks/momentum_box.ipynb#1 | momentum_box | Linear-momentum drift max_k \|P_k - P_k,0\| / Σ\|p_gas\| per cycle in a periodic 3-D box (gas + radiation face momenta). A report, no criterion; here at round-off. |
| star3d_sp_profiles.png | docs/vet_clean/notebooks/star3d.ipynb#5 | star3d_sp | Spherical-polar BSG deck on 4 A100 GPUs: shell-averaged density and radiation energy density at t = 0 and at the end (9.8 d); the atmosphere's initial radiation has streamed out. (Its luminosity plot is dominated by the t = 0 burst, L_top = 1600 L_in; expected: L_top / L_in between 0.96 and 1.07 over the second half.) |
| star3d_cs_lum.png | docs/vet_clean/notebooks/star3d.ipynb#2 | star3d_cs | Cubed-sphere RSG deck on 4 A100 GPUs: L_top, L_bot and L_in against time (left) and their ratios to L_in (right). Expected: L_top / L_in between 0.986 and 0.997 over the second half. |
| star3d_sp_slices.png | docs/vet_clean/notebooks/star3d.ipynb#3 | star3d_sp | Spherical-polar deck at its last dump: density contrast rho/<rho>_shell - 1 and log10 E in a meridional and the equatorial plane. |
| star3d_cs_slices.png | docs/vet_clean/notebooks/star3d.ipynb#4 | star3d_cs | Cubed-sphere deck at its last dump: the same in the planes z = 0 and y = 0. |

Tests without a figure: `dt_retry` and `layout_tile` (bitwise comparisons of two runs; no
notebook).
