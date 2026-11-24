# Notebooks: the clean VET tests with plots

One Jupyter notebook per test of `tests/vet_clean/run_tests.py`. Each notebook runs its test
through the suite (or, if the run is already there, only loads it), plots the result against
the exact solution, the reference or the conservation criterion, and prints the suite's
PASS / FAIL line.

| notebook | suite tests | what it plots |
|---|---|---|
| `shell.ipynb` | `conservation_sp`, `conservation_cs` | drift of total energy, mass and L_z vs cycle |
| `thermal_eq.ipynb` | `thermal_eq` | E(t), T(t) vs the exact ODE, their errors, total-energy drift |
| `order.ipynb` | `order_time`, `order_space` | L1 error vs dt and vs N, log-log with a slope-2 line |
| `wb.ipynb` | `wb_sp`, `wb_cs` | max \|v\|/c_s vs time; radial profiles of `rhse_sp` at t = 0 and t = 10 |
| `nlte.ipynb` | `nlte_eps1e-1`, `nlte_eps1e-4` | J vs tau against the exact solution (`ref/nlte_exact_eps*.txt`) |
| `momentum_box.ipynb` | `momentum_box` | per-cycle table and plot of the linear-momentum drift (a report, no criterion) |
| `restart.ipynb` | `restart_sp`, `restart_cs` | straight vs restarted run: the difference (0) |
| `radwave.ipynb` | `radwave` | wave amplitude vs the analytic decay; phase-speed and damping errors vs the trapping parameter x |
| `wedge_sym.ipynb` | `wedge_sym` | F_phi per phi column across the wedge; the wedge-mean F_phi of each wedge vs the criterion |
| `floor_pay.ipynb` | `floor_pay` | total-energy drift from cycle 1 on; the floor energy paid from the radiation per cycle |
| `idort.ipynb` | `idort_*` (8 tests) | crossing beams (image and cuts vs the exact beams), shadow (image, cuts vs the exact attenuated beams), emitting cylinder (E, F vs exact transport), radiative shocks (rho, T vs Lowrie & Edwards), cooling at large dt (T(t) vs the exact ODE) |
| `radiation.ipynb` | `radshock`, `radshock_m5`, `diffusion`, `marshak`, `advect_pulse`, `advect_uniform`, `shadow` | the classic radiation tests: shock profiles vs Lowrie & Edwards; diffusing pulses vs the exact Gaussian and the diffusion-coefficient error vs tau_cell; Marshak-wave profiles; the advected vs the shifted resting pulse; E maps with and without the clump and the geometric umbra |

## How to run

```bash
pip install numpy matplotlib jupyter                       # nothing else is needed
export ATHENA_EXE=$PWD/build_tests/src/athena              # the CPU build of QUICKSTART.md 3a
cd tests/vet_clean/notebooks
jupyter notebook                                           # and open a notebook, or:
jupyter nbconvert --to notebook --execute --inplace thermal_eq.ipynb
```

Settings (environment variables, read by the first cell of each notebook):
- `ATHENA_EXE`: the binary (default `../../../build_vettest/src/athena`, the one
  `tests/vet_clean/run_tests.sh` builds).
- `M1_NB_WORKDIR`: where the runs go (default `runs/` next to the notebooks). A test whose
  record `WORKDIR/<test>.json` exists is loaded, not rerun: delete `WORKDIR/<test>*` to rerun.
- `M1_NB_THREADS`: OpenMP threads per run (default 4).

With 4 threads per run most notebooks take under a minute; `order.ipynb`, `wb.ipynb`
(the cubed-sphere `wb_cs`), `idort.ipynb` and `radiation.ipynb` take 2-5 min.
`nbutil.py` holds the two shared helpers: `suite()` (runs `run_tests.py --only TEST --keep
--record`) and `verdict()` (the suite's PASS / FAIL check against `ref/references.json`).
The committed notebooks were executed with the CPU binary of commit 30d0e3461 (Python 3.13,
numpy, matplotlib, nbconvert 7); their key plots are exported to `docs/vet_clean/figures/`
(`export_figures.py`). All twelve also run with `nbconvert` on Python 3.6 (numpy 1.19, matplotlib 3.3, nbconvert 6.0). `jupyter lab` needs JupyterLab,
which `pip install jupyter` does not install on Python 3.6; use `jupyter notebook` there.
