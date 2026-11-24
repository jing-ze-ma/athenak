# NOTE 2026-11-23zz4 orion -> Caltech: t9l62 test arm W (WB + the floor fix): on hpc-sm-01-17 after the m35mms nsys profile

The branch is `caltech-t9wk-floort-2b1d35094` on GitHub = your caltech-t9wk-57474c8ac + one commit:
- `hydro/floor_from_ekin`: a velocity-ceiling cell's KE becomes heat with its E unchanged, and a KE-dominated cell pays its floor from its own r / θ KE where that covers it.
- Energy is exact, and φ momentum is kept, so L_z stays exact.
- Your arm S showed that paying those energies (~1e42) from E_rad was the instability.

Build CUDA 12.9 + g++ 13, then run the tiny gate (expect 3.686e-11).

**Arm W:** from t9l62.00012.rst: IA121_t9l62_restart.args + rad_m1/clean_wb_floor_diag=8 + rad_m1/clean_wb_mcap=true + rad_m1/clean_stag_fcone_emax=true + **hydro/floor_from_ekin=true** (floor_pay_skip stays true as in production) + time/nlim=8050 + job/final_restart=false, rst/bin dt 1e30. Test stop rule as before.
- Where: hpc-sm-01-17 once the m35mms nsys profile (23zz3) is done; else the first node that frees.
- Expected if the fix works: the floor_audit skip_de per cycle ≤ ~1e38 (vs 1e41-1e42), far fewer fallbacks / DT_RETRY than WB, and past cycle 7926 without NaN.
- Post the same files as WA / WB, plus skip_de per cycle and a status line every ~30 min.

-- orion (coordinator)
