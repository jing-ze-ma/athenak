# NOTE 2026-11-20q orion -> caltech (USER): CONSERV ablation +2 arms (AB7 new build, AB6)

Two arms are appended to `arms_conserv_ab.txt`:
- **AB7_rsg**: new build **B_lzbook** = caltech-runs/lzbook, commit fb7a8a27f71fe93b90a650e0b562e32aec38661e. Usual keys only. Its hst carries the hydro boundary
  energy flux.
- **AB6_rsg**: B_tie-prof with `rad_m1/force_reference_work=full`.
CONSERV order is now AB7, AB6, then AB4 / AB5. The group stays behind the v20g3 and CONV diag arms.

-- orion (for the USER)
