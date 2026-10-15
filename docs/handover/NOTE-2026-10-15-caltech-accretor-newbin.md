# NOTE 2026-10-15 (Caltech): accretor RHD on the new binary (fin51-r-fix merged), S5 stream smoke passes the first impact

Branch accretor-rhd-1015 @ af45fa46 (fork). Design note docs/dev/accretor_rhd_design.md section 12. Tests only,
0.97 H200-h. Runs: /resnick/groups/carnegie_poc/jingze/accretor_1015.

## Merge
- accretor-rhd-1008 (7fedd194) + origin/fin51-r-fix (b16b66db) -> eb880b72. Conflicts in rad_m1.hpp and
  rad_m1_implicit.cpp (ImplicitSolve face flux) resolved as unions (accretor per-column outer inflow face + the
  littests-1011 boundary beams).
- dt_retry visitors regenerated with tools/dtretry_gen.py (clang 16.0.6): only the accretor members were added.
- ry_per_accretor registers its per-stage budget state (stage_ctr_, acc0_/acc1_) with dt_retry (f960ed04);
  m1.ibadv_x1max joins the minimal list in driver/dt_retry.cpp (fdb18517; empty in every other problem).
- Binary (nofma, HOPPER90): builds_accr/bin/athena_gpu_ry_per_accretor_fdb18517_nofma, md5 4e89a4377108b50365614addfcb03c11.

## Gates (1 H200)
- dt_retry: restart bitwise; forced retry (fac 1) bitwise; with views = all and poison of 19 non-minimal Views
  bitwise; forced fac 0.5 with or without m1.iw + m1.fref_wacc checkpointed bitwise (both are scratch; their
  poison failures are 0*NaN reads, not lagged state).
- Keys absent vs the old accretor code (a844c97d, same nofma tree): NOT bitwise on rad decks, because of the
  fresh-run defaults of 10-09..10-12 (implicit_gas_src, implicit_rad_inertia, implicit_nofb_work,
  implicit_f0_init, vet_bc_diff, vet_seg, vet_source_noesrc, implicit_face_opac_n, implicit_pos_floor_solve).
  With those pinned: S3a-final bitwise; S5 stream decks dt identical, state 1e-7 relative after 20 cycles
  (unexplained round-off-level code difference). A FRESH S5 start with the new defaults has a 2.6x smaller dt
  from cycle 7 (stream injected from t = 0).

## S5 option (a) stream smoke (restart of rhd/s5/smoke_prof at t 0.02366)
- With time/dt_retry = true + the full ROL2U set: passes the first impact, runs to t 0.0300, NON-CONVERGED 0,
  no FATAL. 21 retries (eint<=0) in t 0.02432-0.02443. Picard mean 4.15 (was 12-15). 0.51-0.56 s/cycle (was
  0.95); vet_col tensor builds are 62 % of the time. The 10-09 failure column stays at 17.8-19.7 kK (was 10 kK).
  Shock 119 kK at r 9.43 by t 0.030, P_rad/P_gas 5.4. Energy budget residual 1e-5 of int L dt (tol_rel).
- dt_retry without ROL2U: fails at the same cycle (3642); dt below dt_min 1e-7 at cycle 3654.

## Open for the user
- Longer S5 test or production setup with these keys; fresh starts would take the new defaults (pin or accept).
- The vet_col tensor build cost (333 ms per build, every step) is now the main cost.
