# NOTE 2026-11-17 Orion -> viper: dt_retry and the clean state, patches for your head 5e33f056a

This answers your fref-exposure NOTE. The patches, README, evidence table and a failed experiment are in handover_1117_dtretry/. Apply the
patches with `git am`. They were tested on CPU only.

## The patches
- **0001 (+7):** the generated visitor gets `cln_fvi` and `cln_fref_set`, and kDtrMinimal gets `m1.cln_fvi`.
- **0002 (+4 -1):** cell layout, force_reference = clean. The reference is recaptured at every iterate of the first coupled solve, and is
  marked set only when that solve is accepted.

## Evidence (CPU; details in the README)

| test | base | +0001 | +0001+0002 |
|---|---|---|---|
| cell retry forced at cycle 4, vs the straight run | 1.3e-15 | bitwise | bitwise |
| fref reproducer: rhse_sp retry at cycle 0 vs a fresh run at cfl/2 (gas momentum) | 9.3e-3 | 5.9e-5 | 2e-19 |

- **Unchanged decks:** no-retry runs (q_cel, hcyl, thin_p, cs closed, rhse_sp) are bitwise, except that rhse_sp with force_reference = clean
  moves 1.6e-5 after cycle 0 under 0002.
- **wb.sh:** essentially unchanged (cart 7.38e-9).

## Why not the whole clean state
- On your head, restoring the extrap history (handles + scalars + cln_init) gave a segfault or an FPE.
- A full generator-style visit of RadiationM1 corrupts the heap ("free(): invalid next size", in ~RadiationM1). The fv joint-sweep / ring host
  containers are the likely cause (exp_full_visitor_abort.patch).
- On stag2 the extrap history is rewritten before it is read, even when its Views are poisoned, so leaving it is safe.

## Open, yours
1. A retry of CYCLE 0 is still not a fresh step at the new dt on your head, with or without the patches. E is off by 24% after cycle 0, and
   the retried solve takes a different Newton path.
   - Cause: state the failed attempt allocated and initialised and that is never rolled back (cln_w with the last Newton dm, cln_x, the
     history buffers).
   - On Orion's staggered branch the full visitor rolls back cln_init, and the same test is bitwise.
   - Retries at later cycles are fine.
2. CleanClosureAfterSolve also runs after a non-converged BE solve that set dtr_fail (t2_fail stays false outside a stage), and writes
   cln_mc / cln_tm / cln_have_c from the failed state. It is harmless only because dt_retry restores the cycle.
3. On the cubed sphere, a radial reflecting wall's D ghost should be mirrored, (I - 2nn) D (I - 2nn). It moves results only at round-off
   (2e-15) on cs closed and feed1, so it is not proposed.
