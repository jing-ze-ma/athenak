# dt_retry / force_reference = clean patches for clean-cc-push16 (5e33f056a)

Made by the Orion fix2 worker, 2026-11-16. The patches apply with `git am` on 5e33f056a
(branch clean-cc-push16-fref in the worktree /orion/ptmp/jinma/clean_1115/audit/fix2/wt_v).

| patch | content | +/- |
|---|---|---|
| 0001 | dt_retry visits `cln_fvi` and `cln_fref_set`; `m1.cln_fvi` goes into kDtrMinimal | +7 |
| 0002 | cell force_reference = clean: the reference is captured at every iterate of the first coupled solve and marked set only when that solve is accepted (same rule as clean-cc-stag2-wb 6dc54b2e3) | +4 -1 |

Merge note: 0002's line in CleanSolve does not depend on the layout. stag2-wb 6dc54b2e3 has
the same line with `!cell &&` at the same place. Keep the layout-free one. The staggered
capture there already follows that rule.

## Evidence (CPU, OpenMP, gcc 13 Release, -ffp-contract=off; binaries in bin/)

vbase = 5e33f056a, v1 = + 0001, v2 = + 0001 + 0002. Run dirs: ../t3/runs, ../t3/wbr, ../reg.

| test | vbase | v1 | v2 |
|---|---|---|---|
| q_cel (sp, cell, fv1, extrap1), 8 cycles: forced retry at cycle 4, fac 1, vs the straight run | **not bitwise** (from output 6, 1.3e-15) | bitwise | bitwise |
| no retry, vs vbase: q_cel, hcyl x0.5, thin_p, cs closed, rhse_sp (reg) | - | bitwise | bitwise |
| no retry, rhse_sp + force_reference = clean, clean_init = given | - | bitwise | changes (m1 1.6e-5 after cycle 0) |
| wb.sh (force_reference = clean, t = 4) max v/c_s fv1 / fv / lld / cart | 1.18e-9 / 2.76e-9 / 4.72e-9 / 7.38e-9 | - | 1.21e-9 / 2.96e-9 / 4.75e-9 / 7.38e-9 |
| rhse_sp fref, retry at cycle 0 (fac 0.5) vs fresh run at cfl/2: gas momentum after cycle 0, max diff (max \|mom\| 5e-7) | **9.3e-3** | 5.9e-5 | 2e-19 (by construction: the first coupled solve does not move the gas) |

What the failed attempt left in place (found on clean-cc-stag2-fix2, whose visitor has every cln_* member):
- `cln_fvi`: the start intensities of the fv sweep. That sweep iterates to a tolerance, so a different start gives a different result. On the q_cel cell deck, restoring cln_fvi alone made the forced retry bitwise. Restoring everything except cln_fvi did not.
- `cln_fref_set`: the reference force is captured once, at the first coupled solve. A retried first cycle kept the reference from the failed attempt, which was taken at the larger dt. The momentum row of the table shows the effect.
- Extrap history (`cln_mn/mp/mc`, `cln_have_*`, `cln_tn/tm`): not affected. The handle swap and the time tests keep a retry consistent without restoring anything. On stag2, poisoning these Views after the restore was still bitwise.

## Open on this head (not fixed here)

1. **A retry of the FIRST cycle is not a fresh step at the new dt.** This holds in vbase, v1 and v2: rhse_sp with clean_init = given gives E 24 % off after cycle 0. The failed attempt initialises clean arrays lazily: CleanInit allocates cln_w, cln_x and the history buffers, and cln_w holds the last Newton dm. None of that is rolled back, so the retried solve takes another Newton path ("it 1" against "it 0").
   - On stag2 the full visitor rolls `cln_init` back with every other cln_* value, and the same test is bitwise.
   - On this head, two attempts failed:
     - `cln_init` plus the extrap handles alone: FPE in the retry.
     - The full generator-style visit of this class (exp_full_visitor_abort.patch): heap corruption ("free(): invalid next size"; under gdb a View refcount in ~RadiationM1). The likely cause is a host container of the fv joint sweep or ring (cln_fvls, cln_fvkl, cln_fvdm, cln_ring ...) restored against arrays the attempt rebuilt.
   - Retries at later cycles are bitwise (row 1).
2. **Side note:** CleanClosureAfterSolve also runs after a non-converged backward-Euler solve that set dtr_fail (t2_fail stays false outside a stage). It writes cln_mc, cln_tm and cln_have_c from the failed state. This is harmless only because dt_retry puts the cycle back. Without dt_retry, the accepted non-converged solve feeds the extrapolation.
3. **Cubed sphere, reflecting radial walls:** CleanGhostFill copies the Cartesian D into the ghost. The proper ghost is the mirror (I - 2nn) D (I - 2nn) about the face normal.
   - A kernel for it is in ../cs_mirror_exp.patch (written against stag2-fix2's CleanGhostFill).
   - On cs closed and feed1 it changes the result only at round-off (max 2e-15). The E and L_z budgets are unchanged (L_z is round-off with either ghost), so it was not committed.
   - CLN_KP and CLN_G are read only in the cell itself, never in a ghost (CleanCellOutStress, CleanCellSphSource), so their wall and pole ghosts need no parity.
