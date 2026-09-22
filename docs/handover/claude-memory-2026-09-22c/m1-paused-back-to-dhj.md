---
name: m1-paused-back-to-dhj
description: User decision 2026-09-21 night - after the three VET tensor tests, stop all M1 work; next project is the deep hot Jupiter (dhj)
metadata:
  type: project
---

User, 2026-09-21 night: once the three tensor tests (frozen / tilted / tau-refreshed Eddington tensor in
the implicit solve) are done and recorded, DO NOT TOUCH the M1 radiation module any further; the next
step is to go back to the deep hot Jupiter (dhj) work.

**Why:** the M1 module exists to prepare for a VET scheme whose Eddington tensor will come from short
characteristics (numerical closure, refreshed every hydro step); analytic chi(f) fixes do not serve
that goal. Tonight established: every chi(f) closure (Levermore, Minerbo, Kershaw) is unstable in the
optically thin top of the 2-D implicit solve because (chi, n) are refreshed from the cell's own F
(local algebraic feedback at c dt/dx ~ 7e3); a fixed tensor (Eddington) is clean over 5 turnovers.
No further M1 closure / Newton / optimisation work was asked for.

**How to apply:** finish + commit the tensor tests, update the handover doc and [[rad-m1-design]],
then switch to dhj: start from index-hot-jupiter.md (sub-index) and the latest dhj notes; ask the
user which dhj question comes first if it is not obvious. Do not start new M1 agents or items
(multigrid, halos beyond the running agent, (E,F) Newton) unless the user reopens it.
