---
name: solar-convection-restart-gravity-bug
description: "FIXED and pushed (ef9561e2): solar_convection restarts ran with ZERO gravity and dt collapsed to 1e-135"
metadata:
  type: project
---

Found and fixed 2026-09-04. **Same bug class as [[dhj-restart-loses-gravity-potential]],
which was fixed in the dhj pgen months earlier and never applied here.** Assume any pgen
that fills solver-owned arrays itself has this bug until checked.

**Symptom:** restart a `solar_convection` run and the first cycle reports
`dt = 4.449524e-135`. No error, no NaN message.

**Cause:** `phi0.x1f/x2f/x3f` and `phicc0` are ALLOCATED by Hydro/MHD
(`hydro.cpp:174-177`) but FILLED only by the pgen, in kernels that sat after
`if (restart) return;`. Hydro leaves them zero, so the restart ran with no gravity at all
while the state was still a stratified atmosphere -- the pressure gradient had nothing to
balance it. Misleading detail: the base state and the sponge parameters HAD already been
hoisted above that return, so the pgen looked restart-safe.

**Fix** (`ef9561e2`, branch `fix/solar-convection-restart-gravity`, pushed to the fork,
parented on `origin/polar-average-perf` -- NOT on the local branch, which was 108 commits
behind): the return now happens after the potential kernels, and only the
initial-condition kernels (`probini`, `probwb`) are skipped. The face potential is written
inside `probwb`, so a dedicated `gravfaces` kernel reproduces those writes -- the potential
is just phi = -g*x1 (`probwb`'s theta/lam/phi bookkeeping is dead code, the rotation term
is commented out). A fresh run writes identical values twice and is unchanged.

**Still broken on purpose:** restarting with `wellbalance_static`/`wellbalance_dynamic` is
now a FATAL ERROR rather than silently wrong -- the well-balanced FACE background
(`w0facewb`) is built only inside `probwb` and cannot be reconstructed on the restart path.
Fix that if a well-balanced run ever needs to restart.

**Verify a restart in one line:** dt on the first restarted cycle must match the dt the
previous run ended on (here 0.231 vs 0.238). Do this check on EVERY restart of this pgen.
