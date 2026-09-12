---
name: dhj-restart-loses-gravity-potential
description: "FIXED (30d21859 + 4cfdc329): deep_hot_jupiter_rt restarts used to lose the gravitational potential and the cell-centered field"
metadata: 
  node_type: memory
  type: project
  originSessionId: 85cd606f-82a3-4a02-a9d2-5e410a10595a
  modified: 2026-08-09T10:51:54.082Z
---

**The SAME bug hit `solar_convection` and was fixed separately on 2026-09-04 — see
[[solar-convection-restart-gravity-bug]]. Check any pgen that fills solver-owned arrays.**

Found 2026-08-09. **Affects real science runs — every `run/dhj_*` setup has
`etotgrav = true` and they are long enough to restart.** Nothing to do with the general
EOS: ideal and general are affected identically.

`ProblemGenerator::UserProblem` in `src/pgen/deep_hot_jupiter_rt.cpp` does
`if (restart) return;` (line ~153), immediately after enrolling `user_srcs_func` and
`user_bcs_func`. But the gravitational potential arrays `phicc0` and `phi0.x1f` are filled
**inside** `UserProblem` (lines ~462, ~481, ~691) and nowhere else — `hydro.cpp` /
`mhd.cpp` only `Kokkos::realloc` them, which zero-initialises. So after a restart the
potential is identically ZERO, and `use_etotgrav` (which converts total <-> internal energy
using rho*phi) and the well-balanced scheme both silently operate on phi = 0.

Proven by experiment (hydro, spherical polar, 40 cycles, restart at cycle 20):

| | straight run | restarted run |
|---|---|---|
| `etotgrav = true`  | t = 953.650 | t = **783.851** (diverged) |
| `etotgrav = false` | t = 949.375 | t = 949.375, fields bitwise identical |

So with etotgrav off, restart is exact; with it on, the run silently changes physics at the
restart point. The timestep drops ~36%.

**FIXED in `30d21859`.** The early `return` is gone; only the two genuinely one-off blocks
(the `probini` initial condition and the initial magnetic field) are now guarded with
`if (!restart)` individually. The background kernels `probwb` / `wbgrav` / `wbcc` do not
read `u0`/`w0`, so skipping the initial condition does not affect them. Verified: straight
and restarted runs now agree exactly for hydro + etotgrav under both `eos = ideal` and
`eos = general`, and for the dynamic well-balanced scheme; fresh-start results unchanged.

**Second bug, also FIXED, in `4cfdc329`.** MHD restarts were separately non-exact
(~4e-4 in `bcc2`, ~1e-2 in `velx`). Cause: `u0` and `b0` ARE restored exactly, ghost zones
included (restart files store `nout1 = nx1 + 2*ng`, see `src/outputs/restart.cpp`), but
**`bcc0` is derived and NOT stored**, so it starts at zero. The inner-x1 user boundary
swaps ghost magnetic energy *statefully* —
`u0(IEN) -= 0.5 bcc0^2` ... recompute `bcc0` from `b0` ... `u0(IEN) += 0.5 bcc0^2` —
so with `bcc0 = 0` the first subtraction removes nothing and the magnetic energy is
DOUBLE COUNTED on the first step. A fresh start fills `bcc0` over the whole block in
`pgen_b0`, which is why it only bit on restart. Fix: fill `bcc0` from the restored face
fields on restart, matching the ConsToPrim interpolation.

Now every field is bitwise identical across restart except velocity residuals of one
float32 ULP.

**Diagnostic technique that worked, reuse it:** dump every cycle (`<output>` `dt` smaller
than the timestep), match dumps by TIME not file number, and find the first cycle that
differs. Here the restart-point state was already bitwise identical and the divergence
appeared in the single step after, which immediately pointed at state used during a step
rather than at the restart I/O. Then collapse the difference over axes
(`d.max(axis=...)`) to localise it — it sat at `i = is` with `ix1_bc = user`.

**Measurement floor:** `.bin` output is SINGLE precision, so ~6e-8 relative is the noise
floor and apparent differences of one ULP are not real.

See [[general-eos-project]].
