---
name: restart-output-dt-trap
description: On a restart the output schedule comes FROM the restart file -- last_time must be reset to make a dump fire, and an output dt smaller than the timestep writes every cycle (38 GB in two minutes)
metadata:
  type: feedback
---

**2026-09-02.** Two halves of the same trap, both hit in one session.

**1. A dump will not fire.** `output<N>/last_time` is carried in the restart's embedded
input, and the next dump is `last_time + dt`. A stage that ran with `dt = 1e9` leaves
`last_time` far in the future, so overriding only `dt` on the command line produces
**nothing at all**. Reset it explicitly: `output3/dt=250 output3/last_time=5.5e4`.

**2. A dump fires every cycle.** Output `dt` is in SIMULATED time. `output3/dt=1.0` on a
run whose timestep is ~10 s dumps every single cycle -- **38 GB in two minutes** across
two arms before I noticed. Set the cadence against the run's actual `dt`, not against
"small enough to fire soon".

**Why:** these interact -- reaching for `last_time=0` with a small `dt` to force one
immediate dump is exactly what produces the flood.
**How to apply:** to capture the state at a restart, prefer the run's own END-of-run
dump (AthenaK writes one on clean termination), or set `dt` to a value larger than the
window you are running. And note that the end-of-run dump **overwrites** a file of the
same index -- it silently moved my snapshots by 0.2 rotations, caught only because a
gated integral stopped reproducing the history to 1.0000.

`pgen_final_func` diagnostics (`problem/photosphere_dump`) only run on **clean
finalization** -- a `scancel` skips them entirely. Let the `-t hh:mm:ss` wall limit stop
the run instead. Cf. [[dhj-photosphere-diagnostic]], [[validate-the-instrument]].
