---
name: dhj-conservation-check
description: measured mass/energy drift of the dhj general-EOS + correlated-k run over t=4.32e5 - mass +0.21%, internal energy +0.56%, and how to measure it (hst is unusable)
metadata:
  type: project
---

Measured 2026-08-24 on the `ckrepro` reproducer ([[session-state-2026-08-24]]).

## Getting the numbers at all

Two obstacles, both worth remembering:
1. **`<output1> hst dt = 8.64e6`** in `inputs/mhd/deep_hot_jupiter_rt_eos.athinput` is 10 %
   of `tlim`, so a run that reaches t ~ 7e5 writes **only the t = 0 header line**. Lower
   it (1e4 or so) for any run being watched for a blow-up.
2. **The hst volume weight is Cartesian** -- see [[hst-cartesian-volume-on-spherical]].
   Even with samples, the hst columns are not physical integrals on this mesh.

Recipe that works: restart a snapshot for a couple of cycles into a scratch dir with the
output forced on, then integrate the `.bin` yourself.

```
srun build/src/athena -r <run>/rst/dhj.0000N.rst -d <scratch> \
     output2/dt=1.0 output2/last_time=0.0 time/tlim=<restart_time + ~60>
```
`Driver::Initialize` skips the initial output when `res_flag` is set, so the run **must
take at least one cycle** -- and `time/nlim` cannot be used to bound it, because the
restart's `ncycle` is already far past any small `nlim` and the loop then never runs. Bound
it with `tlim` instead. Scripts kept in `/orion/ptmp/jinma/Athenak/hstprobe/`.

Volume element for the integral: `dV = (r2^3-r1^3)/3 * (cos th1 - cos th2) * dphi`, block
bounds from `filedata['mb_geometry']` (the `.bin` header's x2min/x2max are the literal
`-1.0` from the input, since the pgen sets the real range). `.bin` is single precision --
fine here, it costs ~1e-7.

## Result, t = 13.4 -> 4.32e5 s (general EOS, max_eta 1e14)

| | t = 13.4 | t = 4.32e5 | rel |
| --- | --- | --- | --- |
| mass | 3.18409e26 | 3.19077e26 | **+2.10e-3** |
| internal E | 4.81386e38 | 4.84097e38 | **+5.63e-3** |
| kinetic E | 1.77e29 | 3.77e33 | +2.1e4 |
| magnetic E | 3.89e30 | 7.93e32 | +203 |

KE and ME grow by orders of magnitude but are ~1e-5 of the internal energy, so the total
energy budget is the +0.56 % internal-energy gain -- the atmosphere heating toward
radiative equilibrium, not a conservation failure. Neither drift is blow-up-like.

**Radial structure of the mass change** (this is the interesting part): the inner shells
r < 9.8e9 LOSE mass (-0.6 % at the very bottom), r = 9.9-11.0e9 GAINS (up to +12 % per
shell), and the tenuous outer region r > 11.2e9 is EVACUATED (-85 % per shell, though it
holds little mass). The net +6.7e23 g is a ~30 % residual of those opposing terms, so it
is redistribution plus a small boundary supply -- the ix1 user BC is the natural suspect
and has not been checked. Profile saved at `/orion/ptmp/jinma/Athenak/hstprobe/massprof.npz`.

Densities in the outer region are ~5e-10 at t = 0 and ~9e-11 at t = 4.32e5, both well above
`dfloor = 7.26e-12`, so the floor is NOT the source of the mass drift.
