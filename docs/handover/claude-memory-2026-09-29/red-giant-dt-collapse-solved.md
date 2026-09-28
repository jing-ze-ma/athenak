---
name: red-giant-dt-collapse-solved
description: "The red giant dt collapse at t=1.9e5 s: full control set proving it was the EXPLICIT radiative diffusion at a cooling photospheric cell, and the fix that works - a DEEPER tau blend so the semi-implicit two-stream owns those layers (2026-09-08)"
metadata:
  type: project
---

**THE SYMPTOM (all of 2026-09-08).** Every cubed-sphere red giant run lost its timestep
at **t = 1.89-1.91e5 s**, dt falling from ~30 s to 0.05 s over a few hundred cycles.

**THE CONTROL SET.** The collapse is invariant to everything except the RT handover:

| varied | still collapses at 1.9e5? |
| --- | --- |
| IC: adiabatic column, MLT alpha=1.7, MLT alpha=3 | yes |
| inner BC: open, wall at 0.5 R, wall at 0.03 R (RADIATIVE zone) | yes |
| inner wall leaking vs sealed (`wall_noflux`) | yes |
| top sponge on / off | yes |
| nx1 = 192, 320, 512 | yes |
| **blend tau_lo/hi 1/10 -> 10/100 or 30/300** | **NO - passes 2.05e5, dt flat 30.66 s** |

**THE MECHANISM.** `Conduction::NewTimeStep`'s report (added today) named the cell:
photosphere, r = 3.429e12, tau ~ 10.8, w = 1 (fully diffusive). It had COOLED from
5124 K to 2945 K; kappa_R fell 12x with the H- opacity, kappa_rad rose 10x, and
dt_cond = dx^2 rho c_v/(w kappa_rad) followed: dt1 = 1.79 s while dt2/dt3 (angular) were
3800/4450 s. Purely radial diffusion. NOT a boundary, NOT the IC, NOT resolution.

**THE FIX (the user's idea).** `Conduction::NewTimeStep` skips any cell whose blend
weight is zero, so moving the handover deeper hands those stiff layers to the two-stream,
which is semi-implicit (048dff30) and unconditionally stable. dt at cycle 0:
1/10 -> 18.8 s; 10/100 -> 30.6 s; 30/300 -> 30.6 s. The last two agreeing to six digits
means conduction no longer binds at all -- 30.6 s is the HYDRO timestep.

**THE TRAP that nearly killed the idea.** `nabla` built the initial column radiative
wherever tau < rad_tau_hi, so raising rad_tau_hi ALSO built a different (far hotter,
wildly superadiabatic) star, and the first test looked like deepening the blend made dt
WORSE (18.8 -> 10.5 -> 4.6 s). `problem/ic_tau_rad` (25692d16) decouples them; set it to
10 and vary only the blend.

**Still open at handover:** whether the emergent flux relaxes toward L (the gate was
0.877 -> 0.556 with the old setup, `tools/grid/flux_gate.py`), and whether convection
actually develops. See [[red-giant-session-2026-09-08]].
