---
name: dhj-dt-limited-by-alfven-floor
description: "MEASURED from a dump, not argued: dhj MHD dt is set by the RADIAL Alfven CFL in the floored outer atmosphere (r/Rp ~ 1.28), in a cell sitting exactly on dfloor. Raising dfloor x10 would buy 2.7x."
metadata:
  node_type: memory
  type: project
---

Question 2026-08-18: "in the mhd run, is the timestep limited by the Alfven speed at the top?"
**Yes.** Located exactly, by reproducing `MHD::NewTimeStep` on `L_id13/bin/dhj.mhd_w_bcc.00002.bin`
(ideal, bbot = 10 G, t = 2e6). Reconstruction gives dt = 0.750 against the run's ~0.7 -- good
enough to trust the location.

## The limiting cell

- **Direction is RADIAL** (dt_r = 0.750, dt_phi = 0.925, dt_theta = 12.9)
- **r/Rp = 1.2796** -- 85% of the way up the domain, NOT the very outermost shell
- rho = 5.477e-12, and **dfloor = 5.44e-12, so rho/dfloor = 1.007: the cell is ON the floor**
- |B| = 42 G (amplified 4x from bbot = 10), vA = 1.81e7 cm/s, cs = 6.7e5 -> **vA/cs = 27**

Radial profile of the local limit: 11-30 s below r/Rp = 1.23 (phi-limited, rotation+sound),
then a cliff -- 1.7 at i=48, **0.75 at i=52**, ~0.85 above. Only **1.27% of the volume** sits
within 1% of dfloor and that 1.27% sets dt for the entire run.

With B deleted from the same snapshot dt = 3.42, so **the field costs 4.6x** here. (Not 14x --
the remaining gap to the bbot=0 runs' dt = 11 is that those never develop this atmosphere.)

## Raising dfloor buys dt as sqrt(1/rho), re-evaluated on the snapshot

x2 -> 0.97 (1.29x) | x5 -> 1.48 (1.97x) | **x10 -> 2.01 (2.68x)** | x20 -> 2.43 (3.25x)

This is the exact mirror of [[dhj-ideal-xe-floor-relaxation]], which measured dfloor **x1e-3**
costing **5.9x** by the same mechanism (max vA there went 2.5e6 -> 1.8e7). Note the base run's
max vA has since GROWN to 1.8e7 by t = 2e6 -- that note's 2.5e6 was at t = 1.22e5, one rotation.

**Do not treat "raise dfloor" as a free speedup.** It is a physics change: the floor currently
lives only in the top 12% of the domain in radius, and x10 pushes it deeper. The static
re-evaluation above holds rho's structure fixed; the real run would evolve differently.

## How to redo this

`/tmp/.../scratchpad/dtlimit.py`. Two traps:
1. **The .bin dump's `x2v/x2f/x3v/x3f` are the LOGICAL grid, not physical angles** -- x2 comes
   back as -1..-1 (the input's sentinel) so dx2 evaluates to 0. Reconstruct theta yourself:
   `theta = pi/2*(1+sinh(a*(2*xi-1))/sinh(a))` with `a = f_stretch_theta = 2`, plus the
   Mignone-corrected centres in `coordinates.cpp:695-699`. Radial dx1 = r_r-r_l IS correct
   from x1f.
2. dx1 = dr, **dx2 = r*dtheta, dx3 = r*|sin theta|*dphi** (`coordinates.cpp:697-699`), and
   AthenaK is Heaviside-Lorentz so **vA = B/sqrt(rho)**, no 4pi.

Related: [[dhj-ideal-vs-general-cost]], [[dhj-ideal-xe-floor-relaxation]], [[dhj-highB-outer-bc]].
