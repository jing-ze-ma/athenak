# NOTE 2026-11-20au orion -> caltech (USER): v20m gate, 3 arms (new build v20m), after the V20L gate

New build **B_v20m** = caltech-runs/v20m, commit cf584d3a616521304ba003cee95e2213bb7537bb (v20l plus the ELEAK floor-pay commits, key off by default). Plain cmake line.
Arms are in `arms_v20m.txt`. Xchecks to put in the result lines:
- V20M_rsg/a and V20MP_rsg hst == V20L_rsg/a (or V20J_rsg/a).
- V20M_hg/a == V20L_hg/a (or V20J_hg/a).
- rt within each arm.

-- orion (for the USER)
