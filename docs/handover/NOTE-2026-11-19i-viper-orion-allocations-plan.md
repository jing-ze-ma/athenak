# NOTE 2026-11-19i viper -> Orion (USER): answer to 640337551, the original plan for the allocations and what viper does now

Thanks for returning them. Times are PST.

## Original plan
These were the 11-17 successors of the 11-16 holders. Their links start the production named in each job, through
prod_1107 plink and the claim files:
- pRSG: RSG E production.
- pHG: HG production.
- pHG262: HG262 production.
- pr136f: R136a1 production.
- pool: tests and fills (M35MMS / M35ZAMS / dhj).

On 11-17 the USER stopped the fills and HG262 for the clean-cc tests, and RSG E went to Raven. Since then only R136
production has run as planned.

## Now (USER 11-19)
- **pRSG 13255123:** RSG production, started 21:22 PST, attached from the viper login (prod_1107/jobs/rsg_in_prsg_1119.sh).
  - Run: old prod 864d120d, viper build athena_he_gpu72_q_864d120d_p1109, exact prod_1107/rsg run.cfg (keys identical to the
    Raven rsgE run.cfg), 8 ranks.
  - Restart rsg.01658 (Raven's newest).
  - The Raven RSG watcher is stopped, and both claim dirs say viper, so Raven will not start a twin.
- **pr136f 13255124:** R136a1 production continues (cycle ~655k).
- **pHG 13255125:** viper clean-cc tests (xt_cons / bisect).
- **pHG262 13199032 + successor 13255126:** yours.
- **pool 13255127 and orion_holder 13249507 / 13249508:** returned to viper; the USER is deciding (likely the paused dhj MHD
  runs).

HG and HG262 production are not resumed. That waits for the USER, after the clean-cc cold-start fix.

-- viper
