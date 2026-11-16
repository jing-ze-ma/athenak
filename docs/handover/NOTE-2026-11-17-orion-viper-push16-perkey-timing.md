# NOTE 2026-11-17 Orion -> viper: per-key CPU timing of the push-16 defaults (answer to TASK-2026-11-17-viper-orion-push16-perkey-timing)

## Short answer
- **The push-16 CPU slowdown does not come from the new default keys.** On 0678599eb, with every new key off
  (clean_gas_exact, clean_fv_joint, clean_pc_factor, clean_yres = false, clean_rot_ring = off), the code is as slow as with the defaults.
  It is still +4..+26 % slower than push 15 on 6 of the 8 decks. That arm is bitwise to push 15 on the reg harness (see A arm below), so it
  does the same work: the cost is overhead in code paths that the keys do not gate.
- No key that is a default on 0678599eb costs > 5 % on any deck. clean_gas_exact, clean_pc_factor and clean_rot_ring = live can stay default.
- **clean_fv_joint = true costs +7 % on cs_thin1** (spread 0.3 %), so opt-in is right by our rule. It is also **5..15 % faster** on
  cs_feed1, RT feed and feed_x1. Both heads show the same pattern.
- Where to look for the remaining cost: in cs_closed, the "clean stage steps" wall goes from 1.77 s (push 15) to 2.20 s (all keys off),
  with the same output. The ungated push-16 commits are 57a08c09f (gas eint less d phi, |B|^2/2), 0a3df762d (relative Newton dE norm),
  a97e368cb / c249955c4 (x1-line team-owned blocks, factor store, FV level reduction forms), plus whatever ff74a049e, d589df2d7, b7964ecaa
  and 44c0fedce added on the paths that run with the keys off. We did not bisect commits. Tell us if you want that.

## Setup
- Your message moved the head to 0678599eb while the 5e33f056a run was half done. We cancelled that run (job 230811) and reran everything on
  0678599eb (job 230835), with the arms from your message plus clean_fv_joint = true. The partial 5e33f056a data are at the end.
- Both jobs ran on one dedicated p.exclusive node, orion097. Each deck ran standalone, one after another, with the same decks, threads and
  method as job 230771 (OMP_PROC_BIND=false, hwloc binding none, SLURM_/PMI_/PMIX_/OMPI_ unset in the srun step).
- Each arm ran twice. Round 1 went in the order push 15, def, gas, pcf, rng, off, fvt, and round 2 went in reverse, so drift does not favour
  either end. Job 230835 took 69 min wall in total: about 4.9 min per binary-arm set, plus 2 min for the A arm.
- The keys go into a copy of each -i deck (after `<rad_m1>`) and on the command line, as athenaA does (wrappers S/athenaH_*).
- Figure: cpu time used / cycles (s/cycle), the mean of the 2 runs. ±x % is the half-spread of the 2 runs. **Bold** is % vs push 15.
  [brackets] is % vs the defaults of the same head, which is the cost of the key.
- Cycles per run: flz 40, cs_closed 4, cs_feed1 10, cs_thin1 100, thin_f 150, RT feed 16, feed_x1 16.
- Binaries (recipe gpufix_reg/build.sh: gcc 13, openmpi 4.1, OpenMP, -ffp-contract=off):
  - push 15 dfe21a746, md5 0f714f25;
  - 0678599eb, md5 8916d177;
  - 5e33f056a, md5 9716e1da.
- Noise: the half-spread is ≤ 2 % on most entries and up to 5 % on the short decks (flz_Q, cs_closed, cs_feed1, feed_x1). Push 15 itself
  moved by up to 10 % between this job and job 230771 (feed_x1 1.06 -> 1.16, flz_Q 0.28 -> 0.31). So we treat differences below about 5 %
  as noise.

## 0678599eb (job 230835, s/cycle)

| deck | push 15 | 0678599eb defaults | gas_exact=false | pc_factor=false | rot_ring=off | all off | fv_joint=true |
|---|---|---|---|---|---|---|---|
| flz_q | 1.747 ±1% | 1.969 ±1% **+13** | 1.885 ±2% **+8** [-4] | 1.886 ±1% **+8** [-4] | 1.852 ±4% **+6** [-6] | 1.872 ±3% **+7** [-5] | 1.897 ±1% **+9** [-4] |
| flz_Q | 0.302 ±4% | 0.295 ±0% **-2** | 0.298 ±3% **-1** [+1] | 0.293 ±3% **-3** [-1] | 0.300 ±4% **-1** [+1] | 0.297 ±2% **-2** [+1] | 0.293 ±2% **-3** [-1] |
| cs_closed | 1.293 ±1% | 1.582 ±1% **+22** | 1.622 ±5% **+25** [+2] | 1.663 ±4% **+29** [+5] | 1.562 ±1% **+21** [-1] | 1.624 ±2% **+26** [+3] | 1.682 ±3% **+30** [+6] |
| cs_feed1 | 1.057 ±4% | 1.213 ±1% **+15** | 1.248 ±3% **+18** [+3] | 1.170 ±3% **+11** [-4] | 1.170 ±0% **+11** [-4] | 1.175 ±0% **+11** [-3] | 1.029 ±1% **-3** [-15] |
| cs_thin1 | 0.556 ±0% | 0.585 ±0% **+5** | 0.583 ±1% **+5** [-0] | 0.590 ±0% **+6** [+1] | 0.585 ±0% **+5** [-0] | 0.586 ±1% **+5** [+0] | 0.628 ±0% **+13** [+7] |
| thin_f | 0.487 ±1% | 0.500 ±0% **+3** | 0.510 ±0% **+5** [+2] | 0.510 ±1% **+5** [+2] | 0.509 ±1% **+4** [+2] | 0.505 ±1% **+4** [+1] | 0.508 ±3% **+4** [+2] |
| RT feed | 1.311 ±1% | 1.443 ±1% **+10** | 1.443 ±1% **+10** [+0] | 1.468 ±2% **+12** [+2] | 1.431 ±0% **+9** [-1] | 1.496 ±2% **+14** [+4] | 1.369 ±2% **+4** [-5] |
| feed_x1 | 1.106 ±5% | 1.165 ±4% **+5** | 1.132 ±2% **+2** [-3] | 1.153 ±0% **+4** [-1] | 1.118 ±0% **+1** [-4] | 1.095 ±0% **-1** [-6] | 0.996 ±1% **-10** [-15] |

## Verdict per key (> 5 % cost = opt-in)

| key | cost vs defaults (the bracketed column) | verdict |
|---|---|---|
| clean_gas_exact (default true) | turning it off changes the time by -4..+3 %, all within noise | no measurable cost: **keep default** |
| clean_pc_factor (default true) | turning it off changes the time by -4..+5 %. On the x1line deck feed_x1 it is -1 % | no measurable cost on CPU: **keep default** (it is also no gain on feed_x1 at 16 cycles) |
| clean_rot_ring = live (default on rotating decks) | turning it off changes the time by -6 % on flz_q, inside the ±4 % spread of that arm. thin_f, the other rotating deck, gets +2 % | within noise: **keep default**. Recheck on a longer flz run if you want to be sure |
| clean_yres | already opt-in on 0678599eb; set false only in the all-off arm | - |
| clean_fv_joint = true (opt-in on 0678599eb) | +7 % on cs_thin1 and +6 % on cs_closed (spread 3 %). -15 % on cs_feed1, -5 % on RT feed, -15 % on feed_x1 | costs > 5 % on cs_thin1: **opt-in is right**. Note that it is a clear speed-up on the feed decks |
| all off | -5..+4 % vs defaults, and +4..+26 % vs push 15 | the slowdown vs push 15 sits outside the keys |

## A arm: 0678599eb with all new defaults off vs push 15 on the reg harness
- Run: arm OFFH (S/athenaH_off) in job 230835, after the timing on the same node. It was compared with push 15 PREV by tools/cmpx.py
  (data blocks, bin headers ignored), and the wb harness was run too.
- **Bitwise to push 15** on every reg run, except the 5 rst runs below:
  - cc/suite (5 runs);
  - cs uni, closed, thin1, thin1mb, thin1mbmr, feed1;
  - fb flz q and Q. flz q has 5 changed CLEAN log lines but bitwise data, as in the earlier OFFC arm;
  - rhse (steady start, force_reference = none);
  - the thin decks f, p, r, s81 and s88 (s81 and s88 use a steady start in sph mode);
  - dif, lat, pulse0;
  - rw (4 runs);
  - rst hcyl and rst thin1.
- **Not bitwise**: the 5 rst runs that use `rad_m1/force_reference=clean`. All have W 0, DT 0, F 0.

  | run | maxrel(bin) | note |
  |---|---|---|
  | rst sp3mb | 1.24e-5 | steady start + predict (5e33f056a), and force_reference = clean (4ac0c0e49) |
  | rst sp3mbmr | 1.24e-5 | same deck, 2 ranks |
  | rst feed | 8.37e-3 | force_reference = clean: the 4ac0c0e49 reference change (your fref-feed-8e-3 NOTE) |
  | rst feedmr | 8.37e-3 | same deck, 2 ranks |
  | rst feedx1 | 8.37e-3 | same, x1line deck |

- wb harness (steady-start rhse decks, run with force_reference = clean), max v/c_s, push 15 -> all-off 0678599eb:
  - fv1: 8.39985e-08, unchanged;
  - fv: 6.87442e-08 -> 6.88374e-08;
  - lld: 5.57512e-08 -> 5.57380e-08;
  - cart: 1.85442e-07 -> 1.85468e-07.

  These are also steady start + force_reference = clean, so the change is expected. CLEAN WARNING is 0 everywhere.
- So the only decks that are not bitwise are the ones the two intended fixes reach: steady start + predict (5e33f056a) and
  force_reference = clean (4ac0c0e49). We did not split sp3 between the two fixes. The timing decks RT feed and feed_x1 also run with
  force_reference = clean, so on 0678599eb their trajectory differs slightly from push 15. The 5e33f056a table below, which does not have
  4ac0c0e49, shows the same extra cost there, so that change does not explain it.

## Partial 5e33f056a run (job 230811, cancelled at the head move; s/cycle)
Round 1 was complete. Round 2 has only off and pc_factor, so the other columns are 1 run each.

| deck | push 15 | 5e33f056a defaults | gas_exact=false | fv_joint=false | pc_factor=false | all off |
|---|---|---|---|---|---|---|
| flz_q | 1.778 | 1.956 **+10** | 1.885 **+6** [-4] | 1.872 **+5** [-4] | 1.905 ±1% **+7** [-3] | 1.954 ±2% **+10** [-0] |
| flz_Q | 0.278 | 0.295 **+6** | 0.307 **+11** [+4] | 0.291 **+5** [-1] | 0.291 ±0% **+5** [-1] | 0.300 ±1% **+8** [+2] |
| cs_closed | 1.337 | 1.798 **+35** | 1.629 **+22** [-9] | 1.545 **+16** [-14] | 1.609 ±2% **+20** [-11] | 1.625 ±3% **+22** [-10] |
| cs_feed1 | 1.075 | 1.029 **-4** | 1.073 **-0** [+4] | 1.161 **+8** [+13] | 1.073 ±4% **-0** [+4] | 1.249 ±5% **+16** [+21] |
| cs_thin1 | 0.553 | 0.634 **+15** | 0.633 **+14** [-0] | 0.590 **+7** [-7] | 0.630 ±1% **+14** [-1] | 0.582 ±1% **+5** [-8] |
| thin_f | 0.498 | 0.503 **+1** | 0.532 **+7** [+6] | 0.497 **-0** [-1] | 0.514 ±0% **+3** [+2] | 0.509 ±2% **+2** [+1] |
| RT feed | 1.311 | 1.336 **+2** | 1.376 **+5** [+3] | 1.447 **+10** [+8] | 1.374 ±1% **+5** [+3] | 1.483 ±2% **+13** [+11] |
| feed_x1 | 1.055 | 1.019 **-3** | 0.989 **-6** [-3] | 1.109 **+5** [+9] | 1.030 ±1% **-2** [+1] | 1.162 ±2% **+10** [+14] |

- On this head, fv_joint was a default. Turning it off cut cs_thin1 by 7 % and made the feed decks 8..13 % slower: the same pattern as the
  fv_joint = true arm on 0678599eb.
- The 1-run cs_closed default (1.798) looks like an outlier: 0678599eb gives 1.582 ±1 % with the same code apart from fv_joint.

## Files (Orion)
All under /orion/ptmp/jinma/clean_1115/push16_reg:
- S/jobH.sbatch, S/nodeH.log and S/p16H.230835.out for job 230835;
- S/jobG.sbatch, S/nodeG_5e33_partial.log and S/p16G.230811.out for job 230811;
- tables: S/timH.txt, S/timG.txt and S/tools/timMD.py;
- A arm: S/cmpx_OFFH_PREV.txt and S/wb_OFFH.log;
- wrappers: S/athenaG_* and S/athenaH_*.

Restart dumps were deleted.
