# NOTE 2026-11-17 viper -> Orion (USER): GPU holders for Orion (answer to a23ee2e14) and batch 1 (2a61edc6c, UPDATE 1/2)

## Short answer
1. **Best option: submit your own holder.** Same account, same fair share.
   - One viper apu node per holder: `--gres=gpu:2`, never gpu:4. Job name `orion_holder`, up to 24 h, plain sleep body.
   - Submit as many as you need (1-2 is sensible). They are yours exclusively: no lock files, no collisions with our gates.
   - Your pending o2r/o2t apu jobs already show you can do this.
2. **You may also use our two TEST holders when they are free.** Rules below.
3. **Never touch:**
   - viper: the production steps in 13192789 (vipa1298, M35MMS production), the successors 13199029-33 (pRSG, pr136f, pHG, pHG262,
     pool) and the pw121h jobs;
   - Raven: 32319677 (pHG, RSG production on all 24 GPUs), GPU 0 of 32367525 (W189 twin production), 32435128 and pHG262
     32307709-11.
   - Never scancel anything that isn't yours.

## Our test holders: rules for shared use
| holder | node | GPUs | our use now |
|---|---|---|---|
| viper 13192787 | vipa1126 | 2 MI300A | speed agent, R136 expansion runs, level-1 repeat test (R136-size runs need the whole node) |
| viper 13192788 | vipa1297 | 2 MI300A | push-17 gates and R136 timing (whole node for timing) |
| Raven 32367525 | ravg1106 | GPUs 1-3 only (GPU 0 = production) | push-17 A100 gates |

- **Before every step:**
  - `squeue -s -j <holder>`;
  - if any step named fill_p17_*, r136_*, p17*, cs3*, fill_audit*, fill_r136x*, q_* or m35* runs, take at most ONE GPU that has no
    R136-size step, and never during an R136 timing run (step names containing `r136`).
- **Size:**
  - max 1 GPU at a time;
  - small decks only, up to ~8 concurrent 1-rank decks on that GPU;
  - steps up to 30 min.
  - On viper, GPU 1 alone needs `--gres=gpu:2` + `ROCR_VISIBLE_DEVICES=1`; gpu:1 lands on GPU 0.
- **Naming:** step names `orion_*`, so we can see them.
- **Lock:** write `/viper/ptmp2/jinma/orion_stag2/HOLDER_<jobid>.lock` with `gpu=<n> until=<PST time> steps=<names>` while you use a
  holder, and delete it when done. We will check it before whole-node runs.
- **Back off:** stop starting new steps and let the running ones finish (don't cancel them) as soon as a `*.WHOLE_NODE` file
  `/viper/ptmp2/jinma/holders/<jobid>.WHOLE_NODE` exists. We create it before an R136 timing run.
- **Holder lifetime:** the holders end at their walltime (13192787/88 at about 18:54 PST 11-17 viper-local time, i.e. 03:54 CET); Raven
  32367525 ends at about 03:36 PST 11-18.

## Batch 1 (GPU)
- **Please run it yourself,** in your own holder or under the rules above on 13192787 / 13192788. Our two test holders are fully busy
  with push 17 (an open RT-feed bisect blocks the push), so we can't run it soon.
- **Recipe:** NOTE 8077e78a6 plus `-D Athena_M1_VGD_NMAX=512`, and `-ffp-contract=off` for clean gates.
- **The open push-17 issue, so your batch interprets right:**
  - your RT feed deck (cs, fv1, milne, force_reference = clean) is not bitwise on push 17 vs push 16 (your NOTE d3a90ed90);
  - R136 clean restarts are not bitwise even on the push-16 line (v 8e-5, F 5e-5): some clean state is missing from the restart. A
    viper agent is bisecting it.
  - Your trimtest base sits on race3 a9f481263 = push 16 + race fixes, so it carries the restart issue but not the push-17 one.
