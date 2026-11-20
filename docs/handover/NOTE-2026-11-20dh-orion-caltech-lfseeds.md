# NOTE 2026-11-20dh orion -> caltech (USER): floor_pay RSG seeds (4 arms); restart dumps near the deaths?

## 1. Seeds (LF_rsg died at cycle 2910 with floor_pay, while L1_rsg_off survived 3000)
Four single-node arms in `docs/handover/caltech_runs_1121/arms_lfseed.txt`. They run on B_v202, NR=4, 3000 cycles, **concurrently, right after the current beam arms (TG*/XCP*) finish**, ahead of the G RSG pair and the paused arms:
- Seed a: LFa_rsg (floor_pay on) and L0a_rsg (off), with mesh/cs_vertex_cc_perturb=1e-15.
- Seed b: LFb_rsg and L0b_rsg, with mesh/cs_vertex_cc_perturb=1e-12.
The perturbation key is diagnostic and is read only when named.
Please **keep restart dumps every 250 cycles** for these four on Caltech (no need to push), and post progress as usual.
Note: these four with the G RSG pair make 6 single nodes. If nodes are short, run the 4 seeds first.

## 2. Restart dumps near the deaths?
Do any restart dumps exist near **LF_rsg's death (cycle 592698, around N0+2900)** or **L1_rsg_on's death (592299, around N0+2500)**? If yes:
- Please push one per run, the last before the death (or say where it is on Caltech and keep it), so we can replay the failing step with rad_m1/clean_nl_print.
- If the files are too large to push, keep them, and we will post replay arms instead.

-- orion (for the USER)
