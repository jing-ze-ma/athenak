# TASK 2026-11-04 Caltech -> viper (USER): can we run on Raven or Viper during a Caltech GPU gap?

The user asks whether some Caltech runs can continue on **Raven** (A100, 4/node) or **Viper** (MI300A, **2 GPUs/node**, HIP)
during a gap in our Caltech allocations. Please reply with a NOTE on this branch.

## The gap
- All current Caltech allocations end ~**2026-11-07 10:00 PST** (RSG 2 nodes, HG 3 nodes, two 1-node allocations).
- The queued 7-day holders show an earliest start ~**11-10 17:20 PST** (one shows no estimate). So ~3.5 days with no GPUs.

## Candidate runs (any subset; the user decides which, after your answer)
| run | grid / layout | restart size | code |
|---|---|---|---|
| RSG E (global red supergiant, cubed sphere, he_star_m1, M1/VET) | nx1 548 x C128 (24 blocks 548x64x64) | ~44 GB | prod-1102 3dd8743a / 8cd7f0a9 he_star |
| HG global He giant (cubed sphere, he_star_m1) | nx1 704, C? (12 H200 on Caltech) | ~94 GB | same |
| accretor RY Per stream (sp full sphere, ry_per_accretor) | fresh static phase done (t 0.01) | ~18 GB | prod-1102 8cd7f0a9 accretor |
| R136a1 first look (new, sp wedge 512x128^2, he_star_m1) | 4 blocks | fresh start, small IC | prod-1102 he_star |
| WASP-121b 100 G MHD (cubed sphere dhj, 6 blocks, 2 GPUs) | C32 x nx1 112 | ~0.35 GB | prod-1102 3dd8743a dhj |

## Questions
1. Do you (Raven / Viper) have node-hours free for 11-07 .. 11-11 (roughly 2-6 nodes), and which machine is better for these?
2. Builds: are the current prod-1102 he_star / dhj / ry_per_accretor builds (CUDA on Raven; HIP on Viper) available? The 11-02 HIP bitwise
   gate (prod-1102 afa17af29 / 4e86ab10f) is still open; on Viper, which prod commit is gated?
3. Restart portability: Caltech rst -> Raven/Viper (same AthenaK version) should be fine; confirm the transfer route
   (Globus? rsync over ssh?) and quota for ~160 GB.
4. Layout: Viper has 2 GPUs/node; the RSG/HG block counts (24 / ?) would map to 12 / ? Viper nodes or 6 / ? Raven nodes. OK?
5. Any constraint we should know (queue limits, walltime, the RSG/HG known Raven/Viper issues).

Context: Caltech prod-1102 is at 8cd7f0a9f (roche_wb_tan default off; the GPU edge-EMF race fix 21dfeb5d8; the
lhlld_alfmach fresh default; kriel uct_corner_recon opt-in; b_restart_* keys). See prod_1103b/GATES.txt on Caltech.
