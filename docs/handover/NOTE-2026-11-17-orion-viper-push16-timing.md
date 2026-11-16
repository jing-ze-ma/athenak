# NOTE 2026-11-17 Orion -> viper: push-16 CPU timing (follow-up to f16759de0)

## Setup
- Job 230771 ran on one dedicated p.exclusive node (orion083), with the three binaries one after another on the same node.
- Each deck ran standalone (OMP threads as in the reg harness). The figure is cpu time used / cycles, with the cycle count in brackets.
- Binaries: push 14 d6b6ed4e2 (md5 e0899b8c), push 15 dfe21a746 (0f714f25), push 16 0aee9a208 (1abe7a4a), default keys.

## Results (s/cycle)

| deck | push 14 | push 15 | push 16 | 16 vs 15 |
|---|---|---|---|---|
| flz_q | 2.359 (40) | 1.685 | 1.925 | +14 % |
| flz_Q | 0.362 (40) | 0.280 | 0.296 | +6 % |
| cs_closed | 1.910 (4) | 1.545 | 1.652 | +7 % |
| cs_feed1 | 1.309 (10) | 1.044 | 1.125 | +8 % |
| cs_thin1 | 0.660 (100) | 0.540 | 0.634 | +17 % |
| thin_f | 0.538 (150) | 0.492 | 0.507 | +3 % |
| RT feed (16 cycles) | 1.734 | 1.342 | 1.412 | +5 % |
| feed_x1 (x1line, 16 cycles) | 1.764 | 1.060 | 1.023 | -3 % |

## What it shows
- **The CPU slowdown is gone.** Push 15 dfe21a746 (host line solve, one thread per line) is now 1.1-1.7x FASTER than push 14 on CPU.
  The 1.2-3.2x slowdown we reported was for the pre-fix e3f4881dc.
- **Push 16 is 3-17 % slower than push 15** on every deck except the x1line one (-3 %).
  - Largest: cs_thin1 +17 % and flz_q +14 %.
  - Candidates: clean_fv_joint (fv1 levels over all MeshBlocks), clean_gas_exact, clean_pc_factor.
  - We did not bisect the keys. Tell us if you want a per-key timing (each new default off in turn).
