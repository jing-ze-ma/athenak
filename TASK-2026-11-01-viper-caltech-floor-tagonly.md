# TASK 2026-11-01 viper -> Caltech (USER): one tag-only floor-tracer arm (0.5 rot), to tell whether the 1.2e-5 bar shell difference is a correction or a side effect

Follow-up to your T1 (d35f4a9ad). Fast: a single arm.

## Why ~2/3 of the tag is not drained (viper analysis, C8, 300 cycles)
- At p < 1e-6 bar, 100 % of the tag is in floor-SUPPORTED cells (rho < 2 rho_fl), and tag/mass there is 0.46-0.91.
- That mass is the floor's own reservoir. As B grows, rho_fl = B^2/v_A^2 grows, and draining it just re-floors it (margin 1 showed the pump: drain 2.34e15 = re-floor 2.33e15 g/s).
- So the undrained ~63 % is mostly the top reservoir at p < 1e-6 bar, which is harmless under the accuracy rule as long as it stays there.
- The Mdot_top drop is expected: the drain removes floor-made gas before it leaves.
- **Open question: the 1.1-1.2e-5 bar shells** (rho diff 6.4e-3, 2x the fac spread). There rho/rho_fl ~6, so the drain acts. Either:
  - (A) it is the intended correction: off contains leaked floor mass there;
  - (B) it is a side effect.

## The arm
- Binary: **fork/tracers-defaults-1101 @ 8ccae6636** (includes the x_e floor-dilution fix: floor mass now inherits the cell's x_e). Or 89213fcb, for a like-for-like comparison with your T1 fl arm.
- Deck: your T1 deck. Arm **tagonly** = floor_tracer = true, floor_tracer_fmax = 1e-12 (the tag is measured but essentially not drained), xe_tracer = false. t = 0 to 0.5 rot.
- If the GPU is free, also rerun **fl** on the same binary, so fl and tagonly share a binary.

## Analysis
`docs/handover/scripts/tag_where.py` (pushed with this TASK; it needs only vis/python/bin_convert.py on the path):

```
python3 tag_where.py <fl bin at 0.5 rot> <tagonly bin at 0.5 rot>
```

Per log-p band it gives the fl/tagonly mass ratio - 1 next to the tagonly tag fraction. Please send the output.

## Decision
- If the shell rho difference at 1.1-1.2e-5 bar is about the tag fraction there -> (A). The floor-tracer default stays, and the spread test is the wrong yardstick for that shell.
- If the tag fraction there is << 6e-3 -> (B). We fix it, in this order:
  - margin 4;
  - cap the drain to the tag above the local floor-support level;
  - a per-time drain rate.
- The USER decides the default after this.
