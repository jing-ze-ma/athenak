---
name: carnegie-hpc-usage-policy
description: Carnegie's Caltech HPC policy (from the user, 2026-09-25) -- per-quarter overuse thresholds (GPU 7,500 / 15,000 GPU-h, CPU 500k / 1M CPU-h, storage 50 / 100 TB) and the required publication acknowledgement
metadata:
  type: reference
---

Carnegie's allocation on Caltech HPC. Policy text pasted by the user on 2026-09-25 ("Caltech HPC Policies", last updated 2023-10-11). This is the budget that matters, not Caltech's $ rates in [[caltech-cluster-facts]].

**Overuse thresholds, per quarter:**

| resource | Tier 1 | Tier 2 |
|---|---|---|
| GPU | 7,500 GPU-h | 15,000 GPU-h |
| CPU | 500,000 CPU-h | 1,000,000 CPU-h |
| storage | 50 TB | 100 TB |

- Tier 1: SciComp / Carnegie HPC may start a dialog and ask you to curtail use.
- Tier 2: use may be curtailed for the rest of the quarter; they suggest applying to large centres instead.

**How to apply:**
- Plan campaigns in GPU-hours against 7,500/quarter. For example, 4 arms x 1 H200 x 24 h = 96 GPU-h/day, which reaches Tier 1 in about 78 days.
- Report the projected GPU-h per quarter before proposing production.
- Keep the group directory well under 50 TB.

**Publication acknowledgement (required):** "The computations presented here were conducted through Carnegie's partnership in the Resnick High Performance Computing Center, a facility supported by Resnick Sustainability Institute at the California Institute of Technology"
