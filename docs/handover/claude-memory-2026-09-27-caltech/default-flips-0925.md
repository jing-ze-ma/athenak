---
name: default-flips-0925
description: Handover first code task (09-25 PDT) -- ck_impl_kkt_row default true and ck_impl_xstep default 8 (every>1) MERGED + pushed 504d8e30; mg_gf Cartesian default HELD on branch default-flips f3703090 until H200 M1 timing
metadata:
  type: project
---

**Merged into rt-integration and pushed as 504d8e30 (user: "merge the two ck flips now, hold mg_gf")**
- `ck_impl_kkt_row` now defaults to true.
- `ck_impl_xstep` now defaults to 8 when `ck_impl_every > 1` and the stored operator exists. An explicit value always wins.
- Gates:
  - (a) old value explicit = old binary, bitwise;
  - (b) key absent = new value explicit, bitwise;
  - (c) problems not using the feature, bitwise;
  - restart gate passed.
- Test suite: `--cpu` 220 pass / 41 fail, identical on old and new commits (the failures already existed).
- GPU smoke on an H100: clean.

**HELD: `implicit_precond = mg_gf` as the Cartesian default** (branch `default-flips`, commit f3703090).
- At the input tolerances it differs from the old default by 4-6x the round-off spread. At tol 1e-12 it is within noise, so the gap comes from solver tolerance.
- Gain measured on MI300A only: about -12 % (1 GPU) and -8 % (2 GPUs) vs rbgs_fwd.
- Its host LU and MPI_Allreduce per application may cost more on a discrete GPU.

**How to apply:** decide on mg_gf after the Caltech M1 H200 lever timing (m1-port agent). The handover table row 4 says HELD.

**DROPPED 2026-09-26 (user):** mg_gf Cartesian default abandoned (H200: -2.3 % on 1 GPU, +33 % on 2 GPUs). Branch default-flips deleted. ck_beam_par = true added to both WASP-121b inputs (user 09-26).
