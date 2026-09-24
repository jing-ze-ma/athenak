---
name: ck-x1-one-block-deferred
description: Lifting the one-meshblock-in-x1 RT constraint (dhj) is DEFERRED by the user 09-24; design options A (column transpose) / C (affine interface solve) recorded
metadata:
  node_type: memory
  type: project
  originSessionId: 05818bc0-259d-4b2f-b317-f82d1e083443
  modified: 2026-09-23T23:37:13.170Z
---

deep_hot_jupiter_rt.cpp:839 fatals unless mesh nx1 == meshblock nx1 (ck RT sweeps, beam tau, Newton tridiagonal are
sequential in r). User 09-24: "worry about x1 constraint later" -- not needed for the 1024-equator target at nx1 = 256
(1536 blocks of 256x16x16; ck_implicit + sweep_form 1 capped at n1 <= 264).

**Why:** multi-node plan for high-res cs dhj MHD; weak-scaling test at nx1 = 256 launched 09-24 (/viper/ptmp2/jinma/cs_weak_0924).

**How to apply:** if revisited, recommended option A = transpose/pencil gather of whole columns within a radial rank group
(kernels unchanged, bitwise with 1 x1 block, keep radial groups on one node); option C = per-block affine (I_out = a I_in + b)
+ SPIKE interface solve (removes n1 caps, rewrites every sweep form). Pipelined sweeps (B) rejected. Revisit only if
weak scaling shows memory/load-balance trouble or nx1 > 264 is wanted.
