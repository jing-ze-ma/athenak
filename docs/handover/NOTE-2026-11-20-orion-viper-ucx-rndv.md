# NOTE 2026-11-20 orion -> viper (USER-approved): on MI300A the whole XENV slowdown is UCX_RNDV_THRESH=inf; ^rocm_ipc alone is bitwise and ~36 % faster

Measured on the pRSG loan 13255123 (vipa[1180,1182,1184,1185]), RSG restart rsg.01510 (cycle 684504), 2 nodes x 2 ranks;
HG hegiant.00184 on 4 nodes. Binaries: ours = staggered clean 9ba4919a7 (key off = e970303f8 + 8fef68815);
prod = athena_he_gpu72_q_864d120d_p1109 (RSG). All runs 0 CLEAN WARNING, rc 0.

**Speed (mid s/cycle, 5 cycles; ours unless noted)**

| UCX env | RSG | HG |
|---|---|---|
| default | 7.03 | 3.78 |
| UCX_TLS=^rocm_ipc only | 7.70-7.73 | - |
| UCX_RNDV_THRESH=inf only | 12.03-12.26 | - |
| full XENV (both) | 11.73-12.07 | 6.96 |
| prod RSG, default / XENV | 2.47 / 3.44 | - |

Cause: RNDV_THRESH=inf sends the MB-size halo messages (~1100 per cycle in our clean solver) as eager bounce-copy
fragments, ~5 ms each, whether the buffer is on host or device. ^rocm_ipc alone costs ~+0.7 s/cycle.

**Determinism with UCX_TLS=^rocm_ipc only (no RNDV_THRESH=inf), 20 cycles:**
- ours x2 on different node pairs: hst md5 identical (hydro 8b4e42e3, user acf44f70); default UCX gives the same md5.
- prod x2 on different node pairs: hst md5 identical (hydro a18ad43b, user c6183876).
- ours restart 1 cycle + rst + 4 vs straight 5: final dump and hst BITWISE.
- 5-cycle runs: XENV, RNDV-only, ^rocm_ipc-only and default all give bitwise-identical hst.

**Recommendation:** keep `UCX_TLS=^rocm_ipc`, drop `UCX_RNDV_THRESH=inf`. Before switching production, please confirm
with your own longer repeat + restart-continuation runs (the 2026-10-21 non-determinism that led to XENV,
rsg-rst-fix-1021, may have been a send-buffer race that eager sends mask; our runs above show none). Expected gain:
~-36 % for ours on RSG; prod also slowed under XENV (2.47 -> 3.44 mid), so it should gain similarly.

Side note: 5-cycle medians exaggerate ours/prod — over 20 cycles (^rocm_ipc) ours mean 5.61 vs prod 5.65 s/cycle;
by median 5.36 vs 3.03 (prod has periodic 10-22 s spikes, ours has a 4-cycle restart transient).

Our steps on 13255123 continue until 18:30 PST at the latest; a done NOTE follows.
