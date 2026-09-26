---
name: next-prod-ck-c2
description: User decision 09-24 -- the next dhj production uses implicit ck T4 with the c2 lever combination (ck-fast merged 90b01f2f)
metadata:
  node_type: memory
  type: project
  originSessionId: 05818bc0-259d-4b2f-b317-f82d1e083443
  modified: 2026-09-23T23:44:57.891Z
---

User 09-24: "use c2 levers for the next production". Applies to the next deep_hot_jupiter (cs) production launch, not the
running prod4 / hyd4 (untouched). ck-fast merged into rt-integration as 90b01f2f (levers default OFF, so the input must set them).

Keys (tests_ck_implicit/fast/gpu/common.sh, T4 + c2):
T4: problem/ck_implicit=true ck_impl_arat=1e30 ck_impl_frozen_op=true ck_impl_lin=true ck_impl_lin_thr=1 ck_impl_fuse=true ck_impl_jac_lin=true ck_impl_cvsec=true
c2: problem/ck_impl_once=true ck_impl_xstep=2 ck_impl_jreuse=0.2 ck_impl_pred=true ck_impl_pred_fac=0.5 ck_impl_nosync=true

**Why:** A/B 11955220 on the production state: c2 0 non-converged, 3.07 passes, day rms 3e-4 (= lever 1 error), kinks 716
vs t4 737; RT cost 2.03x hydro (T4 alone 7.24x). c4 marginal (4.56 passes), c8 fails (103 non-conv).

**How to apply:** copy the keys from common.sh into the production athinput when the next production is set up; check
the ck_sweep_form = 1 + ck_implicit compatibility (prod4 header says sweep_form 1 "refuses ck_implicit"; the weak-scaling
agent in /viper/ptmp2/jinma/cs_weak_0924 is checking it). Related: [[ck-x1-one-block-deferred]].

**Update 09-24 (user):** next production = T4 + c2 + problem/ck_impl_every=4 (no guard; ck-cadence merged e89954e2;
RT 0.94x hydro, kinks 733 vs c2 716). Blockers before launch: (1) MHD task list never calls user_split_func, so c2
once/cadence silently skip the RT on MHD (branch mhd-split, agent running); (2) ck restart not bitwise, cross-call Newton
+ cadence state not saved (branch ck-restart; its PART 2 = 3-rotation A/B every=4 vs c2 on hydro); (3) nx1 > 136 dies on
MI300A with scratch OUT_OF_RESOURCES (branch ck-scratch). sweep_form 1 + ck_implicit works (the prod4 header note is stale).

**Update 09-24 evening:** ALL BLOCKERS CLEARED. mhd-split merged 0f196ca4, ck-scratch e28fb4f4, ck-restart 7c16ff86 (restarts bitwise incl. MHD c2+every4). 3-rot A/B: e4-c2 equals c2p-c2 noise (chaotic), kinks +3.5 % over noise run; verdict every=4 OK. Weak scaling 1->16 nodes eff 1.02. Next production setup needs only the user go-ahead.

**Production goal (user 09-24, final):** 300 ROTATIONS (P_rot = 3.05e5 s = 3.53 d, so ~9.2e7 s = ~1060 Earth days; the user also calls it "1000 days"). Estimate for the 1024-around-equator MHD grid on 16 apu nodes: ~25 min/rotation at dt ~18 s -> ~125 h (~5.2 days of 16 nodes, ~6 chained 24-h links) if dt holds; longer if dt shrinks.

**User 09-26: production adds ck_impl_kkt_row = true (in wasp121_0925 x1 and 10x inputs) and ck_impl_xstep = 8 (sparc_0925/run.sub C2; was 2, a no-op).** Needs a binary >= 427f9f58. Fresh starts only; the running sponge arms keep their old keys (binary e3a8442e). Backups: *.bak0926.
