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

**User 09-27: production spin-up WITHOUT the bottom sponge** (top sponge on). Set up and smoke-tested in /viper/ptmp2/jinma/w121prod_0927; the chain is not launched until the user picks nodes and start.

**User 09-27 (updated: run 1x and 10x IN PARALLEL): the 10x WASP-121b run with the SAME updated setup, only 1x -> 10x swapped** (ck tables ckdata10, CE table 10x, IC ic_w121.txt, and the 10x grid if it differs). NO extra 10x solver settings (drop maxit 12, rsec 20 etc.; use the same defaults as 1x). Expect some NOT-CONVERGED from stiff day-side columns (the Broyden/banded Jacobian is not built); report it, don't add settings.

**User 09-27 go:** 1x production spin-up LAUNCHED on its own node, 2 GPUs, TROT 60, chained apu links (/viper/ptmp2/jinma/w121prod_0927, binary 041fac8f md5 7ea1537a; seam defaults, kkt_row, xstep 8, beam_par, top sponge only). 10x: ALLOW the 10x solver settings (tol 1e-7, maxit 12, rsec 20); smoke first, then ask for the go; separate node, 2 GPUs. No floor switches, no ck_dif_dtau. Port the remap later.
- 09-27 00:36: 10x spin-up LAUNCHED as 11992453 (6.5 h, 2 GPUs, TROT 60, 10x solver keys; smoke showed 10-50 % NOT-CONVERGED at residual <= 4e-7). Stop it and build the denser-Jacobian fix if residuals grow well beyond 1e-6 or dt falls. 1x = 11992260 (4 h).
- 09-27 07:05: both spin-ups at rotation 60 (1x 2 h 46 min, 10x 4 h 30 min, 0 collapses). **User: run the LOW-RES grid to 300 rotations before any remap to high res.** In progress: a verbose ck check at rot 60 + a relaxation check; the continuation to 300 is prepared, submission awaits the go.
- 09-27 ~07:50: CONTINUATION to rot 300 submitted: 1x 11995741 (15 h), 10x 11995742 (23:55), same binary, from dhj.00120.rst; STOP files removed. At rot 60 the deep layer was still drifting (100 bar +13-16 K per 10 rot); 10x has 45 % NOT-CONVERGED (max resid 8e-6).
