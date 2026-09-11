---
name: rt-semi-implicit-changes-dhj-answer
description: Orion's semi-implicit two-stream RT source (048dff30, merged 2026-09-09) is NOT neutral on the cs dhj configuration despite its commit message; KE differs 3e-4 at cycle 1 and 3 % at cycle 20. Flag problem/rt_semi_implicit added; ablation arms must set it false
metadata:
  type: project
---

A/B 2026-09-09 (viper, CPU, 20 cycles, cs_ablate/ctl input at 16x16 panels): pre-merge ca616d1a
vs merged cb8c6c2f are NOT bitwise identical. Mass exact, tot-E 3e-6, 1-mom 7e-5, KE -3 %
by cycle 20 (-3.6e-4 at cycle 1). A third build at 1d99e22e (the RT extraction into
src/utils/{correlated_k,two_stream_rt,atm_column}.hpp) is byte-identical to pre-merge, and
ctl/nocond diverge identically, so the culprit is 048dff30 alone (semi-implicit
de = (src/lambda)(1-exp(-lambda dt)), lambda = 4E/e, unconditional). Its commit claims
"bitwise identical over 200 cycles" on deep_hot_jupiter_rt -- false for rt_ck + ck_pcut_bar 10.
Cost +4.5 % on CPU. Every other merged change is default-off or diagnostic; 6ef9492d (cycle-0
conduction dt) does not bite here; 5c0b98e4 (WB restart cache) affects restarts only.

**Why it matters:** the cs vertex dt-collapse ablation needs a ctl arm on the physics that
collapsed. The first restart arms (cs_ablate_r/ctl, nocond, merged binary) did NOT collapse
where the original did -- uninterpretable until the flag arms and the old-binary controls
(cs_ablate_r/old_cache1, h_old_cache1) are read.

**How to apply:** `problem/rt_semi_implicit` (default true = orion's behaviour; false = the
pre-048dff30 explicit source). Set false in every ablation / reproduction input. Whether the
semi-implicit source is wanted for dhj production is a separate decision (it is a stability
improvement for the optically thin top); do an A/B at production length before adopting.
See [[cs-vertex-dt-collapse-0907-defaults]], [[inflight-2026-09-09-viper]].

**Residual audit 2026-09-10 (CPU, flag false, pre ca616d1a vs 42c0a6a2): BIT-IDENTICAL** on
MHD twin 16x16 40 cycles, hydro FULL 128x32x32 15 cycles, hydro 16x16 200 cycles. The only
added default is `<hydro>/rad_blend_radial = true` (= old algebra). Two conduction changes can
move a number only when conduction OWNS dt: the cycle-0 weight build (pre-merge cycle 0 ran
with NO conduction constraint; harmless here, but a run with rad_kappa_fac=1e8 shows 1163x)
and the flux-limited diffusivity in NewTimeStep (keff <= kappa, so dt can only rise). In the dhj
runs the radial hydro CFL owns dt. Latent: mesh.cpp's "dt COLLAPSE" owner report reads only
phydro, so it never prints for an MHD run.

**REINTERPRETED 2026-09-10 (GPU bisect + per-cycle A/B).** The "3 % KE" was CHAOTIC AMPLIFICATION
of a ROUND-OFF-level change: at cycle 1 the 1-KE difference is 6.6e24 erg against tot-E 6.6e38
(1e-14), mass/momenta/tot-E/dt unchanged to 1e-6; it grows geometrically to 10 % by cycle 50
(the cs 1-ulp amplification, [[cs-ulp-amplification]]). Same magnitude on CPU (flag on vs off)
and on GPU (old f6f0d6f4 vs new). GPU bisect (HIP builds of ca616d1a, 1d99e22e, 048dff30):
viper's guards and the RT extraction are bit-clean; 048dff30 alone reproduces the new binary.
On GPU the flag does NOT restore bitwise identity because 048dff30's rt_Em accumulation
inside the chain sweep kernels runs unconditionally and changes hipcc's FMA contraction of
the neighbouring flux sums. And on GPU semi-implicit ON vs OFF (cb8c6c2f vs ae300379) are
byte-identical over 6 rot: the branch is effectively inactive on dhj (lambda*dt <= 1e-4).
GPU runs ARE deterministic across jobs and nodes (old binary reproduces the original lhllc
history byte-for-byte). Bottom line: orion's merge is semantically neutral on dhj; nothing
orion did wrongly impacted the runs; the ORIGINAL collapses cannot be reproduced by any
binary that differs at round-off, so the collapse is a chaos-selected rare event.

**CORRECTION 2026-09-10 02:00:** the "rt_Em FMA contraction" explanation above is RETRACTED.
The GPU binary never read the flag ([[gpu-build-race-flag-not-compiled]]); every GPU arm ran
semi-implicit ON, so "old vs new on GPU" = OFF vs ON, and the flag has NOT been tested on GPU.
What stands: the semi-implicit source moves 1-KE by 3.65e-4 at cycle 1 (1e-14 of tot-E), CPU
and GPU alike, and the cs flow amplifies that to 10 % in 50 cycles. Semantically small,
chaotically decisive for reproducing a given realization. GPU determinism across jobs/nodes
stands (old binary reproduced the original lhllc history byte-for-byte).
