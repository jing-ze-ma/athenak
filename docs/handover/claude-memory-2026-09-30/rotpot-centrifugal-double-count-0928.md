---
name: rotpot-centrifugal-double-count-0928
description: ROOT CAUSE of w121 100-bar warming / 1-10 bar cooling: horizontal centrifugal energy work double counted when rot_potential + etotgrav (dhj pgen SourceFunc); fix branch dhj-rotpot-energy in progress
metadata:
  type: project
---
Found 09-28 (deepdrift_0927/RESULTS.md, verified by main session in code): etotgrav adds rho*v*Phi on x1/x2/x3 faces
with Phi = TotPotAt (gravity + centrifugal -Omega^2 (r sin th)^2/2), so the full centrifugal work is already in the
energy; deep_hot_jupiter_rt.cpp SourceFunc removes only the RADIAL centrifugal work (cs ~3334, sp ~3262) -> the
horizontal part is added twice: ~-1e31 erg/s in 0.1-10 bar, ~+1.1e31 below 10 bar (net ~1e29 = the "unexplained"
global residual). Fix test S2: 100-bar drift +28 -> +0.2 K/10 rot. Deep CZ destroyed (sub-adiabatic to the wall).
Affects every rot_potential+etotgrav run: w1x, w10x, likely prod3/prod4 (deep profiles after rot ~20 suspect).
User approved 09-28 ~01:00: implement + gates (adiabatic E conservation, static balanced state, momentum-side check
of the x2/x3 WB reconstruction, bitwise when off, 10-rot w1x restart) on branch dhj-rotpot-energy
(/viper/ptmp2/jinma/rotpot_0928); main session merges after review. Then decide restart vs clean start of w121
and whether the fixed-entropy bottom BC is still needed. Related: [[w121-bottom-wall-energy-0927]],
[[session-0927-night-agents]].

**09-28 ~05:30 USER DECISIONS:** merge dhj-rotpot-energy (tst dhj_ck_mpicpu -> tolerance compare; base already
rank-dependent 1e-16..3e-14) and cs-resist-stretch via the merge coordinator (/viper/ptmp2/jinma/mergecoord_0928,
order: box-wall-truephi, energy-guards, dhj-rotpot-energy, cs-resist-stretch). **FRESH START of w1x and w10x** from
the final merged tip (NOT a restart: the bug froze a sub-adiabatic layer below ~60 bar, diffusion time ~1e6 rot).
10-rot fixed-code check: 100-bar drift +28 -> -1.7 K/10 rot; old 1-10 bar temperatures were ~100-190 K too low.

**09-28 ~08:20 USER: WASP-121b fresh start ON HOLD** ("let's not do the wasp121b for now. make sure the code works").
Tip after the night's merges: 08216395 (m1-mhd2, box-wall-truephi, energy-guards, dhj-rotpot-energy + tst tolerance,
cs-resist-stretch, ck-determinism). ck-lin2-hip V1 NOT merged (shifts GPU results vs tip; only restores pre-ck-lin2
rounding) -- recommended to drop. Running: validate_0928 (full tst + production smokes), audit_m1_0928 / audit_dhj_0928
(Caltech TASK-2026-09-28-speed-default-audit: keep default only if >= ~3 % on MI300A AND <= 1.2x noise),
ckrev_0928 (HIP noise gate + structural review of ck-next/ck-lin2; decides keep vs revert: on MI300A ck-lin2 -2.7 %,
ck-next +1.5 %), cknewton_0928 (10x Newton 133/500 non-converged), topleak_0928 (top-shell energy residual 1-3e29).

**BOTTOM-BC VERDICT 09-28 (botbc_0928, fresh 1x on 11c9a5be, 25 rot in 64 min on 2 GPUs = ~24 rot/h):** no fixed-entropy
BC needed. Deep stays on the starting adiabat (s at 100/300 bar -5.685/-5.690 vs IC -5.689/-5.685; old run drifted
+0.1-0.4). 50-230 bar weakly stable (nabla-nabla_ad -0.001..-0.004), only 300-475 bar convects; wall exact
(Mdot 0, Etot_bot = Lrad_bot). Thermal time below 10 bar ~8e5 rot -> the IC deep adiabat IS the interior entropy:
choose it deliberately. Recheck s(100)-s(300) at ~rot 100 (+0.001/10 rot). T(1 bar) 2006 -> 1801 K (upper layers still
relaxing).
**USER 09-28 evening: KEEP THE CURRENT IC ADIABAT** (the deep adiabat of the 1-D ck RCE IC with the input's T_int)
for the WASP-121b fresh start; it sets the interior entropy for any feasible run length.
Fresh-start input keys: 1x ck_impl_maxit 16; 10x ck_impl_maxit 24 + ck_impl_dtmax 0.25; ck_impl_xstep 8 named for
spin-up; flux_hst_floor optional (merged). Tip with everything: 6d690e09. Launch still needs the user's go.
