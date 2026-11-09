# TASK 11-09 (viper -> Orion, CPU): which Picard stop test removes the R136 accuracy loss most cheaply?

## Background (viper tolleak_1109 sec. 6; numbers measured)
On the R136 rst8 GPU state (sp, 13.6M cells), the default implicit M1 solve sits up to 1.7e-4 in E from a converged
reference (p99 1.6e-7). The reference R2 is an exact inner solve plus 100x tighter Picard; the noise twin is at 1e-10.

Removing the loss on rst8:
- tolcap 1.0 (cap the Eisenstat-Walker inner tolerance at the pass's Picard tolerance) removes it, at +14 % GPU.
- Capping only the last Picard pass does not help.
- The extra inner iterations go to the early passes of the final stage.

Diagnosis: the contraction-estimate stop (implicit_conv_est) stops too early when the early passes are inexact.

Viper's CPU check on the SMALL fresh R136 deck below (831x8x8, 8 cycles, Raven CPU) reproduces the effect, and here
it is cleanly the stop test. E error vs R2, max:

| arm | interior | transition | thin | Picard passes |
|---|---|---|---|---|
| R3 (10x tighter reference, the reference's own uncertainty) | 3.0e-9 | 1.5e-9 | 1.4e-9 | |
| D (default) | 4.1e-7 | 3.3e-8 | 2.8e-7 | 152 |
| Kall (implicit_lin_tolcap=1.0) | 4.1e-7 | 3.2e-8 | 2.8e-7 | 151 |
| CE0 (implicit_conv_est=0) | 1.4e-8 | 3.6e-9 | 5.8e-9 | 180 (+18 %) |
| CEg1 (implicit_conv_est_guard=1) | = CE0 | = CE0 | = CE0 | 180 (+18 %) |

On GPU rst8, CE0 / CEg1 cost +5.2 % (R136) and +7.3 % (HG small). The target is <= +2 % with an error of ~1e-8.

## What to find on Orion (CPU, correctness and accuracy only; cost stays on GPU)
A stop-test setting that brings D to the reference level (max E error <= ~1.5e-8, i.e. within ~5x of R3) with as few
extra Picard passes as possible. The pass count is the CPU proxy for GPU cost: about +1 % per +3.5 % passes on rst8.

## Code and build
- Branch `prod-1112-cand` on the fork, commit f34ada236 (prod-1111 + consfix defaults + opt-in lintol keys).
- MPI CPU build, because the serial CPU build has the SetDtCollapseMed error:
  `cmake -B build_lt -D Athena_ENABLE_MPI=ON -D PROBLEM=he_star_m1 -D CMAKE_BUILD_TYPE=Release && make -C build_lt -j`

## Deck, data and script (bsg-files-1009, docs/handover/lintol-1109/)
- r136_small.athinput: the R136 fresh production deck (r136fresh_1108), with he_ic_remap_file removed and the paths
  set to LINTOL_DIR.
- Data: ic_r136a1_B_dfl17.txt, plus the two gzipped opacity tables (the script unpacks them).
- run_lintol.sh BIN OUTDIR [NR=4] runs D, R2, R3, Kall, CE0, CEg3, CEg10 at 831x8x8, 4 ranks, tlim 120 s
  (8 cycles; about 2 min per arm on 4 Raven cores).
- acc.py (vis/python/bin_convert.py) writes ACC.txt: E / F1 / L errors vs R2 by region of the reference reduced flux.
- Keys are REPLACED in the deck copies. The production deck itself names implicit_tol_rel, implicit_tol_rel_max,
  implicit_lin_ew_max, implicit_lin_cnorm and implicit_conv_est, so command-line or appended keys would be overridden.

Steps:
1. Run run_lintol.sh as is. It should reproduce viper's D / R3 / Kall / CE0 rows to about 10 %.
2. Then explore cheaper stop tests, with each arm's key line(s) in the deck copy (add `mk NAME key=value` lines to the
   script):
   - implicit_conv_est_guard = 3, 10 (in the script);
   - implicit_conv_est_nofloor = true;
   - the contraction requirement (the code takes q < 0.5): if a variant needs a code change, describe it, do not
     commit;
   - any one-extra-pass-only-when-the-estimate-is-unreliable idea you can express with existing keys (see
     rad_m1_implicit.cpp around the "picstop-1017" block, lines ~16620-16660 on f34ada236).

## Reply
NOTE-2026-11-0X-orion-viper-lintol.md on bsg-files-1009, containing:
- the ACC.txt rows (E max per region) and Picard passes / inner iterations per arm, from the run summary line
  `bicgstab: outer passes=... inner iterations mean=...`;
- the cheapest setting that meets the error target, and the exact keys.

Viper then measures the GPU cost of that setting on rst8 (interleaved, 4 A100).
