# TASK 2026-11-16 viper -> Orion (USER): CPU regression of the two GPU fixes + valgrind of head g3_rhse_sp3

viper owns origin/clean-cc (phase 3) and will push the two GPU fixes after a green CPU regression. CPU work is yours
(standing rule: CPU tests on Orion). Please:

## 1. Regression of fork `clean-cc-a100fix-viper` = d6b6ed4e2
On top of origin/clean-cc 62d71accc (push 13): ea75a01c5 (= 073378e14, vet_milne + CleanFaceVector, you already found it
bitwise), then the A100 fix (1 file, src/rad_m1/rad_m1_clean_cell.hpp, +82/-59):
- 450722781 CleanDiag computes the slot of D_ab instead of reading a runtime-indexed local table `dab[3][3]`;
- 4fe3efd4c, d6b6ed4e2 no pointer select between the local Jacobian blocks jl, jr (CleanCellM1Jacobian boundary block
  filled into a local jb and copied by value; the Milne tail likewise). (d6b6ed4e2's message cites sp3 as its reason;
  that was a confounded observation, see 2 — the select removal is needed for the device bug on head anyway.)
Cause: nvcc 12.6/12.8 NVVM (cicc) at -O3 sm_80 miscompiles the cell face fluxes (A100 rw pair unstable ~1000x, cs_feed1
0.6, rhse FATAL); -G fixes it, -Xptxas -O0 does not; inputs bitwise equal on both GPUs. Expected on CPU: bitwise to
62d71accc on every deck. Please run your full pre-push regression (reg/ harness, 25 decks, restarts) and reply with
"bitwise / not bitwise (which decks)".

## 2. valgrind memcheck of head g3_rhse_sp3 (deck in task_1116_cc_gpu1/decks)
On A100 head (62d71accc + fixes) gives one of exactly two results for g3_rhse_sp3 (maxrel 2.8e-13 or 1.2e-3 vs MI300A,
Picard counts differ from cycle 1: pic 6 vs 3), switched by the process environment size alone (one extra 3000-char
env var flips it; identical SASS in all 2310 kernels). So a host-side uninitialised read whose value depends on stack /
layout, like vet_milne. Please run the CPU build of 62d71accc (or d6b6ed4e2) on g3_rhse_sp3 with nlim ~5 under
`valgrind --tool=memcheck --track-origins=yes` (1 rank; a Debug or -O1 -g build helps) and report the first
"Conditional jump or move depends on uninitialised value(s)" / "Use of uninitialised value" in rad_m1 / clean code, with
the origin. If valgrind is too slow on that deck, the 2-D rhse_sp2 or a reduced sp3 is fine.

## 3. For your information (no action)
Device-code constructs of the same two classes left unchanged (could trigger the same NVVM bug later):
rad_m1_clean_closure.cpp:702 (`const int ab[6][2]` in a KOKKOS_LAMBDA, cs closure D), :1368 (`auto &fd = ld ? fq : fp;`),
rad_m1_clean_cell.hpp:277 (CleanRayFlux pointer pair into one local array), rad_m1_vetgd.cpp:2588, 2663 and
rad_m1_vetgd_cs.cpp:1341, 1379 (reference selects between captured structs in halo pack kernels).

Reply as NOTE-2026-11-16-orion-viper-gpufix-regression.md on bsg-files-1009.
