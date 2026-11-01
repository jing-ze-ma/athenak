# NOTE 2026-11-01 Caltech -> viper (USER): rkg_halo_ien_only v2 vs rkg_lean on r1m: lean stays (H is +43 % vs L at R15, +18 % at R14). On r1m every slot is exchanged, because valfven_noinfall_fac writes the momenta.

Answers TASK bf6ae9e3.

## Setup
- Binary: rkghalo-test-1101 @ 3503c166, athena_gpu_dhj_rkghalo2_3503c166_fmadoff, md5 5c4892fa86ff71823cfa940652c7e66c. Incremental build in our 0594b145 tree.
- Start: keep/off_0p50.rst (developed 10 G r1m, t = 0.5 rot, cycle 8774). r1m_trc deck, production keys (rkg_dt_safety 0.5, outflow1, fac 4, ...).
- 1 rank on ONE dedicated H200 (hpc-sm-01-05 GPU 3, nothing else on it), 300 cycles per run.
- Interleaved L, N, H x 3, first at max_eta 1e15 (R15), then at 1e14 (R14).
- s/cycle is taken from the run.log stamps at cycles 8800-9000 (start-up excluded). The repeats agree to <= 0.2 %.

| max_eta | arm | s/cycle (3 repeats) | mean | vs L | RKG stages |
|---|---|---|---|---|---|
| 1e15 | L (rkg_lean) | 0.0672 / 0.0672 / 0.0672 | 0.0672 | - | 16 |
| 1e15 | N (non-lean) | 0.0947 / 0.0947 / 0.0947 | 0.0947 | +40.9 % | 16 |
| 1e15 | H (non-lean + rkg_halo_ien_only) | 0.0962 / 0.0962 / 0.0962 | 0.0962 | **+43.1 %** | 16 |
| 1e14 | L | 0.0431 / 0.0430 / 0.0430 | 0.0430 | - | 5 |
| 1e14 | N | 0.0505 x 3 | 0.0505 | +17.3 % | 5 |
| 1e14 | H | 0.0509 / 0.0509 / 0.0508 | 0.0509 | **+18.2 %** | 5 |

H histogram lines:
- R15: `rkg_halo_ien_only: 4494 super-stages >= 2, 4494 with every slot exchanged; per slot: 0:4494 1:4494 2:4494 3:4494 4:4494`
- R14: `rkg_halo_ien_only: 1199 super-stages >= 2, 1199 with every slot exchanged; per slot: 0:1199 1:1199 2:1199 3:1199 4:1199`

## Why every slot on r1m: valfven_noinfall_fac
- The same H run with `mhd/valfven_noinfall_fac=0`, 50 cycles: `744 super-stages >= 2, 0 with every slot exchanged; per slot: 0:744 4:744`. That is {IDN, IEN} only, as you saw.
- So the r1m production key fac 4 (no inward radial momentum in floor-propped cells) makes the ConToPrim write the momenta on every stage.
- Even at fac 0 the 1-rank saving is small (one repeat each, 300 cycles, same GPU):

| fac 0, max_eta 1e15 | s/cycle | vs N |
|---|---|---|
| L | 0.0650 | -27.9 % |
| N | 0.0901 | - |
| H ({IDN, IEN}) | 0.0889 | -1.3 % |

- On one rank the halo exchange is device-local, so the slot count hardly matters. Viper's -17 % (s10, MI300A) presumably comes from multi-rank MPI exchanges. Our r1m jobs are 1 rank per H200.

## Verdict for the USER's rule ("H within ~10 % of L -> recommend H")
- **H is not within 10 % of L**: +43 % at R15, +18 % at R14.
- On r1m, H is even slightly slower than N (+1.6 % / +0.8 %), because fac 4 forces the full exchange and the mask reduction is extra work.
- **Lean stays for r1m.** H would only pay on multi-rank runs without noinfall.
- I skipped the optional L vs H 0.5-rot difference: H fails on cost, so the noise-level question does not decide anything.

## Side note on our G4 numbers (NOTE ...-gates-results, 0594b145)
- The G4 cost runs on hpc-sm-02-03 GPU 3 gave ~0.26 s/cycle for non-lean, vs 0.095 here. That node's CPUs were shared with 5 other runs and a build, and these runs are launch-bound.
- The G4 off/on RATIO (+1.7 %, interleaved) still holds and matches H vs N here (+1.6 %). The absolute numbers in that NOTE are not representative.

Run dirs: /resnick/groups/carnegie_poc/jingze/trc_1101/runs/lnh15_*, lnh14_*, lnh15_H_fac0, lnhf0_* (script gates_1101/lnh.sh).
