# NOTE 2026-10-10 (viper): vgdspeed-1009 merged into rt-integration

## What merged
- fork/vgdspeed-1009 ec0757c1 is merged into rt-integration f39e26b1 (non-ff, no conflicts).
- One commit was added on top: 89693b59. It changes only the summary text (see "lin_tol print" below).
- vgdspeed-1009 touches only src/rad_m1/{rad_m1.hpp, rad_m1.cpp, rad_m1_vetgd.cpp, rad_m1_vetcol.cpp}.
- Details: /viper/ptmp2/jinma/vgdspeed_1009/VGDSPEED.md.

## Speed-ups (bitwise, on by default)

| key | default | what it does |
| --- | --- | --- |
| `rad_m1/vet_gd_halo_list_mb` | 4096 MB per rank; on a GPU capped at half the free device memory at the first build; 0 = off | Per (pass, shell), the compact halo masks become gather/scatter lists. Rebuilt only when the direction set rotates. Later sweeps do: pack, MPI, zeros on the other pass's entries, scatter. The old dense path (flag / scan / pack / expand / dense unpack) runs only when the band state is not provably the same. |
| `rad_m1/vet_gd_shell_list` | true | The shell kernel runs only over the (m, k, j, d) of the pass's branch, so half the threads no longer idle. |

Both keys are read only when they are named in the input.

Measured gains, against b2f2e897:

| machine | case | base s/cycle | new s/cycle | change |
| --- | --- | --- | --- | --- |
| Raven 4 A100 | AG Car A hr | 0.527 | 0.403 | -23.6 % |
| Raven 4 A100 | BSG hr reduced | 0.391 | 0.299 | -23.6 % |
| viper 2 MI300A | AG Car A hr | 1.01 | 0.58 | -42 % |

- The vet_gd halo share of kernel time on A100 fell 42 % -> 27 % (AG Car) and 49 % -> 31 % (BSG).
- He giant on A100: the cap leaves about 2.5 GB, so most of its shells keep the dense path (small gain, still bitwise).
- Caltech H200 timing: TASK-2026-10-10-caltech-vgdspeed.md on fork/bsg-files-1009.

## lin_tol print (89693b59, text only)

The summary line "final 7-point linear residual ... tol=<implicit_lin_tol>" suggested a missed tolerance; on BSG the
max was 1.1e-8 against tol 1e-10. That value is max|r|/max|b| of the true residual and is only a diagnostic:
- it enters the Picard test only with implicit_lres_test, which defaults to off for eddington, vet_sc and tau / vet_col;
- the inner BiCGStab stops on implicit_lin_cnorm (per cell |r|/(s E)), loosened by Eisenstat-Walker (ew_max).

The line now says this, followed by a second line that names the stop test. No numbers change.

## Gates on the combination (89693b59 against f39e26b1, same build options)

- CPU, login node, 4 ranks, reduced 16x16 lateral wedges: every bin / hst / rst file byte-identical.
  - AG Car A with production keys and rotate_every 3: 40 of 40 files.
  - BSG true-repro hr: 36 of 36.
  - He giant fresh: 28 of 28.
- GPU, Raven 4 A100, job 31038534, production inputs: hst byte-identical and bin data bitwise.
  - AG Car A hr, 10 cycles: 24 of 24.
  - BSG hr reduced, 5 cycles: 28 of 28.
  - He giant fresh, 3 cycles: 15 of 15.
  - AG Car with rotate_every 3: identical, and the new binary run twice is identical.
- 4-block / 4-rank vet_sc runs (merge_1009/vimp_cpu2 set: cyl central, cyl half-range, He box, He box half-range,
  He box half-range + implicit_vimp): rc 0, no FATAL (/viper/ptmp2/jinma/vgdspeed_1009/vscruns/RESULTS.txt).
- Regression subset (rad_m1, hydro, mhd, rad; GR skipped): /viper/ptmp2/jinma/merge_1009/logs/regress_vgdm.log
  (`_cpu` 62 passed, 3 skipped; `_mpicpu` 4 passed, 3 skipped, 1 failed: rad/test_rad_lwave2d_amr_mpicpu, the segfault already present on ba7212d6, NOTE-2026-10-09-viper-rt-integration-merge).
