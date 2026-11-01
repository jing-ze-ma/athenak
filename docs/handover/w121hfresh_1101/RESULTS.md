# w121hfresh_1101: fresh-start (t = 0) HYDRO WASP-121b, `<hydro>/am_vector = true`, on viper MI300A

## Verdict
**GO on viper.** It is fast and healthy. 1 apudev node (2 MI300A) gives **3.7 min per rotation**, so 300 rot take
about 18.5 h on 1 node if dt holds (only the first 2.7 rot were measured). No FATAL, no eos_fail, no fofc and no vceil events in 6 runs.
**Cost of am_vector: +1.5 % per cycle** (same node: 0.0277 vs 0.0273 s/cycle). dt and memory are the same in both arms.
am_vector does act: the spin reservoir is nonzero in every active cell of every restart (see below).
Placement on viper or Caltech is the user's call.

## Code / build
- fork/fin51-r-fix = 182477553, detached worktree `/viper/ptmp2/jinma/w121hfresh_1101/wt` (with the kokkos submodule).
- Built with the `dhj_gpu72` recipe of `builds/build_inc_viper.sh` (gcc/16, rocm/7.2, hipcc 7.2.4, Kokkos_ENABLE_HIP,
  Kokkos_ARCH_AMD_GFX942_APU, Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF, MPI, Release, PROBLEM=deep_hot_jupiter_rt),
  full build 358 s. Script `build.sh`, logs `build_cmake.log` and `build_make.log`.
- Binary `bin/athena_dhj_gpu72_18247755_hfresh` (copy in `/viper/ptmp2/jinma/builds/bin/`), md5 38c3fa6cc0c27720b5643ee7c7b02bfe.
- fmad-1031: the CMake change appends `--fmad=false` only when Kokkos_ENABLE_CUDA is on and CMAKE_CXX_FLAGS contains
  -ffp-contract=off. This is a HIP build with neither, so it did nothing here. The compile lines carry `-O3 --offload-arch=gfx942`
  and no -ffp-contract or fmad flag, which leaves the hipcc default contraction as in every earlier viper dhj build.

## Deck: `deck/hfresh.athinput` (ON) and `deck/hfresh_off.athinput` (OFF)
Made from `/viper/ptmp2/jinma/w121_s10fresh_1027/pkg/s10fresh.athinput`. Run dirs symlink `ck` and `ic` to that package.
Changes vs the MHD deck (`diff` shows nothing else):
- `<mhd>` -> `<hydro>`; `rsolver lhlld -> lhllc` (lhllc_cb_sensor left at its default; the dhj default
  `<hydro>/lhllc_x1_phi_min = 1` is confirmed in the log).
- dropped MHD-only keys: `lhlld_x1_phi_min`, `ohmic_resistivity`, `max_eta`, `use_rkg_sts`, `valfven_max_cgs`;
  and `problem/bbot_gauss` (the B seed).
- added `<hydro>/am_vector = true`.
- added `problem/am_transport = false`. **Required:** with am_vector on, the pgen's default problem/am_transport
  (true in a rotating run) calls `AmModule::YieldToPgen`, which stops with a FATAL because the correction would
  be applied twice. Side effect: the pgen still prints "### WARNING: problem/am_transport = false ... L_abs drifts".
  With am_vector on, that warning is misleading.
- `<output3>/variable mhd_w_bcc -> hydro_w`; the comment line and header text.
- OFF arm: `am_vector = false` and no problem/am_transport line, so the pgen correction is on as in the MHD deck.
  The recipe log then says "<hydro>/am_transport = false (default yields: problem/am_transport is on)".
- The parse (`-n`, runs/parse_on/parse.log) confirms the production grid and keys: C32 (nx2 = nx3 = 32, MeshBlocks 112x16x16, 24 MBs),
  nx1 112, x1max 2.1e10, rk2, cfl 0.3, ck_impl_every 4, ck_impl_xstep 8, ck_impl_once, ck_nquad 2, cs_seam_flux positive,
  sponge_top on, sponge_bottom off, etotgrav true.
- Recipe: no `<job>/recipe` is set, so this fresh run gets the new defaults. The startup block prints
  `<hydro>/am_transport = true (default on the cubed sphere ...)` and `<time>/dt_retry = true (default)`. Nothing else is
  defaulted (no rad_m1). am_vector is named in the deck, so it gets no recipe line.

## Runs (apudev, 1 node x 2 GPU, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1, script `job.sh`)
| job | runs | node | result |
|---|---|---|---|
| 12523640 | parse_on, smk_on (40 cyc), parse_off, smk_off (40 cyc) | vipa1326 | rc 0, no FATAL (smoke) |
| 12523832 | on_a 10 min | vipa1326 | rc 0 |
| 12523833 | off_a 10 min | vipa1001 | rc 0 |
| 12523834 | off_b 10 min | vipa1001 | rc 0 |
| 12523835 | on_b 10 min | vipa1326 | rc 0 |
| 12524231 | on_c 5 min then off_c 5 min, SAME node | vipa1001 | rc 0 |

Numbers from `ana.py`, saved in `ana_out.txt`. s/cycle excludes the start-up and the first 100 cycles. 1 rot = 110154 s.

| run | cycles | rot reached | s/cycle | mean dt [s] | cyc/rot | min/rot (1 node) |
|---|---|---|---|---|---|---|
| on_a | 21100 | 2.623 | 0.0279 | 13.69 | 8044 | 3.7 |
| on_b | 21200 | 2.636 | 0.0278 | 13.69 | 8044 | 3.7 |
| on_c | 10500 | 1.296 | 0.0277 | 13.59 | 8106 | 3.7 |
| off_a | 21600 | 2.682 | 0.0273 | 13.68 | 8054 | 3.7 |
| off_b | 21600 | 2.682 | 0.0274 | 13.68 | 8054 | 3.7 |
| off_c | 10600 | 1.305 | 0.0273 | 13.55 | 8128 | 3.7 |

- dt settles at 13.79 s in all runs. A fresh start may still develop flows that lower dt later; these runs reached only 2.7 rot.
- **4 apu nodes (8 GPUs, 3 MBs per GPU): NOT measured, extrapolated.** Ideal scaling gives 0.93 min/rot. With only 3 MBs
  per GPU, communication will take a larger share, so expect clearly less than 4x. 1 node already gives about 18.5 h per 300 rot.
- Memory: GTT peak 20.5 GiB per GPU in both arms. Node MemTotal-MemAvailable rose by 43.7 GiB during the run
  (peak 59.7 GiB used on the 2-APU node). rocm-smi VRAM shows only 16 MB, as expected on an APU. Slurm MaxRSS (smoke) is 1.5 GB per task.
- Events (summed over the dhj.log rows of the 10-min runs; ON / OFF):
  eos_dfloor 4.87e8 / 4.98e8, eos_efloor 1.10e9 / 1.13e9, eos_tfloor 3.65e8 / 3.77e8; eos_vceil, eos_fail, c2p_it, fofc = 0.
  The two arms agree within 2 %. For comparison, the MHD fresh smoke 12405174 had 3.9e8 dfloor in only 600 cycles.
  The floor energy column Efloor reaches about 1e23 in both arms.
- Angular momentum (dhj.user.hst, first row / last row at about 2.6 rot), ON(on_a) vs OFF(off_a):
  Lrel 9.247e37 -> -1.009e37 vs 9.274e37 -> -1.015e37; Iw 8.0513e46 -> 8.0508e46 vs 8.0507e46;
  Tcor -5.70e32 vs -5.53e32; Tcc -1.614e34 vs -1.613e34; Tbnd -1.06e25 vs -1.04e25; Tcons -1.49e25 vs -1.56e25; Twall 0.
  |Lrel| / (Iw Omega) is about 2e-6, and the two arms agree within 0.5 %. The am_hst columns do not report the spin reservoir.

## Does am_vector act? Yes.
- The spin reservoir `am_spin` is written to the restart (restart.cpp, kAmSpinRstMagic) only when am_vector is on.
  ON restarts are 8.909 MB larger than OFF restarts, which is 24 MBs x 116x20x20 x 8 B.
- `spincheck.py` reads the last slab of each MeshBlock record (`spin_out.txt`). At rot 1, 2 and 3 it is nonzero in all
  688128 active cells and zero in the ghosts. max|s| grows 9.2e34 -> 1.8e35 -> 2.4e35, and on_a and on_b are bitwise
  equal at rot 1 and 2.
- The solutions also differ ON vs OFF (Lrel, Tcor, floors), so the key is not silently ignored.

## Why it is ~5-7x cheaper per rotation than the MHD runs
The grid and the ck cadence are the same (C32 x nx1 112, ck_impl_every 4, xstep 8, once, nquad 2), so neither explains the gap.
The difference is hydro vs MHD:
- Per cycle: MHD smoke 12405174 (same grid, 6fac5c7f) took 0.059 s/cycle (cycles 400-600) vs 0.0275 here, about 2.1x.
  MHD adds HLLD, face fields and CT, eos ohmic resistivity, and the RKG2 STS sub-cycle.
- dt: MHD 5.8 s (Alfven limited; the valfven_max_cgs 2e6 floor applies at the top) vs 13.8 s hydro, about 2.4x.
- Together that is about 5x per rotation against the MHD fresh smoke. The rot-320 MHD runs were not measured here.
  Their ~7x also includes whatever dt they have at a developed state.

## Files
deck/, job.sh, ana.py, spincheck.py, ana_out.txt, spin_out.txt, runs/<arm>/{run.log,dhj.log,dhj.user.hst,mem.log,rst}, logs/.
