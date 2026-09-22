# cs_wellbalanced_src under `<hydro>` (2026-09-22/23)

What was done: `<hydro>/cs_wellbalanced_src` is now read for pure-hydro cubed-sphere runs,
so hyd4 is exactly "prod4 minus B". Two gates: the MHD path is bitwise unchanged, and the
flag really acts in hydro. Then off-vs-on tests to answer whether the flag **helps pure
hydro**, and whether it can be made default-on. Nothing is committed or launched.

## 1. The patch (`tests_hyd4/cs_wbsrc_hydro.patch`, uncommitted in the working tree)

* `src/coordinates/coordinates.cpp`, in the `<hydro>` branch of the constructor (line 53):
  `cs_wellbalanced_src = pin->GetOrAddBoolean("hydro","cs_wellbalanced_src",false);`.
  The `<mhd>` branch is untouched and still takes precedence.
* `src/coordinates/coordinates.hpp`: comment only.
* **No kernel change was needed.** The hydro source path is `hydro_tasks.cpp:536`
  `SrcTermsGnomonicEquiangle`. It calls the same `SrcTermsGnomonicEquiangleImpl` as MHD
  (coordinates.cpp:1153, `is_mhd = false`), and that function branches to
  `SrcTermsCurvilinearWB(..., is_mhd=false, ...)` when `cs_wellbalanced_src` is set
  (coordinates.cpp:1494). `SrcTermsCurvilinearWB` already handles `mhd_ == false`, and
  `wb_geom` is built lazily. The only thing missing was the parameter read.

## 2. The binary: `/viper/u2/jinma/ATHENAK/bench/cs_hyd4_prod/athena`

md5 `fef759d9e02fb499572fa433363c6023`. Built from a `git archive HEAD` snapshot of
1d22a118 (`src CMakeLists.txt config.hpp.in`, kokkos symlinked) with the patch applied.
The recipe is prod4's: HIP GFX942_APU, MPI, Release, PROBLEM=deep_hot_jupiter_rt. cmake
and make both returned 0, with 18 warning lines, the same count as prod4. Changes in
`src/` since prod4 (395db5bc): `mesh/build_tree.cpp` (the per-rank restart fix, 886bd677),
`rad_m1/*` and `src/CMakeLists.txt` (rad_m1), and this patch. Details are in
`bench/cs_hyd4_prod/BUILD_COMMIT.txt`. Before building, an inode check created 3000 files
without error.

## 3. Gates (apudev, 1 node, 2 GPUs; job 11943423, `bench/cs_hyd4_prod/gate/`)

**(a) The MHD path is unchanged: PASS.** The prod4 input was run from scratch for 300
cycles, once with the prod4 binary (MO) and once with the new binary (MN). All 4 files
(bin 00000/00001, rst 00000/00001) have identical payload md5 after `<par_end>`
(`gate/payload_md5.sh`; bin 00001 = 0ced9a32..., rst 00001 = 18a8251e...).

**(b) The flag acts in hydro: PASS.**
* HF (new binary, hyd4 with `hydro/cs_wellbalanced_src=false`, 300 cycles): all 4
  payloads are bitwise identical to the earlier smoke H2 (prod4 binary, job 11943303). So
  the hydro path is otherwise unchanged.
* HW (flag on, 300 cycles): bin 00000 and rst 00000 are identical (the initial state), but
  bin 00001 and rst 00001 **differ** (8256cca1 / 8d9b38ba against db19920c / 6d38bafc).
  dt at cycle 200 is 28.33297 against 28.33253. The source term is active. It prints no
  startup line.
* HL (flag on, `-t 7:40`, the same wall limit as smoke HL): rc 0, no FATAL. It ran 6880
  cycles to t = 1.7035e5 s and wrote rst 00000-00002. Warnings are the same as the smoke
  run: one "explicit radiative source clipped in 2 cells by rt_de_max". No eos_fail. dt at
  cycle 3900 was 24.82 against 25.20, and 26.7 at 6800 against 26.0 (smoke). Counters per
  interval (dhj.log) are the same order as the smoke run: dfloor 7.0e6 -> 9.1e7, tfloor
  1.8e6 -> 6.1e6, efloor 1.8e5 -> 2.5e5. Mass drift was -0.139 % and energy -0.284 % by
  1.70e5 s (smoke: -0.139 % and -0.288 % by 1.75e5 s).

## 4. Does it help pure hydro? Verdict: **NO. It is neutral on dhj and hurts on the exact tests.**

The same binary was used throughout; the only difference between arms is
`hydro/cs_wellbalanced_src` off or on.

### (1) Rigid rotation, exact steady state (CPU, `bench/cs_hyd4_prod/wbtests/`)

The problem is cs_test iprob 3 (`rot_hyd.athinput`): plm, hllc, nghost 3, t = 1, and joint
refinement n = nx1 = nx2 = nx3 = 16/32/64. The region band is scaled 2/4/8 cells, so it
has a fixed physical width. Four flows: the axis runs through panel centres (z) or cube
vertices (v), with omega 0.2 (Mach ~0.3) or 1.0 (Mach ~1, strong rho v v). The table gives
on/off ratios at n = 64 (`wbtests/tables.txt`, logs in `wbtests/logs/`):

| flow | L1(v) INT / SEAM / VTX | L1(v_t) VTX | L1(p) INT / SEAM / VTX |
|---|---|---|---|
| z, 0.2 | 1.077 / 1.040 / **0.898** | **0.853** | 1.076 / 1.326 / 1.557 |
| v, 0.2 | 1.021 / 1.005 / 1.019 | 1.019 | 0.953 / 1.242 / 1.355 |
| z, 1.0 | **1.530** / 1.225 / 1.075 | 1.017 | 1.213 / **2.055** / 1.799 |
| v, 1.0 | **1.227** / 1.024 / 0.969 | 0.947 | 1.100 / 1.669 / 1.350 |

* In 3 of 4 flows the vertex velocity error drops by 3-15 %. Every other metric gets
  worse, most of all the pressure at seams and vertices (1.2x to 2.1x) and the
  interior velocity at Mach 1 (1.2x to 1.5x).
* The ratios grow with n in several columns (z/1.0 interior: 1.25 -> 1.41 -> 1.53). The
  flag does not improve the constant in front of the error at any resolution.
* The orders stay ~2 either way. With the flag on, the interior velocity order is
  slightly lower on z/1.0: 2.11/2.09 against 2.28/2.21.
* Linf(v) is unchanged to within 1-9 %. dt is identical (the cycle count is the same).
  CPU cost is +7-11 % (n=64: 464 -> 504 s, 697 -> 765 s).

### (2) Hydrostatic atmosphere at rest (cs_test iprob 13, `strat_hyd.athinput`, t = 1)

* **With the flag off, the tangential velocity is round-off:** L1(v_t) = 1e-16 in every
  region at every n. The cell-centre source is already exact for isotropic pressure, as
  the existing code comment says.
* **With the flag on, spurious tangential flow appears:** L1(v_t) = 9e-9 / 1.5e-10 /
  2.3e-12 at n = 16/32/64. The interior is the largest region and the vertex the
  smallest. It is tiny and converges fast, but it is nonzero where the default gives
  zero.
* The radial error (the hydrostatic imbalance) and L1(p) are bitwise identical to 4
  digits.

### (3) hyd4 dhj, equal simulated time (GPU, job 11943521, `bench/cs_hyd4_prod/gate2/`)

WF (off) and WT (on) both ran from scratch to tlim = 1.2e5 s (0.39 rotation), with bins
every 3e4 s. Both reached t = 1.2e5 at cycle 4925, and dt was the same within noise
(25.89/25.87 at cycle 1000, 23.40/23.83 at cycle 4500). The floor counters at the last
interval were dfloor 8.460e7/8.480e7, tfloor 5.55e6/5.54e6 and efloor 2.70e5/1.98e5, with
no eos_fail in either run. The warnings are identical. Final mass 3.59562e26 and total E
6.62959e38 are identical to 6 digits. Radial KE was 4.595e32 against 4.697e32 (+2 %);
tangential KE was 1.456e33/1.171e33 against 1.457e33/1.164e33.

Mass-weighted rms speed by region at t = 1.2e5 (`gate2/regions_ke.py`: 4-cell band on
32 cells per panel; TOP = the upper 16 of 128 radial cells):

| | INT vr / vt | SEAM vr / vt | VTX vr / vt |
|---|---|---|---|
| TOP, off | 2.416e4 / 1.789e5 | 2.477e4 / 2.015e5 | 2.610e4 / 2.261e5 |
| TOP, on | 2.407e4 / 1.788e5 | 2.479e4 / 2.014e5 | 2.620e4 / 2.261e5 |
| ALL, off | 1.537e3 / 3.412e3 | 1.508e3 / 3.168e3 | 1.657e3 / 3.005e3 |
| ALL, on | 1.551e3 / 3.406e3 | 1.551e3 / 3.151e3 | 1.581e3 / 3.024e3 |

At the top, where the vertex flow lives, off and on agree to <0.4 %. The full-column
vertex vr changes by -5 % at 1.2e5 s but +1 % at 6e4 s (bin 00002), which is the
difference between two chaotic runs, not a trend. **No vertex signature either way.**

GPU cost, same binary and node: 15.25 against 15.24 cycles/s over cycles 500-4500
(WF/WT). The 300-cycle pair HF/HW gave 15.35 against 14.99 (-2.4 %). With one repeat
each, the cost is 0-2.4 %.

### Earlier evidence (MHD, tests_cs_regions/README.md, tests_cs_vertex_cell/README.md)

* On MHD `strat` the flag made the vertex error about 8 % **worse** at every resolution.
  It was neutral on `rot` and `loop`.
* With wenoz, the long-time spurious velocity got worse (1.85e-2 against 1.13e-2).
* The MHD reason for the flag is different: the low-beta Maxwell-stress cancellation
  (the coordinates.cpp:1171 note). It has no hydro counterpart.

## 5. Should it be default-on for both `<hydro>` and `<mhd>`? **Not for hydro.**

* **Hydro:** it hurts on the exact steady tests and adds tangential flow to an atmosphere
  at rest that the default balances to round-off. It is neutral on dhj and costs 0-2.4 %
  GPU. Recommendation: keep the default off for `<hydro>`. hyd4 keeps it **on** only so
  that it is an exact "prod4 minus B" control. On dhj the effect is within run-to-run
  noise, so the control stays valid.
* **MHD:** these tests say nothing new. The case for it there is low-beta stability, and
  prod4 sets it explicitly anyway. A default-on for `<mhd>` would be a separate decision.
* **What a default flip would change bitwise.** Every cubed-sphere input that does not
  set the key:
  * `inputs/tests/`: cubed_sphere_{toroidal, rigidrot, resist_smr_narrow, uniform, mhd,
    ffdecay, blast, blast_open, resist_smr, mhd_smr, raddiff, mhd_conv, smr, resist,
    mhd_blast}, two_stream_sph_thick.
  * `inputs/hydro/`: red_giant_1d_dilution, red_giant_fofc, red_giant_cs, he4_presn_cs.
  * `tst/inputs/`: fofc_cs, fofc_mhd_cs, restart_cs_mhd. These feed
    `tst/test_suite/hydro/test_hydro_fofc_cs_cpu.py`,
    `mhd/test_mhd_fofc_cs_cpu.py` and `hydro/test_restart_bitwise_cpu.py`.
  * cubed_sphere_raddiff feeds `rad/test_rad_cs_raddiff_cpu.py` and
    `rad/test_rad_cs_implicit_ang_cpu.py`.

  I did not run these tests with a flipped default, so whether their thresholds still
  pass is unknown. The cs_regions_* and cubed_sphere_mhd_strat inputs set it to false
  explicitly and would not change.
* **Side effect of this patch:** any old hydro input that had `cs_wellbalanced_src = true`
  under `<hydro>` changes behaviour with the new binary. Before, the key was silently
  ignored (e.g. `cs_hyd_rs/base.athinput` in bench). In `inputs/`, only hyd4 has it.

## 6. The production directory is prepared but NOT submitted: `bench/cs_hyd4_prod/`

* `submit.sh` and `chain.sh` are prod4's, with the job name `cs_hyd4_prod` and the
  header changed. Same settings: apu, 2 GPUs on 1 node, `-t 23:50:00`, restart from the
  newest rst, STOP-file chain logic.
* `deep_hot_jupiter.athinput` is `inputs/production/deep_hot_jupiter_cs_hyd4.athinput`
  with only the banner changed. Its H3 now says `cs_wellbalanced_src = true` under
  `<hydro>`.
* Launch with `cd bench/cs_hyd4_prod && ./chain.sh . 4`. It starts from scratch; there
  is no rst/ yet.
