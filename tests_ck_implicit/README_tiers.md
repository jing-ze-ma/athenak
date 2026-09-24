# ck-tiers: no compile-time column tiers in the ck / grey RT kernels

Branch `ck-tiers`, from rt-integration d2572b4f. Worktree: /viper/ptmp2/jinma/wt_cktiers.
Work dir: /viper/ptmp2/jinma/cktiers_0924. It holds the binaries, scripts, logs and scratch
tables. The build dirs are deleted.

## 1. What changed (src/utils/two_stream_rt.hpp only)

ck-scratch (README_scratch.md) moved the chain kernel's column arrays into the device
buffer `rt_ckscr_ptr`. This branch does the same for the three kernels that still had
column tiers (72 / 136 / 264 / 520):

| kernel | array that moved | rows per thread |
|---|---|---|
| `rt_chain_ck_lin` (ck_impl_lin_thr = 4) | `Sc[RT_NB][NN]` | RT_NB |
| `rt_chain_ck_lin1` (ck_impl_lin_thr = 1, production T4 / c2) | `Sc[NN]` | 1 |
| `rt_chain_grey` (rt_grey: red_giant, box_convection) | `I_down[2][NN]`, `I_upb[2][NN]` | 4 |

- **Layout.** The layout is the same wave-tiled one as the chain kernel. A new helper,
  `CkScrRows`, gives row r, element i of lane l at `tile + (r*n1 + i)*64 + l`. Each
  64-thread wave owns one contiguous tile, and the thread number is the par_for
  flattening. `Sc[cc][i]` and `I_down[q][i]` keep their syntax, so the kernel bodies are
  textually unchanged.
- **Buffer.** The kernels share `rt_ckscr_ptr` (`CkScrEnsure`).
  - ck_lin and ck_lin1 each need nch·n1 Reals per column. That is 1/7 of the chain
    kernel's groups in semi / T4, so at the call's size they never grow the buffer.
  - The grey chain needs 4·n1 Reals per column.
- **Templates.** The kernels no longer take an NN tag. `launch_ck_lin_tier` is now a
  plain lin_thr switch, and each kernel is instantiated once.
- **Cap removed.** The `n1 <= 520` fatal on the correlated-k path is gone.
- **Registers.** The small per-chain state (ss, sfc, suc and the other [NC] arrays)
  stays in registers.
- `jlin` (`rt_chain_ck_jlin`) never had a column array, since it writes lpf / lps / lpj
  Views. It is unchanged.

## 2. Scratch per lane (bytes, `.private_segment_fixed_size`)

Measured with `tools/scratch_all.py` / `tools/rtscr.py`. Tables: `scrall_{base,new}.txt`.

| kernel | tier 72 | 136 | 264 | 520 | new (one instantiation, any n1) |
|---|---|---|---|---|---|
| rt_chain_ck_lin (VGPR 123-125 -> 127) | 2320 | 4368 | 8464 | 16656 | **0** |
| rt_chain_ck_lin1 (VGPR 89 -> 89) | 592 | 1104 | 2128 | 4176 | **0** |
| rt_chain_grey (VGPR 128 -> 128) | 2688 | 4736 | 8832 | 17024 | **352** |
| rt_chain_ck (ck-scratch, unchanged) | - | - | - | - | 964-1984 (identical) |
| rt_chain_ck_jlin | no column array | | | | unchanged |

## 3. Gates

### CPU (serial, login node, base d2572b4f vs new)

Scripts: `cpu/gates.sh`, `cpu/rgate.sh`, `cpu/mhdgate.sh`. Outputs: `cpu/gates.out`,
`cpu/gates3.out`, `cpu/rgate/gate.out`, `cpu/mhdgate/gate.out`.

| gate | arms | result |
|---|---|---|
| well-posed A2 (100 calls): rec.txt, hst, bin, rst | semi, t4, c2 | **BITWISE** |
| full hydro, 6 cycles from the wp IC (README_fast sect. 2) | semi, t4, c2, t4 with lin_thr = 4 (ck_lin; jac_lin off) | **BITWISE** |
| restart R: 8 straight vs 4+4 and 3+5; 6 straight vs 3+3 | semi, t4, c2, c2 + every 4 | **all BITWISE** (the c2 3+3 lands on the xstep rebuild) |
| restart D: base vs new, 6 cycles | semi, t4, c2, c2 + every 4 | **BITWISE** |
| restart O: a base-written restart read by new | same | **BITWISE** |
| MHD, 12 ranks, c2 + every 4, 6 straight vs 3+3 (README_restart) | prod4 dhj.00148.rst, rotation 74.00 (t = 2.25700e7 s, P = 3.0501e5 s), read in place | **BITWISE** (rst, hst, bin) |
| MHD, the same run: base vs new, 6 cycles straight | same | **BITWISE** |
| grey chain: red_giant_1d_dilution, 4 cycles (red_giant binaries) | rt_grey | **BITWISE** (hst, bin, rst) |

### GPU gate at nx1 = 128 (job 11966892)

Setup:
- apudev, 2 ranks, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1.
- Restart: bench/cs_hyd4_prod/rst/dhj.00183.rst, read in place. t = 2.79075e7 s,
  **rotation 91.50**, cycle 1444516.
- Input: prod_fast.athinput. tlim = T0 + 3000 s, which is 140 cycles at the current dt.
- Arms run interleaved fwd / rev / fwd.

| arm | base r1 / r2 / r3 | new r1 / r2 / r3 | median base -> new | base vs new outputs |
|---|---|---|---|---|
| semi | 50.51 / 50.35 / 50.42 | 49.49 / 50.55 / 50.57 | 50.42 -> 50.55 (+0.3 %) | hst + rst **BITWISE** |
| c2 | 44.27 / 44.36 / 44.12 | 44.72 / 44.20 / 44.36 | 44.27 -> 44.36 (+0.2 %) | hst + rst **BITWISE** |

Both arms are within the 3 % limit. semi is bitwise because its kernel (the chain kernel)
is unchanged. c2 is bitwise here although ck_lin1's Sc moved: in this case the move did
not change FP contraction. That is unlike the chain kernel under ck-scratch, where semi
went round-off.

## 4. Large nx1 (cs dhj hydro from the analytic IC, 1 node, 2 GPUs)

Setup:
- Arm: production RT (T4 + c2 + ck_impl_every = 4, verbose), 30 cycles.
- Input: `gpu/big_hyd.athinput`. This is prod_fast plus the `scratch_level` and
  `ck_impl_every` keys, so the command line can set them.
- 16x16 blocks. nx2 = nx3 = 64 gives 48 blocks/GPU for 256 and 320. nx2 = nx3 = 48 gives
  27 blocks/GPU for 512 and 640.
- hydro/scratch_level = 1 at nx1 >= 320, in both arms. The hflux team scratch is the
  reason.
- ms/cycle is measured over cycles 5-30.
- Memory is the peak rocm-smi "VRAM Total Used" of GPU0, sampled every 2 s.

| nx1 | blocks/GPU | base ms/cyc | new ms/cyc | NOT-CONV | peak mem GB base / new | jobs |
|---|---|---|---|---|---|---|
| 256 | 48 | 175.9, 226.1, 210.9, 187.7 (median 199.3) | 197.1, 224.7, 231.9, 186.3 (median 210.9) | 0 / 0 | 85.1 / 85.1 | 11967655, 11967725 |
| 320 | 48 | 228.0, 213.0, 213.3, 256.8 (median 220.7) | 221.8, 215.9, 209.4, 224.2 (median 218.9) | 0 / 0 | 103.0 / 100.9 | same |
| 512 | 27 | 253.3 | 246.0 | 0 / 0 | 95.2 / 93.2 | 11967718 |
| 640 | 27 | **cannot run** (fatal: n1 = 644 > tier 520) | 405.2 | - / 0 | - / 111.5 | 11967655 |

Wall ms/cycle on these 30-cycle runs scatters by ±15 % from run to run in both binaries.
The kernel-time profile below is the reliable comparison. Every run printed 8
ck_implicit lines and 0 NOT-CONVERGED.

Per-kernel GPU time (rocprofv3 kernel traces, 20 cycles, window of 19 cycles, rank 0;
jobs 11967656 and 11967719; `gpu/<job>/prof_*/prof_summary.txt`). Values in ms/cycle:

| nx1 | kernel | base | new |
|---|---|---|---|
| 256 | rt_chain_ck | 36.46 | 36.11 |
| 256 | rt_chain_ck_lin1 | 10.09 | **8.24 (-18 %)** |
| 256 | rt_chain_ck_jlin | 5.60 | 5.61 |
| 256 | all kernels | 124.05 | 122.03 |
| 512 | rt_chain_ck | 82.21 | 81.16 |
| 512 | rt_chain_ck_lin1 | 8.34 | **6.94 (-17 %)** |
| 512 | rt_chain_ck_jlin | 4.39 | 4.38 |
| 512 | all kernels | 175.50 | 172.90 |
| 640 | rt_chain_ck | - | 109.46 (split into 2 launches by the 16 GB rt_ckscr_gb cap) |
| 640 | rt_chain_ck_lin1 | - | 8.28 |
| 640 | rt_chain_ck_jlin | - | 5.62 |
| 640 | all kernels | - | 219.96 |

Other limits found along the way:
- MHD cannot run nx1 >= 320 on either binary. `mhd::CalculateFluxes` uses level-0 team
  scratch sized by ncells1, and it fails with "could not find a valid team size" (job
  11966893). mhd_fluxes.cpp hard-codes `scr_level = 0`, where hydro has the
  `hydro/scratch_level` key. So the big-nx1 runs are hydro.
- 512 at 48 blocks/GPU is out of memory in both binaries: "failed to allocate 6.5 GiB
  (ck_lpf)" at 118 GB used.

## 5. What still has a cap

- **Picket-fence split path** (`rt_chain`: rt_ck = false, rt_grey = false,
  rt_split = true). It uses `I_ir_down_c[NC][RT_NNC]`, with RT_NNC = 72 at compile time.
  It fatals when n1 > RT_NNC. It is a test path and was not converted.
- **Monolithic picket-fence path** (`launch_grey_rt`, `2stream_rt` reduce, rt_split =
  false). Its tiers run 72 to 1032, and it fatals above 1032. The 1032 tier is 66.9 kB/lane.
- **Legacy RT functions in src/pgen/deep_hot_jupiter_rt.cpp** (`constexpr int NN = 270`,
  two old `2stream_rt` kernels). They are not reached by the ck / grey paths.
- **rt_chain_ck.** Its buffer is capped at `rt_ckscr_gb` = 16 GB per launch. That cap
  splits the launch; it does not cap n1.
- The production paths have no n1 cap: the ck chain, ck_lin, ck_lin1 and jlin, and the
  grey chain. The limits there are memory (lP / lpf) and, for MHD, the flux team scratch.

## 6. Files

In /viper/ptmp2/jinma/cktiers_0924:
- Binaries, each with a `.commit` file:
  - `athena.{cpu,gpu}.base` (d2572b4f)
  - `athena.{cpu,gpu}.new` (6732bde9; same kernel code as this commit)
  - `athena.cpu.{base,new}.red_giant`
- `build.sh`
- `cpu/{gates,rgate,mhdgate}.sh`
- `gpu/{common.sh,a128.sub,big.sub,prof.sub}`
- `tools/{scratch_all,rtscr,big,ana_prof}.py`
