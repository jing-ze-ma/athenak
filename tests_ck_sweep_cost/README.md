# tests_ck_sweep_cost — the cost of the correlated-k sweep, and what was recovered

Follow-on to `bench/prof_0922/README.md`, which measured the correlated-k sweep kernel
(`picket_fence_two_stream_RT_pass`, the `rt_chain_ck` chain + its pre-sweep) at **70.6 %**
of GPU time on the cs deep-hot-Jupiter production configuration, RT total 82 %, against
37 % on 09-12.  This directory asks where the extra went and takes back what could be
taken back **without changing a result**.

Everything below runs the same thing `bench/prof_0922/plain` ran: restart
`bench/cs_mhd_prod3/rst/dhj.00567.rst` (rot 283, `ncycle = 4429831`), the production
input with `nlim = 4430131` (300 cycles), 2 ranks / 2 GPUs on `apudev`, 24 MeshBlocks.
The reference number is `cpu time used = 23.04 s` for those 300 cycles, i.e.
**13.02 cycles/s**.

## 1. The redundancy, counted from the source

`rt_chain_ck` runs **four** full column passes per thread once `ck_spherical` is on:

| pass | what it is | added by |
| --- | --- | --- |
| P1 | plane-parallel DOWN probe `d0`, parked in `I_down` | `3c5d4846` (`ck_spherical`) |
| P2 | UP probe carrying the face mixing against `d0`, `u_below` parked in `Cmx` | `3c5d4846` |
| down | the real down-sweep, forming `c` from `Cmx` and overwriting it | — |
| up | the real up-sweep, reusing `c` | — |

and, since `d8f00d49` (the centre-to-centre layers), **each pass crosses two half layers
per cell** instead of one whole-cell layer.

The decisive observation — already written into the source, in the `ck_impl_frozen_op`
note above `step` — is that **every half layer any of the four passes ever crosses is the
half of ONE cell taken at ONE chain's `mu`**, so a single `(kappa; e0, c_in, c_out)` per
`(cell, chain)` describes the entire operator.  The implicit path exploits that through
the global `ckc0_g/ckci_g/ckco_g/ckkro_g` triple.  **The explicit path — which is what
production runs — did not.**  It rebuilt the operator from the k-table in every pass, and
twice within each pass (once as the lower half of cell `i` in iteration `i`, once as the
upper half of cell `i` in iteration `i-1`):

| per (cell, chain) | plane-parallel | `ck_spherical` | after `ck_sweep_cache` |
| --- | --- | --- | --- |
| `ck_kappa` look-ups | 2 | **4** | **1** |
| `expm1` (double) | 4 | **8** | **1** |
| double divides in the layer coefficients | 12 | **24** | **3** |
| `BFace` calls | 4 | 8 | 8 (unchanged) |

With 88 chains per column (11 bands x 8 g-points, `ck_nquad = 1`) and RT_NB = 4 chains
per thread, the spherical form therefore doubled the transcendental and k-table work of
the kernel on top of doubling the number of recurrences.

## 2. What was changed — `problem/ck_sweep_cache`

`src/utils/two_stream_rt.hpp`, the `rt_chain_ck` kernel only (plus the one `GetOrAddInteger`
in `src/pgen/deep_hot_jupiter_rt.cpp`).  The kernel gains a compile-time tag `CCH` next to
the existing `SPH`/`BSP`/tier tags, so every setting is its own instantiation and the ones
that are off cost nothing.  Three settings, because the two kinds of reuse buy and cost
different things:

| mode | what it reuses | per (cell, chain), `ck_spherical` | extra private storage |
| --- | --- | --- | --- |
| **0** | nothing — the code as it was | 4 look-ups, 8 `expm1`, 24 divides | 0 |
| **1** | **carry**: within ONE pass the second crossing of a half layer reads the triple the first left in `cry0/cryi/cryo`, `NC` registers, exactly like `kfar` | 4 look-ups, **4** `expm1`, **12** divides | **0** |
| **2** | **column**: the first pass over a cell (P1 under `ck_spherical`, the down-sweep otherwise) also parks kappa and the triple in `Kpc`/`Cc0`/`Cci`/`Cco`, `NC x NN` each, and the other three passes load them | **1** look-up, **1** `expm1`, **3** divides | `4 x NC x NN` = **17408 B/thread** at the production tier (`NN = 136`, `NC = RT_NB = 4`, `RtF = Real = double`), on top of the 8704 B `I_down` + `Cmx` already use |

Mechanically: `step` takes two flags, `fill` (the one call per PASS that computes fresh —
mode 1) and `first` (the one call in the whole COLUMN that first visits this half layer —
mode 2); a new `kapof`, which `krof` now calls, takes `first` alone.  `Kpc` caches kappa
and **not** kappa rho, so the beam's `kap` and every `kappa rho` are still formed by the
expression that formed them.

Nothing else moved: no reassociation, no hoisting across a rounding boundary, no change to
`ck_nquad`, the bands, the beam, the geometry or the Jacobian assembly, and the
`ck_impl_frozen_op` global store/load path is untouched and still takes priority.  The
numbers stored are the numbers the recomputation would have produced, from the same
operands in the same order, so **the kernel is bitwise identical at every setting** —
which is the A/B that gates it (section 4).

Default: **2 when `ck_spherical` is on, 1 otherwise** — mode 2 pays for its private memory
only when there are four passes to amortise it over (section 3.2).

## 3. GPU timing

### 3.1 How much `ck_spherical` costs today (step 1 of the brief)

`a1_plane`, job **11933405**: the **same binary** `bench/prof_0922/plain/athena` and the
same 300 cycles as the `plain` arm of the profile, with only `problem/ck_spherical`
flipped to false (the repaired plane-parallel path):

| arm | `ck_spherical` | `cpu time used` | cycles/s |
| --- | --- | --- | --- |
| `bench/prof_0922/plain` (job 11932897) | true | 23.04 s | 13.02 |
| `a1_plane` (job 11933405) | false | **17.47 s** | **17.17** |

**The spherical form costs 1.32x on the whole loop.**  Charged to the sweep alone, which
the profile puts at 13.05 s of that 23.04 s loop, the 5.57 s difference makes the
spherical sweep **1.74x** the plane-parallel one — far more than the 1.22x per-cycle
figure carried over from `tests_gpu_ab_0922`, and consistent with the count in section 1:
the probe passes double the recurrences AND doubled the transcendental/table work.

### 3.2 The cache (step 3 of the brief)

`a2_gpu`, job **11933619**: five 300-cycle runs in one job, same restart and input.
`athena.base` is the pre-edit tree, `athena.opt` the edited one; both built from
`bench/cksw_0922/src_snapshot` with the `bench/prof_0922` recipe (hipcc 6.3.4,
`Kokkos_ENABLE_HIP`, `Kokkos_ARCH_AMD_GFX942_APU`, MPI, Release).

| run | binary | `ck_spherical` | `ck_sweep_cache` | `cpu time used` | cycles/s | vs its control |
| --- | --- | --- | --- | --- | --- | --- |
| r0 | base | true | (absent) | 29.886 s | 10.04 | — |
| r1 | opt | true | off | 29.782 s | 10.07 | **+0.3 %** (the template tag is free) |
| r2 | opt | true | **on** | **27.836 s** | **10.78** | **+7.0 %** |
| r3 | opt | false | off | 21.859 s | 13.72 | — |
| r4 | opt | false | on | 22.152 s | 13.54 | **-1.3 %** |

Two things to read off this table.

**r0 vs r1 is the control and it is clean**: adding the extra template tag and the two
flags to `step`/`kapof`, with the cache off, changes nothing (0.3 %, inside run-to-run
scatter).  So r1 -> r2 is the cache and nothing else.

**r2: +7.0 % of the WHOLE loop** on the production (spherical) configuration.  Charged to
the sweep, 1.95 s off a ~17 s sweep is about **12 % of the kernel**.  The spherical
overhead itself falls from 1.363x (r1/r3) to **1.257x** (r2/r3).

**r4: the cache LOSES 1.3 % when `ck_spherical` is off**, which is exactly the trade
section 5 predicts: with only two passes there are 3 of 4 `expm1` and 1 of 2 look-ups to
save, and that no longer pays for 17.4 kB/thread more scratch.  That is what motivated the
third setting, measured next.

### 3.3 The three settings — `a3_gpu`, job 11933921

Same job structure, binary `athena.opt2` (the shipped three-mode kernel), in-job control
in every pair:

| run | `ck_spherical` | mode | `cpu time used` | cycles/s | vs mode 0 at the same `ck_spherical` |
| --- | --- | --- | --- | --- | --- |
| r0 | true | 0 | 29.821 s | 10.06 | — |
| r1 | true | 1 (carry) | 29.371 s | 10.21 | **+1.5 %** |
| r2 | true | **2 (column)** | **28.706 s** | **10.45** | **+3.9 %** |
| r4 | false | 0 | 22.020 s | 13.62 | — |
| r3 | false | **1 (carry)** | **20.700 s** | **14.49** | **+6.4 %** |

So the two modes split cleanly along the number of passes, and the default follows the
measurement:

| `ck_spherical` | best mode | gain | why |
| --- | --- | --- | --- |
| **true** (production) | **2** | **+3.9 % here, +7.0 % in 3.2** | four passes amortise the 17.4 kB |
| **false** | **1** | **+6.4 %** | two passes; the carry is free, the column store is not (-1.3 %, a2 r4) |

`problem/ck_sweep_cache` therefore defaults to **2 under `ck_spherical`, 1 otherwise**,
which is what `deep_hot_jupiter_rt.cpp` sets.

One honest caveat on the size of the spherical win.  Mode 2 measured **+7.0 %** in
`a2_gpu` (binary `athena.opt`, two settings, 32 instantiations) and **+3.9 %** in `a3_gpu`
(binary `athena.opt2`, three settings, 48 instantiations).  The mode-0 controls agree
between the two jobs to 0.13 % (29.782 vs 29.821 s), so this is not run-to-run scatter:
the three-mode build's mode-2 kernel really is ~3 % slower than the two-mode build's,
presumably register allocation moving under the extra (dead, in mode 2) carry registers
and the larger instantiation set.  **Open**: a build that drops mode 0 — which exists only
to run this A/B — may recover it.  The conservative claim is **+4 %** and the measured
best is +7 %.

Note on the absolute level: r0/r1 at 29.9 s are 1.30x slower than the 23.04 s the
`bench/prof_0922/plain` arm recorded for the identical 300 cycles.  That gap is NOT this
change — it is present in `athena.base`, which is `rt-integration` HEAD 446b7b30 plus the
uncommitted Jacobian-reuse work, against the profile's a544a761 — and it is why every
comparison above is made against a control inside the same job.  **It is worth a
bisection of its own** and is not explained here.

## 4. Bitwise status — verified, both settings of `ck_spherical`

`cpu_ab/run.sh`: `inputs/tests/dhj_ck_spherical.athinput` (nx1 = 128, nx2 = nx3 = 8,
CPU build `build_cksw`), 20 cycles, history every cycle and a `hydro_w` binary dump every
cycle, run four times — `ck_spherical` x `ck_sweep_cache`:

| `ck_spherical` | pair | `.hst` | the 22 binary dumps |
| --- | --- | --- | --- |
| false | mode 0 vs 1 | **bitwise identical** | **0 of 22 differ** |
| false | mode 0 vs 2 | **bitwise identical** | **0 of 22 differ** |
| true  | mode 0 vs 1 | **bitwise identical** | **0 of 22 differ** |
| true  | mode 0 vs 2 | **bitwise identical** | **0 of 22 differ** |

("differ" above means outside the header.  `cmp -l` on every dump reports exactly **5
differing bytes**, at the same offset in every file, and they are the five characters of
`ck_sweep_cache = 0/1/2` inside the parameter block the dump header echoes.  No field byte
differs anywhere.)

`gate.sh`: the `tests_ck_sph` README section 0.2 thermal gate, re-run at `t = 0` on the
spherical form at all three settings.  The three one-shot column dumps (`ck_dump_file`,
`mu0 = 0.922`) are **byte-identical to each other**, and the residuals on that column are

```
THERMAL  max per-cell |V dep - dPhi|/max|Phi| = 2.808e-10      (gate: 1e-10 round-off)
THERMAL  column budget                        = 3.871e-09      (gate: 1e-9  round-off)
BEAM     absorbed power sum V Q = 1.31810e+29, /(A_top mu0 F*) = 0.37814,
         /(A(tau=1) mu0 F*) = 1.03761, top-face beam flux ratio = 0.99133
```

— the beam line reproducing `tests_ck_sph/README.md` section 0.2 digit for digit.

No reassociation was needed anywhere, so there is no 1-ulp caveat to state.

## 5. Occupancy / launch geometry

From `bench/prof_0922/prof/prof_kk/rank_0_kernel_trace.csv` (the dominant kernel,
`Kernel_Id = 100`, 600 dispatches, 13052.6 ms):

```
workgroup      256 threads          grid  264 workgroups x 256 = 67584 threads / GPU
                                          = 12 MeshBlocks x 22 chain blocks x 16 x 16
private segment  9792 B/thread      group segment 0
```

67584 threads is **1056 waves per GPU**.  On a GFX942 APU (228 CUs, 4 SIMDs each) that is
**~1.16 waves per SIMD** — the kernel has essentially no other wave to hide a latency
behind, which is the regime the 09-12 note (`bench/prof/rt_new/NOTES.md`) already
identified and which makes per-thread serial work, not VALU throughput, the thing that
costs.  It is also why this work targeted the count of `expm1`/divides/table look-ups on
the critical path rather than the flop count.

Register usage, read out of the code object metadata of both HIP binaries
(`llvm-objcopy --dump-section=.hip_fatbin`, `clang-offload-bundler --unbundle`,
`llvm-readelf --notes`): every `picket_fence_two_stream_RT_pass` instantiation reports
**vgpr 128, sgpr 108, max_flat_workgroup_size 1024** in BOTH binaries — the spherical form
did not raise register pressure, it is pinned at the 128-VGPR cap either way, so occupancy
is not register-limited and was never going to be improved by shrinking the kernel.

What the cache does cost is private memory: the production tier goes from **9840 B** to
**~27.2 kB per thread** (the `base` binary's whole set of private-segment sizes reappears
unchanged in `opt` as the `CCH = false` instantiations, with the larger `CCH = true` ones
added alongside).  That is the trade the timing below settles.

## 6. What was looked at and NOT taken

- **The pre-sweep is not the problem.**  `rt_pre_geom`, `rt_pre_tp`, `rt_pre_cut`,
  `rt_pre_opac` and `rt_pre` are per-column, chain-independent kernels; the k-table is
  read there once per cell per band, not once per chain, and the profile's three
  `picket_fence_two_stream_RT_pass` specializations are dominated by the one carrying the
  chain (13052.6 ms of 13352.0 ms).  The per-stage rebuild that IS visible in the profile
  is `Conduction::BuildRadWeights` (601 calls, 594.3 ms, 3.2 % of GPU time), and it reads
  the current `w0`, so it is genuinely stage-dependent and cannot simply be hoisted.
- **`BFace`/`s_l`/`s_u`/`s_f` are also built four times per face per chain** (8 `BFace`
  calls and 4 divides per cell pair).  `BFace` early-returns `b_far` on any smoothly
  varying column, so the arithmetic left is small; caching `s_f` would need three more
  `NC x NN` columns, and section 3.2 shows that private memory is the binding cost.  Not
  taken.
- **`BTF` is evaluated twice per face** (once in P2, once in the down-sweep) with the same
  arguments.  It is outside the chain loop, i.e. 1/4 of a chain's weight, and it is two
  divides.  Not taken.
- **`Bb_g` is loaded ~8 times per (cell, chain)** instead of once.  The layout is
  `(m,b,i,k,j)` with `j` fastest, so those are coalesced L1-resident loads.  Not taken.
- **Occupancy is the elephant and it is not fixable here.**  The kernel runs 67584 threads
  = 1056 waves per GPU, ~1.16 waves per SIMD (section 5).  The one knob that would change
  that is `RT_NB` (4 chains per thread; `RT_NB = 1` would quadruple the thread count), and
  it is a compile-time constant that changes the ORDER in which chains accumulate into
  `Src_g`/`Fb_g`/`Em_g` — so it is NOT bitwise, it is a round-off-level change to the
  answer, and the source already records that `RT_NB = 2` was "faster on rt_chain, slower
  on the total".  Flagged, not taken: it needs its own A/B with a physics sign-off.
- **The 1.30x gap between `athena.base` (29.9 s) and `bench/prof_0922/plain` (23.0 s)** on
  the identical 300 cycles is unexplained and predates this work.  Worth a bisect.

## 7. Reproducing

```
# GPU (apudev, 2 ranks / 2 GPUs, 300 cycles from the rot-283 restart)
tests_ck_sweep_cost/a1_plane/submit.sh    # job 11933405, ck_spherical off, prof binary
tests_ck_sweep_cost/a2_gpu/submit.sh      # job 11933619, base vs opt, cache off/on
tests_ck_sweep_cost/a3_gpu/submit.sh      # job 11933921, the three modes
# CPU correctness (build_cksw: cmake -B build_cksw -D PROBLEM=deep_hot_jupiter_rt)
tests_ck_sweep_cost/cpu_ab/run.sh         # 6 x 20 cycles, bitwise A/B
tests_ck_sweep_cost/gate.sh               # the tests_ck_sph thermal gate, 3 modes
```

HIP binaries and the source snapshot they were built from:
`bench/cksw_0922/{src_snapshot,build.sh,athena.base,athena.opt,athena.opt2}`.
`athena.base` is the tree before this change, `athena.opt` the two-setting version,
`athena.opt2` the shipped three-setting one.
