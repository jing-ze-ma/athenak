# tests_ck_kernel_tags — compile-time tags on the correlated-k chain kernel

Fixes what `bench/bisect_cost/README.md` found: `problem/ck_impl_frozen_op` (commit
`90c33ce5`) was a RUNTIME bool inside `rt_chain_ck`, so both sides of ~8 branches per
(cell, chain) and five stored-operator Views were compiled into the production sweep with
the switch off.  The standing rule is the one `ck_sweep_cache` already follows: **a
default-off switch inside the sweep kernel must be a compile-time tag.**

## 1. The change

`launch_ck_chain` takes a sixth tag, `fop_tag` -> `constexpr bool FOP`, next to
`NN / SPH / BSP / CCH / FRM`.  Inside the kernel the three selectors become

```
const bool ckfst = FOP && ckfst_;   // this pass STORES the operator
const bool ckfus = FOP && ckfus_;   // this pass RE-APPLIES it
const bool ckfcf = FOP && ckfcf_;   // ... including the coefficient triple
```

so at `FOP = 0` every frozen branch is dead code.  The tag is dispatched in
`launch_ck_cache`, i.e. it multiplies only the cached forms: the probe-free
`ck_sweep_form` recurrences are refused under `ck_implicit`, which is the only way the
frozen operator can be on.  Instantiations of the chain kernel: 64 -> 112 (16 -> 28 per
radial tier); the pgen translation unit costs +4 min on the HIP build (6:33 -> ~10:30 for
a full build, 4:00 for the TU alone) and 1:36 on the CPU build.

Nothing else changed: no flag semantics, no defaults, no input files.

## 2. GATES — bitwise, 20 cycles, hst + bin dumps (`gate.sh`, `md5.txt`)

Reference `athena.cpu.head` = HEAD `0172f3bc` unmodified, new = this tree, both
`cmake -D PROBLEM=deep_hot_jupiter_rt`, Release, CPU.

| gate | input | settings | result |
| --- | --- | --- | --- |
| a | `dhj_ck_spherical.athinput` | as written (`ck_spherical = false`) | **BITWISE** (23 files) |
| b | `dhj_ck_spherical.athinput` | `ck_spherical = true` (-> cache 2, form tm) | **BITWISE** (23 files) |
| c | `dhj_ck_implicit.athinput` | `ck_implicit = true`, `rt_use_cons = true` | **BITWISE** (23 files) |
| d | `dhj_ck_implicit.athinput` | c + `ck_impl_frozen_op = true` | **BITWISE** (23 files) |
| e | `dhj_ck_implicit.athinput` | d + `ck_spherical = true` | **BITWISE** (23 files) |

23 files = 1 `.hst` + 22 `hydro_w` binary dumps per arm; the dumps were deleted after the
comparison, the per-arm md5-of-md5s is in `md5.txt`.  cpplint: the file carries one
pre-existing `whitespace/blank_line` note and no new one; `count.py` is flake8-clean.

## 3. Code object — instructions / scratch per instantiation

`./dis.sh <obj>` (= `bench/bisect_cost/dis5.sh` with a tag-aware counter; `count.py`
decodes the substitution-compressed tags out of the raw Itanium name, which is what the
old `top2.py` got wrong).  Full tables: `head.counts.txt`, `tag.counts.txt`.

Radial tier 136 (the production tier), matched instantiation by instantiation:

| SPH | BSP | CCH | FRM | HEAD instr / scratch | FOP=0 instr / scratch | instr |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0 | 2 | 1 (**production, tm**) | 15544 / 1099 | **15372 / 1005** | -1.1 % (scratch -8.6 %) |
| 1 | 0 | 2 | 2 (sd) | 15560 / 1245 | 13189 / 843 | -15.2 % |
| 1 | 0 | 2 | 0 (four-pass) | 24203 / 1638 | 24270 / 1815 | +0.3 % |
| 1 | 0 | 0 | 0 | 28389 / 1946 | 27974 / 1686 | -1.5 % |
| 0 | 0 | 1 | 0 (**plane-parallel default**) | 20323 / 1452 | 17427 / 1347 | -14.3 % |
| 0 | 0 | 0 | 0 | 21990 / 1474 | 18225 / 1287 | -17.1 % |
| 0 | 0 | 2 | 0 | 19328 / 1381 | 18766 / 1391 | -2.9 % |

The `FOP = 1` instantiations sit where HEAD's kernel sat (e.g. tier 136, SPH=1, BSP=0,
CCH=2, FRM=0: 23863), i.e. the frozen path is unchanged, it has simply been moved out of
the explicit kernel.

**A retraction of one of `bisect_cost`'s numbers.**  That note read the production
instantiation as the 24124 -> 28391 pair.  With the tags decoded properly, 28391 is the
`CCH = 0` spherical kernel and the production one today (`CCH = 2, FRM = 1`) is 15544.
Removing the frozen-op branches from it recovers 1.1 % of its instructions but 8.6 % of
its scratch traffic — and the kernel is scratch/latency-bound, which is why the timing
below moves by much more than the instruction count does.

## 4. GPU — job 11941488 (one apudev job, 4 arms, 2:12 wall)

`bench/prof_0922/plain` verbatim: restart `cs_mhd_prod3/rst/dhj.00567.rst` (rot 283),
production input, 2 ranks / 2 GPUs, `time/nlim = 4430131` = 300 cycles.  Binaries built
from `bench/cktag_0922/src_snapshot` with the `bench/cksw_0922/build.sh` recipe.

| arm | binary | setting | `cpu time used` | cycles/s | vs HEAD |
| --- | --- | --- | --- | --- | --- |
| r0 | HEAD `0172f3bc` | defaults (spherical, tm) | 21.68 s | 13.84 | — |
| r1 | **FOP tag** | defaults (spherical, tm) | **19.66 s** | **15.26** | **-9.3 %** |
| r2 | HEAD | `ck_spherical = false` | 20.95 s | 14.32 | — |
| r3 | FOP tag | `ck_spherical = false` | 20.71 s | 14.49 | -1.1 % |

The production default is now *faster* than the `a544a761` control of `bisect_cost`
(23.36 s, 12.84 cycles/s): 15.26 cycles/s, 1.19x.  Most of that is the `ck_sweep_form =
tm` default (`0172f3bc`); the tag adds the last 9.3 %.  The plane-parallel arm gains 1 %
in time although it loses 14 % of its instructions — that instantiation is not the one
the production mesh runs, and it is not latency-bound in the same way.

## 5. Audit of the other switches that reach `rt_chain_ck`

| switch | how it enters the kernel | verdict |
| --- | --- | --- |
| `ck_spherical` | `SPH`, compile-time tag | already right |
| `ck_beam_sph` | `BSP`, compile-time tag | already right |
| `ck_sweep_cache` | `CCH` (0/1/2), compile-time tag | already right |
| `ck_sweep_form` | `FRM` (0/1/2), compile-time tag — *not* a runtime int | already right |
| `ck_impl_frozen_op` (+`_cof`) | was a runtime bool; now `FOP` | **converted here** |
| `ck_impl_reuse_jac` | `ckjacp_`, runtime bool -> `jck`, gates the tridiagonal assembly inside the sweep | see below |
| `ck_impl_colskip` | `ckskip_`, runtime bool, ONE early-return test per column | see below |
| `rt_layer_legacy` | `layer_legacy`, runtime bool selecting the whole staggered sweep vs the centre-to-centre one | see below |
| `ck_impl_tau_min`, `ck_impl_arat`, `ck_impl_seed` | never reach `rt_chain_ck`: they are read in the `ck_impl_res` / `ck_impl_tri` kernels of `two_stream_column_ck.hpp`, which are only launched under `ck_implicit` | no cost when off |

`ck_impl_tau_min` and `ck_impl_arat` are per-cell discriminants and must stay runtime
Reals in any case.

PROBE (a throwaway build, `bench/cktag_0922/athena.probe`, not part of this tree): the
three remaining runtime bools forced to `constexpr false` together, measured on the same
instantiations.

Tier 136, matched instantiations (full table `probe.counts.txt`):

| SPH | BSP | CCH | FRM | FOP | this tree | probe (all three off at compile time) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0 | 2 | 0 | 0 | 24270 instr / 1815 scratch | **12339 / 757** (-49 % / -58 %) |
| 1 | 0 | 0 | 0 | 0 | 27974 / 1686 | 20393 / 1148 (-27 % / -32 %) |
| 1 | 1 | 2 | 0 | 0 | 24563 / 1660 | 12636 / 700 (-49 % / -58 %) |
| 1 | 1 | 0 | 0 | 0 | 28227 / 1663 | 20803 / 1170 (-26 % / -30 %) |

So the three of them together are **larger than everything the frozen-op tag removed**:
about half of the cached spherical kernel is code that only `ck_implicit` (the in-sweep
tridiagonal assembly behind `ckjacp_`/`jck`) or `rt_layer_legacy` (the whole staggered
sweep, kept as an `if/else` beside the centre-to-centre one) can ever execute, and both
are off in production.  `ck_impl_colskip` is one branch and cannot be more than noise in
that total.

Two caveats before anyone acts on it: (i) the probe forces all three at once, so the
split between `ckjacp_` and `layer_legacy` is not measured; (ii) with the branches gone
several instantiations become identical code and are folded (60 symbols instead of 112),
which is why only `FRM = 0` rows survive the match — the numbers above are the ones that
still pair up unambiguously.  **No timing arm was run for the probe**, and the same
lesson as in section 3 applies: instruction count is not time.  The natural next step is
an `IMP` tag (0 = no implicit / 1 = implicit, with `FOP` nested inside it so the
instantiation count goes 2x -> 3x rather than 4x) plus a `LEG` tag for the legacy layers,
gated on a timing job of its own.

## 6. Files

- `gate.sh`, `md5.txt` — the bitwise gate and its checksums; `*_ref` / `*_new` keep the
  `run.log`s only.
- `dis.sh`, `count.py`, `head.counts.txt`, `tag.counts.txt` — the code-object counts.
- `probe.counts.txt` — the section-5 probe (its build directory was deleted afterwards;
  `bench/cktag_0922/src_snapshot` holds this tree, not the probe source).
- `gpu.log.11941488` — the timing job; scripts and binaries in `bench/cktag_0922/`.
