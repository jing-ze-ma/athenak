# MPI-rank and multi-node invariance gate for the next cs dhj production (2026-09-22)

Gate of the changes intended for the next cubed-sphere deep-hot-Jupiter production, none of
which had been checked for rank-count or multi-node invariance:

| change | commit | default | covered by |
| --- | --- | --- | --- |
| cs corner-slot seam resample (`bvals_cc`, buffer slots 48-55) | `ec0a9724` | ON | A1, A2, B |
| `mesh/cs_vertex_fill_cc` | `mesh.cpp:143` | OFF | A1, A2 (off **and** on) |
| restart caches | `0eb7e2c4`, `544cb36b` | ON | A4, B3 |
| `problem/ck_spherical`, `ck_beam_sph` (column-local) | — | on / off | A3, B |
| `hlld_bx_zero_tol` | `0aaf7f39` | 1e-4 | A1, B |
| `rad_angular` | — | — | compiled in, exercised by B |

**Snapshot and binaries.** `git archive HEAD` of `rt-integration`
`446b7b30c3a049009b2facdf34e19c09e711f0d8` into
`/viper/u2/jinma/ATHENAK/bench/mpi_gate_0922/snap` (the live working tree is dirty and was
never used for a build; `kokkos/` is symlinked). Three binaries, all
`CMAKE_BUILD_TYPE=Release`, `Athena_ENABLE_MPI=ON`:

| binary | modules | compiler | md5 |
| --- | --- | --- | --- |
| `athena_cpu_cs_test` (`PROBLEM=cs_test`) | gcc/14 openmpi/5.0 cmake/4.0 | `mpicxx` | `be6f82a843827d019d6b939cbd1e7026` |
| `athena_cpu_deep_hot_jupiter_rt` | gcc/14 openmpi/5.0 cmake/4.0 | `mpicxx` | `22a0f969ee1207ca317c77b4bba7ecef` |
| `athena_gpu` (HIP, `Kokkos_ARCH_AMD_GFX942_APU`) | gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0 | `hipcc` | `78724a0106241802ac4618f1955e53c4` |

`strings -a athena_gpu` confirms `ck_spherical`, `ck_beam_sph`, `hlld_bx_zero_tol`,
`cs_vertex_fill_cc`, `rad_angular`, `ck_implicit` are all compiled in.
The GPU md5 differs from `bench/prof_0922`'s `2f2c07df...` (that build is HEAD `a544a761`),
so `prof_0922/plain` was **not** reused as the 2-rank reference; it was rerun (arm `g2a`).

**What "rank-invariant" can mean here.** Read of the source first, because the criterion
depends on it:

* History global sums are per-rank Kokkos `parallel_reduce` partials combined by a plain
  `MPI_Reduce(..., MPI_SUM, ...)` — `src/outputs/history.cpp:414-419`. There is **no**
  fixed-order or per-gid deterministic summation anywhere in that path, so `.hst` columns
  are only expected to agree to summation round-off, not bitwise.
* The timestep reduction is `MPI_Allreduce(..., MPI_MIN, ...)` — `src/mesh/mesh.cpp:1005`.
  `MIN` is order-independent, so `dt` *is* exactly rank-invariant.
* Binary dumps are written contiguously in global `gid` order
  (`src/outputs/binary.cpp:8`, offsets at `:234-244`), so their byte layout does not depend
  on the rank count. **Dumps are therefore the valid bitwise gate**; `.hst` is the weaker one.
* Restarts default to a single shared file (`single_file_per_rank`, default `false`,
  `src/outputs/outputs.cpp:245`), and `Mesh::BuildTreeFromRestart` re-load-balances onto the
  *current* rank count (`src/mesh/build_tree.cpp:420-465`), so N-written / M-read is legal.

---

## A. CPU, exact

Login node `viper13`, `mpirun -n 1 / 2 / 4`, all runs from the same input, `nlim = 20`.
Comparison script `cmp.py` (md5 of each `.bin`, plus every `.hst` column) and `cmpbin.py`
(decoded per-cell `float64` payload via `vis/python/bin_convert.py`, used where header text
legitimately differs). Run tree: `bench/mpi_gate_0922/cpu/`.

### A1-A2. cs blast, 24 MeshBlocks (2x2 blocks per panel)

`inputs/tests/cubed_sphere_mhd_blast.athinput` and `cubed_sphere_blast.athinput` with
`<meshblock>/nx2 = nx3 = 16` against `<mesh>/nx2 = nx3 = 32` -> 6 panels x 4 = **24
MeshBlocks** (verified: `MeshBlock-cycles = 480` = 24 x 20). Blast centre is the cube vertex
`(1,1,1)`, i.e. the seam-and-vertex-hardest placement, and `hlld` with the new
`hlld_bx_zero_tol` runs on the MHD arm.

| test | `cs_vertex_fill_cc` | 1 vs 2 ranks | 1 vs 4 ranks |
| --- | --- | --- | --- |
| cs MHD blast | off (default) | dumps **bitwise**, all 13 hst cols **bitwise** | dumps **bitwise**, all 13 hst cols **bitwise** |
| cs MHD blast | **on** | dumps **bitwise**, all 13 hst cols **bitwise** | dumps **bitwise**, all 13 hst cols **bitwise** |
| cs hydro blast | off (default) | dumps **bitwise**, all hst cols **bitwise** | dumps **bitwise**, all hst cols **bitwise** |
| cs hydro blast | **on** | dumps **bitwise**, all hst cols **bitwise** | dumps **bitwise**, all hst cols **bitwise** |

Both the `t = 0` dump and the cycle-20 dump are byte-identical (raw md5 of the whole file,
header included) at every rank count.

**The `cs_vertex_fill_cc = on` arm is not a no-op** — control: at fixed 1 rank, vf0 vs vf1
differ (`1-ME` rel 1.15e-04, `3-ME` rel 5.13e-06, dumps differ), so the ON arm really
exercises the widened cell-centred corner buffers, and it is rank-invariant.

So `ec0a9724`'s corner-slot seam resample and the cell-centred cube-vertex fill both
produce identical MPI exchanges however the 24 blocks are split across ranks.

### A3. dhj CPU test, 8 MeshBlocks

`inputs/tests/dhj_ck_spherical.athinput` does not decompose into >= 4 blocks as shipped
(`<meshblock>/nx3 = 4` against `<mesh>/nx3 = 8` gives 2, and the shipped `nx3 = 2` split is
rejected by `mesh.cpp:581`, "MeshBlock must be >= 4 cells in each active dimension").
Variant used: `<mesh>/nx2 = 8, nx3 = 16`, `<meshblock>/nx1 = 128, nx2 = nx3 = 4` -> **8
MeshBlocks** (`MeshBlock-cycles = 160` = 8 x 20), with `problem/ck_spherical = true` so the
new spherical correlated-k two-stream is the one under test.

| quantity | 1 vs 2 ranks | 1 vs 4 ranks |
| --- | --- | --- |
| `t = 0` and cycle-20 `bin` dumps (per cell, float64) | **bitwise identical** | **bitwise identical** |
| `.hst` cols `time dt mass 1-mom 3-mom tot-E 1-KE 2-KE 3-KE` | **bitwise** | **bitwise** |
| `.hst` col `2-mom` | rel 5.3e-01 | rel 4.3e-01 |
| `dhj.log` event counters | identical (all zero over 20 cycles) | identical |

**The `2-mom` difference is a reduction-order effect, not an exchange defect**, and the
dumps prove it: every cell of the state is bit-identical, so the only thing that can differ
is the summation. `2-mom` is the theta momentum, which cancels by symmetry — it sits at
~2e13 against `3-mom` ~1.2e27 and `mass` ~3.5e26, i.e. **14 decades below the other
components of the same sum**. A catastrophically cancelling sum reordered by
`MPI_Reduce(MPI_SUM)` (history.cpp:414-419) is exactly this, and the `.hst` file only prints
6 significant digits. No other column moves at all.

### A4. Restart invariance (CPU)

Write 10 cycles at N ranks, restart the resulting `rst`, run to cycle 20 at M ranks, compare
the final state against the uninterrupted 1-rank 20-cycle run, per cell in float64.

| test | write ranks -> restart ranks | restart-leg `t=0` dump | final (cycle 20) state |
| --- | --- | --- | --- |
| cs MHD blast, `cs_vertex_fill_cc = on` | 1 -> 4 | **bitwise** | **bitwise** |
| cs MHD blast, `cs_vertex_fill_cc = on` | 4 -> 1 | **bitwise** | **bitwise** |
| dhj (`ck_spherical = true`) | 1 -> 4 | **bitwise** | **bitwise** |
| dhj (`ck_spherical = true`) | 4 -> 1 | **bitwise** | **bitwise** |

So the restart caches `0eb7e2c4` / `544cb36b` restore correctly at a rank count different
from the writer's, in both directions. (The raw file md5 of the restarted run's dump differs
from the uninterrupted run's only because the binary header embeds the effective input text,
which legitimately records different parameters after a restart; the decoded payload is
identical.)

### A5. DEFECT: restarting from a one-file-per-rank restart file aborts

Not one of the gated changes, but found while covering the `rank_%08d/` layout mentioned in
`src/main.cpp:235-251`.

Writing works: with `<output>/single_file_per_rank = true` a 4-rank run correctly produces
`rst/rank_00000000/ .. rank_00000003/`, each with its own `.rst`. **Reading back fails
immediately, at any rank count:**

```
[viper13:00000] *** An error occurred in MPI_File_get_position
[viper13:00000] *** MPI_ERR_INTERN: internal error
[viper13:00000] *** MPI_ERRORS_ARE_FATAL
```

Cause, located exactly:

```
src/mesh/build_tree.cpp:426  void Mesh::BuildTreeFromRestart(ParameterInput *pin, IOWrapper &resfile,
src/mesh/build_tree.cpp:427                                                       bool single_file_per_rank) {
src/mesh/build_tree.cpp:429    IOWrapperSizeT headeroffset = resfile.GetPosition();
```

`GetPosition()` is called **without** the `single_file_per_rank` argument the function
already has in hand, so it takes the default `false`
(`src/outputs/io_wrapper.hpp:59`) and dispatches to `MPI_File_get_position(fh_, ...)`
(`io_wrapper.cpp:615-620`) on a handle that in the per-rank path is a plain `FILE*` opened
with `fopen` (`io_wrapper.cpp:42`). The two other call sites do pass the flag correctly
(`src/main.cpp:264`, `src/pgen/pgen.cpp:395`), so `build_tree.cpp:429` is the single
offender; `resfile.GetPosition(single_file_per_rank)` is the one-token fix.

**Impact on the production: none.** `cs_mhd_prod3` and every arm here write the default
single shared restart file, which is the layout gated in A4 and B3 and which works. Note
also that even once fixed, the per-rank layout is by construction tied to the writer's rank
count (`main.cpp:249-251` rebuilds the path from `my_rank`), so it can never give N -> M
restart invariance. Not fixed here (the repo was not to be touched).

---

## B. GPU, HIP + MPI

`athena_gpu`, restart `bench/cs_mhd_prod3/rst/dhj.00567.rst` (rot 283, ncycle 4429831),
input = `bench/prof_0922/plain/deep_hot_jupiter.athinput` (the production cs config:
24 MeshBlocks of 128x16x16, MHD + general EOS + correlated-k RT, `ck_spherical = true`,
`cs_gs07_emf = true`, `cs_lowbeta_fallback = 0.5`) with `nlim = 4430131` = **300 cycles**.
Every arm terminated on the cycle limit with `MeshBlock-cycles = 7200` = 24 x 300.

| arm | job | layout | partition | wall | cpu time | zone-cycles/s | per GPU | parallel eff. |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `g1` | 11933697 | 1 rank / 1 GPU / 1 node | apudev | 55.7 s | 47.46 s | 4.972e6 | 4.972e6 | 1.00 (ref) |
| `g2a` | 11933694 | 2 ranks / 2 GPUs / 1 node (**= production**) | apudev | 39.2 s | 29.81 s | 7.913e6 | 3.957e6 | 0.796 |
| `g2b` | 11933695 | 2 ranks / 2 GPUs / 1 node (repeat, for the spread) | apudev | 36.8 s | 29.72 s | 7.939e6 | 3.969e6 | 0.798 |
| `g4` | 11933699 | **4 ranks / 2 nodes**, 2 GPUs+2 tasks per node | apu1 | 37.4 s | 23.31 s | 1.012e7 | 2.530e6 | 0.509 |

The 1-rank arm was **not** skipped: 24 MeshBlocks of 128x16x16 is 786k cells, trivially
inside one MI300A's unified memory, and it ran (55.7 s, the slowest arm as expected).
`apudev` has `MaxNodes=1`, so the 2-node arm went to `apu1` (user-approved), `--nodes=2
--ntasks=4 --ntasks-per-node=2 --gres=gpu:2`; script verified with
`scontrol write batch_script 11933699`.

### B1. Same-binary run-to-run spread — there is none

The criterion was "differences at the level of the same-binary spread". Measured on the
`g2a` / `g2b` pair, at full float64, on the restart file each arm writes at finalize
(`rst/dhj.00568.rst`, 112 605 006 bytes of payload after `<par_end>`):

**the two 2-rank runs of the same binary are byte-for-byte identical after 300 cycles.**
The spread is exactly zero, so the gate becomes the strictest possible one.

### B2. Rank-count and multi-node invariance — bitwise

| comparison | float64 restart payload | `mhd_w_bcc` bin dump (per cell) | `dhj.log` counters | NaN / divergence |
| --- | --- | --- | --- | --- |
| `g2b` vs `g2a` (2r vs 2r) | **bitwise identical** | all 8 vars identical | identical | none, all finite |
| `g4` vs `g2a` (**4 ranks / 2 nodes** vs 2 ranks / 1 node) | **bitwise identical** | all 8 vars identical | identical | none, all finite |
| `g1` vs `g2a` (1 rank vs 2 ranks) | **bitwise identical** | all 8 vars identical | identical | none, all finite |

`.hst` (one row, written at finalize — the restart's stored output cadence means the
`3.05e3` / `1.0e2` interval does not re-trigger inside 300 cycles) agrees in all 13 columns
at every arm: `t = 8.640522e7`, `dt = 1.61819e1`, `mass = 3.55308e26`, `tot-E = 6.21521e38`.

Event counters, identical on all four arms (cumulative from the restart, so the increments
match too):

```
 cycle    eos_dfloor eos_efloor eos_tfloor eos_vceil eos_fail c2p_it fofc efloor_de eos_tclamp eos_tset vceil_de
 4430131  41672216   2292792    91292      0         0        0      0    0.0       8          0        0.0
```

Zero `eos_fail`, zero `fofc`, zero `eos_vceil` on every arm.

**This updates the `rt-integration` memory note** ("same-binary spread ~1e-2 by 200
cycles"): at HEAD `446b7b30`, in this configuration, there is no run-to-run spread at all
over 300 cycles, and the results are bitwise across 1 / 2 / 4 ranks and across nodes. The
`MPI_MIN` dt reduction plus dumps written in gid order leaves nothing rank-dependent, and no
kernel in this path introduces atomics-order nondeterminism.

### B3. Restart cross-check (GPU)

* **2-rank-written restart read by 4 ranks on 2 nodes**: this *is* arm `g4` —
  `dhj.00567.rst` was written by the 2-rank production chain, and `g4` restarted it on 4
  ranks across 2 nodes and reproduced `g2a` bitwise over 300 cycles. Passed.
* **4-rank-written restart read by 2 ranks**: `gr2` (job 11933978, apudev, 2 ranks) and
  `gr4` (job 11933979, apu1, 4 ranks / 2 nodes) both restart `g4`'s own
  `gpu/g4/rst/dhj.00568.rst` and run 50 more cycles (`nlim = 4430181`).

| arm | job | layout | partition | wall | cpu time | cycles | zone-cycles/s | per GPU |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `gr2` | 11933978 | 2 ranks / 2 GPUs / 1 node | apudev | 14.5 s | 5.341 s | 50 | 7.362e6 | 3.681e6 |
| `gr4` | 11933979 | **4 ranks / 2 nodes** | apu1 | 18.3 s | 4.441 s | 50 | 8.853e6 | 2.213e6 |

Both terminated on the cycle limit at `cycle = 4430181`, `MeshBlock-cycles = 1200` = 24 x 50
(the small `cpu time` is by design — these are 50-cycle legs, not 300). No NaN, no fatal,
`rc = 0` on both.

| comparison | whole-file md5, `rst/dhj.00569.rst` (112 626 809 B) | `bin/dhj.mhd_w_bcc.00144.bin` (25 189 989 B) | `.hst` (13 cols) | `dhj.log` counters |
| --- | --- | --- | --- | --- |
| `gr4` (4 ranks / 2 nodes) vs `gr2` (2 ranks / 1 node) | **byte-identical** (`7b922 7cf0…`, `cmp` rc 0) | **byte-identical** (`89108 4657…`, `cmp` rc 0) | identical to the last printed digit, max rel diff **0.0** in every column | identical |

Counters on both restart legs (incremental from `dhj.00568.rst`, `eos_tclamp` reset to 0):

```
 cycle    eos_dfloor eos_efloor eos_tfloor eos_vceil eos_fail c2p_it fofc efloor_de eos_tclamp eos_tset vceil_de
 4430181  6601670    609067     35819      0         0        0      0    0.0       0          0        0.0
```

So a restart file written by a 4-rank / 2-node run is read back identically at 2 ranks on one
node and at 4 ranks on two nodes, and the two continuations stay bitwise for another 50
cycles. **Caveat, stated because it weakens the label:** `g1`, `g2a`, `g2b` and `g4` all write
a *byte-identical* `dhj.00568.rst` (md5 `b477793a3e3baf6775690204affcde62` on all four), so
"4-rank-written" and "2-rank-written" are literally the same bytes here. The writer side is
therefore proven rank-independent (that md5 agreement *is* the proof), and `gr2`/`gr4` test
the reader side at a rank count different from the writer's. What is **not** covered on the
GPU is restart-versus-uninterrupted: every GPU arm starts from a restart, and no 350-cycle
continuous arm was run. That leg is covered on the CPU by A4, in both directions.

### B4. Strengthened statement of B1-B2

The comparison in B1/B2 was made on the decoded float64 payload; the raw whole-file md5s are
in fact identical too, header included, across all four arms:

| file | `g1` | `g2a` | `g2b` | `g4` |
| --- | --- | --- | --- | --- |
| `rst/dhj.00568.rst` | `b477793a…` | `b477793a…` | `b477793a…` | `b477793a…` |
| `bin/dhj.mhd_w_bcc.00143.bin` | `74e14659…` | `74e14659…` | `74e14659…` | `74e14659…` |
| `.hst` max rel diff vs `g2a`, all 13 cols | 0.0 | — | 0.0 | 0.0 |
| `dhj.log` counters | identical | ref | identical | identical |

### B5. Scaling

Per-GPU throughput, 300-cycle arms, 24 MeshBlocks of 128x16x16 (786 432 cells; 24 / 12 / 6
blocks per GPU at 1 / 2 / 4 ranks):

| arm | GPUs | nodes | zone-cycles/s | per GPU | speed-up vs `g1` | parallel eff. |
| --- | --- | --- | --- | --- | --- | --- |
| `g1` | 1 | 1 | 4.972e6 | 4.972e6 | 1.00 | 1.00 |
| `g2a` | 2 | 1 | 7.913e6 | 3.957e6 | 1.59 | 0.796 |
| `g2b` | 2 | 1 | 7.939e6 | 3.969e6 | 1.60 | 0.798 |
| `g4` | 4 | **2** | 1.012e7 | 2.530e6 | 2.04 | 0.509 |

2 ranks on 1 node -> 4 ranks on 2 nodes costs a further **factor 1.56 in per-GPU throughput**
(3.957e6 -> 2.530e6, efficiency 0.796 -> 0.509), i.e. going off-node is about as expensive
again as the first doubling was. The 50-cycle restart legs show the same shape
(`gr2` 3.681e6/GPU vs `gr4` 2.213e6/GPU, ratio 1.66), with a larger start-up share.

This is a strong-scaling curve on a deliberately small problem: at 4 ranks each GPU holds
only 6 MeshBlocks of 128x16x16, so both the per-GPU work and the surface/volume ratio are
unfavourable and the loss is the expected one, not evidence of a defect. Nothing in it is a
correctness statement — B2/B3 already show the answers are bitwise whatever the layout.

---

## Verdict, per item of the gate

| gate item | verdict | evidence |
| --- | --- | --- |
| **Rank invariance, CPU** | **PASS, bitwise.** cs hydro and cs MHD blast (24 blocks, blast on the cube vertex), `cs_vertex_fill_cc` off *and* on, 1 vs 2 vs 4 ranks: dumps byte-identical, all `.hst` columns bitwise. dhj `ck_spherical` (8 blocks): every cell bitwise; only `2-mom` moves in `.hst`, and that is `MPI_Reduce(MPI_SUM)` reordering a sum that cancels 14 decades, proven by the bitwise dumps. | A1-A3 |
| **Restart invariance across rank counts, CPU** | **PASS, bitwise, both directions.** 1 -> 4 and 4 -> 1 for cs MHD (`cs_vertex_fill_cc = on`) and for dhj: restart-leg `t = 0` dump and the cycle-20 final state are bitwise against the uninterrupted 1-rank run. | A4 |
| **Restart invariance, GPU** | **PASS as far as tested.** N-written / M-read is bitwise (`gr2` vs `gr4`, and the writers agree byte-for-byte across 1/2/4 ranks). Restart-vs-uninterrupted is *not* tested on the GPU — no continuous arm exists; it is covered on the CPU. | B3 |
| **Multi-node GPU** | **PASS, far inside the spread — the spread is zero.** `g2a` vs `g2b` (same binary, same layout, 300 cycles) are byte-for-byte identical, so the tolerance collapses to bitwise, and `g4` (4 ranks / 2 nodes) and `g1` meet it. | B1, B2 |
| **Scaling** | Acceptable and understood, not a defect. 0.80 parallel efficiency on-node at 2 GPUs, 0.51 at 4 GPUs across 2 nodes; 6 blocks/GPU at 4 ranks is simply too little work per device. | B5 |
| **Defects found** | **One, unrelated to the gated changes and with no impact on the production:** restarting from a one-file-per-rank restart aborts in `MPI_File_get_position`. Localised to a single missing argument, `src/mesh/build_tree.cpp:429` — `resfile.GetPosition()` should be `resfile.GetPosition(single_file_per_rank)`. Writing that layout works; the two other call sites pass the flag. Not fixed here (repo untouched). | A5 |

**Overall: the gate passes.** The corner-slot seam resample (`ec0a9724`), `cs_vertex_fill_cc`,
the restart caches (`0eb7e2c4`, `544cb36b`), `ck_spherical` / `ck_beam_sph`,
`hlld_bx_zero_tol` (`0aaf7f39`) and `rad_angular` are all bitwise rank-invariant on CPU and
GPU, across nodes, and across restarts at a changed rank count. The next cs dhj production
may run at any of these layouts and at 2 nodes without changing its answer.

---

## Reproducing

```
cd /viper/u2/jinma/ATHENAK/bench/mpi_gate_0922
./build_cpu.sh cs_test ; ./build_cpu.sh deep_hot_jupiter_rt ; ./build_gpu.sh
./run_cpu.sh csmhd_vf1 athena_cpu_cs_test csmhd_vf1.athinput 4   # etc
./run_rst.sh csmhd athena_cpu_cs_test csmhd_vf1.athinput 1 4
python3 cmp.py cpu/csmhd_vf1/r1 cpu/csmhd_vf1/r4 "1 vs 4"
python3 cmpbin.py <a.bin> <b.bin> label
(cd gpu/g4 && sbatch submit.sh)
```

Nothing under `bench/cs_mhd_prod3`, `bench/prof_0922` or `athenak/run` was written to; the
repo working tree was not modified (this README is the only file added).
