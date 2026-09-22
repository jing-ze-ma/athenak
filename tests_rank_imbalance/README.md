# Rank imbalance of the ck sweep, cs deep hot Jupiter, 2 GPUs (2026-09-22)

The question: why rank 1's two-stream sweep costs ~20 % more than rank 0's in the prod4 form
(`ck_sweep_form = 1`, `ck_spherical`, `ck_beam_sph`; `tests_tm_prof_growth/README.md` (7)),
and whether a cheap change to the rank assignment fixes it.

**Verdict.**
1. The imbalance comes from the pseudo-spherical beam's *twilight* (mu0 < 0, grazing) columns
   on the terminator panels. It is **not** a dayside surplus.
2. The GPU sweep kernel's time is **not additive in MeshBlocks**. It depends strongly and
   non-monotonically on how many blocks, and which ones, a rank holds.
3. An explicit contiguous split of **10 / 14** blocks balances the sweep, which was
   11.7 / 14.4 s at 12 / 12 and is 10.7 / 10.6 s at 10 / 14. It cuts the wall time by
   **12.6 %**: 12.87 -> 14.73 cycles/s, same node, same binary.
4. The dumps are **bitwise identical** for every split tested.
5. The user-proposed leading/trailing split does divide cleanly at this layout. But it
   needs non-contiguous rank ownership, which AthenaK does not support, so it was not
   measured.

## Setup

* Snapshot: `git archive HEAD` (`3b6e0e52`) in `bench/rank_imb_0922/snap` (kokkos
  symlinked), plus the patch `bench/rank_imb_0922/lb_nmb_eachrank.patch` (not committed).
* Binary: `bench/rank_imb_0922/athena_gpu`, md5 `bf15148f0665f6a6d065ba2445694f4d`.
  Recipe as in `tests_mpi_gate_0922`: gcc/14 rocm/6.3 openmpi_gpu/5.0, hipcc,
  `Kokkos_ARCH_AMD_GFX942_APU`, Release, MPI ON, `PROBLEM=deep_hot_jupiter_rt`.
* Every arm uses the same binary and the same input as `tm_prof_growth/prof_tm`
  (prod4 form).
* Every arm restarts from `bench/cs_mhd_prod3/rst/dhj.00567.rst` and runs 300 cycles
  (nlim 4430131) on apudev with 2 ranks / 2 GPUs.
* `P*`/`Q*` arms run under `rocprofv3 --kernel-trace --stats`. `U*`/`V*`/`S*` arms are
  unprofiled.
* Jobs 11942109 (vipa1001), 11942119, 11942154 and 11942167 (vipa1327). Per-arm run dirs
  are in `bench/rank_imb_0922/runs/`.
* Scripts: `submit.sh` (the last job; earlier versions in `submit*.bak`), `sumprof.py`
  (per-rank kernel sums) and `geom.py` (column census). The bitwise check is
  `bench/mpi_gate_0922/cmpbin.py`.

### The switch (prototype, default OFF)

`<mesh>/lb_nmb_eachrank = 10,14` sets an explicit MeshBlock count per rank. Blocks stay
contiguous in gid, so only the split point moves.

* It is read in the `Mesh` ctor and applied in `Mesh::LoadBalance`, so it works from
  scratch and on restart.
* It is ignored with a warning unless it has one entry per rank, every entry is >= 1 and
  the entries sum to the block count.
* The parameter must be in the input file. A command-line override alone fails with
  "parameter not found", because the restart's embedded input does not have it.
* cpplint reports no new errors in the patched lines.

## (a) Is the current assignment day-heavy on one rank? No: rank 1 holds the terminators.

`gid` is panel-major, 4 blocks per panel. The panels are 0 = +x (substellar,
`mu0 = cx`), 1 = +y, 2 = -x (antistellar), 3 = +z, 4 = -y and 5 = -z. With 12 / 12:

* **rank 0** = panels 0, 1, 2: the whole substellar panel, one terminator panel, and the
  antistellar panel, whose centre is fully dark;
* **rank 1** = panels 3, 4, 5: three terminator panels.

Column census from `geom.py` (the rcut value is assumed, 1e10 cm):

| rank | day (mu0 > 0) | twilight (mu0 < 0, lit at top) | dark |
|---|---|---|---|
| 0 | 6 blocks (4 at mean mu0 +0.82, 2 at +0.34) | 2 blocks at -0.34, plus 4 x 168 cols at -0.82 | 352 cols |
| 1 | 6 blocks at +0.34 | 6 blocks at -0.34 | 0 |

The dayside is split evenly. What rank 1 has in excess is terminator twilight. There the
`ck_beam_sph` ray is lit deep, and each face runs the tangent search plus the two-leg path
integral, which is O(N^2) in the column.

Measured with beam on/off in the same binary on vipa1001 (per-rank sum of the
`picket_fence_*` kernels, ms per 300 cycles):

| arm | rank 0 sweep | rank 1 sweep | r1/r0 |
|---|---|---|---|
| P12, beam_sph on (prod4) | 11441 | 12992 | 1.136 |
| PB0, `ck_beam_sph = false` | 10091 | 9267 | **0.918** |
| beam cost (difference) | 1350 | **3725** | 2.8 |

Without the spherical beam, rank 1 is the *lighter* rank, matching the old 4-pass
observation. The spherical beam puts 2.8x as much work on rank 1 as on rank 0.

No per-block timer was built. The per-block picture comes from the census and the split
scan below, and the scan shows that per-block costs are not additive anyway.

## (b) The leading/trailing split

The substellar-antistellar meridian plane is y = 0.

* Panels 1 (+y) and 4 (-y) lie wholly on one side each.
* The plane cuts panels 0, 2, 3 and 5 exactly along their x2 mid-line, which is a block
  boundary at the production 2 x 2 per-panel layout. **No block straddles it**, and the
  hemispheres are 12 blocks each.
* Every column has a mirror partner (y -> -y) with the same mu0, so the beam geometry is
  identical on both sides.

It is still not cheap. Within each cut panel, the Z-ordered gids alternate between the two
hemispheres (for panel 0: gids 0, 2 trailing; 1, 3 leading). The split therefore needs
**non-contiguous** rank ownership. AthenaK assumes one contiguous gid range per rank:
`gids_eachrank` and `nmb_eachrank` are used at 67 sites in 20+ files, including the
binary and restart output offsets and bvals. The alternative is a gid renumbering, which
changes the dump byte order and the restart mapping. Neither is a prototype-sized change,
so this split was **not measured**.

The scan below also shows that balancing geometric cost does not predict GPU time here.
A mirror-balanced split could still land on a slow configuration.

## (c) Per-rank sweep and cycles/s, same binary, same node (vipa1327)

Profiled, main sweep kernel plus all `picket_fence_*` kernels, ms per 300 cycles:

| arm | split | GPU map | r0 sweep | r1 sweep | r0 total | r1 total |
|---|---|---|---|---|---|---|
| Q12 | 12/12 (HEAD) | 0,1 | 11652 | **14433** | 16511 | **19367** |
| Q12s | 12/12 | swapped | 12125 | 12818 | 17007 | 17757 |
| Q10 | **10/14** | 0,1 | 10690 | 10562 | 15304 | 16318 |
| P13 | 13/11 | 0,1 | 12534 | 13072 | 17995 | 17921 |
| P11 | 11/13 | 0,1 | 11125 | 14888 | 15894 | 20407 |
| P14 | 14/10 | 0,1 | 10524 | 11706 | 16263 | 16376 |

In the dominant sweep kernel (`picket_fence...<136>`), rank 1 takes 13687 ms with gids
12-23 but **9671 ms with gids 10-23**. Adding two cheap antistellar blocks at the front of
its pack made the whole launch 29 % faster. The kernel time is set by
scheduling/occupancy on the GPU, not by the sum of column work.

The swapped-GPU arm also shows a hardware or placement component. The same 12 blocks run
4-13 % slower on GPU 1 than on GPU 0.

Unprofiled `cpu time used` (loop wall), 300 cycles:

| split | runs (s) | mean (s) | cycles/s | vs 12/12 |
|---|---|---|---|---|
| 12/12 (HEAD) | 23.15, 23.34, 22.65, 22.99, 24.17, 23.60 | 23.32 | 12.87 | - |
| **10/14** | 20.39, 20.50, 20.28, 20.31 | **20.37** | **14.73** | **+14.5 %** |
| 8/16 | 21.62 | 21.62 | 13.88 | +7.9 % |
| 14/10 | 21.26, 21.16 | 21.21 | 14.14 | +9.9 % |
| 16/8 | 21.98 | 21.98 | 13.65 | +6.1 % |
| 13/11 | 23.32, 23.56, 23.59, 22.46, 22.81 | 23.15 | 12.96 | +0.7 % |
| 11/13 | 24.18 | 24.18 | 12.41 | -3.6 % |
| 9/15 | 26.85 | 26.85 | 11.17 | -13 % |
| 15/9 | 28.65 | 28.65 | 10.47 | -19 % |

Two more 12/12 runs on vipa1001 took 21.90 s and 22.23 s. 10/14 was not run on vipa1001.

The response to the split point is non-monotonic. Odd per-rank counts (9, 11, 13, 15) are
bad and 10/14 is best, which is not what a smooth cost model gives. A likely cause is
launch geometry: the workgroup-to-XCD mapping of the heavy (m, chain) workgroups on
MI300A. This is not verified.

**Bitwise check.** The decoded payload of the final dump (`dhj.mhd_w_bcc.00143.bin`,
6291456 values) is identical to 12/12 for every split tested:

* 10/14, 14/10, 13/11, 11/13, 8/16, 16/8, 9/15, 15/9;
* 12/12 with the GPUs swapped;
* the last `.hst` row is also identical.

File md5s differ only because the header embeds the input text, which now has the
`lb_nmb_eachrank` line. The 12/12 md5 `b2f14b64...` is the same in every run.

## Recommendation

* For the 2-GPU prod4 layout (24 blocks, 2 x 2 per panel, this binary), put
  `lb_nmb_eachrank = 10,14` in `<mesh>`: about +14 % cycles/s at zero change to the
  answer.
* The number is empirical and specific to this block layout and kernel build. Re-scan
  (about 30 s per arm on apudev) after any change to the sweep kernel, the block size or
  the rank count.
* Not tested: writing a restart under the override. Restarts are gid-ordered and
  re-balanced on read, so this is expected to be safe.
* The larger follow-up is inside the kernel. Rank 1's 14-block launch runs at ~1.15 ms
  per block per call against ~1.90 ms at 12 blocks, which suggests the 12-block launch is
  badly scheduled. Reordering the chain/workgroup mapping so the heavy twilight beam work
  does not bunch up would help every split, and would not depend on a tuned split.
