# MeshBlock-layout invariance suite

Checks that results do not depend on the MeshBlock decomposition or on the MPI rank count.
It was written after the cubed-sphere seam-halo bug of 11-03, in which the seam stencils were
clamped to the sending MeshBlock instead of the panel. That bug went unnoticed for months
because conservation still held.

## How it works

Each case in `cases.py` runs the same problem several times from identical fresh analytic ICs:

- the **reference layout** (the first one, normally one MeshBlock per panel);
- the reference layout again with a **1-ulp perturbation** (default `time/cfl_number` moved
  to the next double);
- every **other layout**: more blocks per panel, a phi or theta split, a radial split,
  1 to 8 MPI ranks.

Every run writes double-precision `.bin` output (`bin_double`, `dcycle = 1`).
`libin.py` assembles the MeshBlocks into global `(panel, k, j, i)` arrays. Cubed-sphere gids
are panel-major. The comparison is done cell by cell for every listed variable at cycle 1 and
at cycle `ncyc`.

The 1-ulp run measures how far plain round-off spreads by that cycle. The verdict for each
(layout, cycle, variable) is one of:

| verdict  | relative max difference vs the reference layout            |
|----------|-------------------------------------------------------------|
| BITWISE  | 0                                                           |
| ROUNDOFF | <= max(10 x ulp-ref, 100 eps)                               |
| SUSPECT  | in between                                                  |
| FAIL     | > max(1000 x ulp-ref, 1e-10)                                |

For SUSPECT and FAIL the report also gives the cell with the largest difference, and the
fraction of differing cells that lie within `ng` of an INTERNAL block end of the test layout.
A value near 100 % is the signature of this bug class. For scale, the 11-03 clamp gives
about 1e-5 after one cycle, where round-off is about 1e-16.

A case can carry `expect = 'FAIL'`. This keeps a known layout dependence visible: the
`cs_seam_blockend = clamp` default and the r1m clamp deck are both kept this way. The summary
flags any case whose outcome differs from its expectation.

## Running

```
# bins.cfg: "family = /path/to/athena" for bip (built_in_pgens), cs (cs_test), sp (sp_test),
#           dhj (deep_hot_jupiter_rt), he (he_star_m1), acc (ry_per_accretor); CPU builds, MPI on.
python3 li.py --list
python3 li.py --bins bins.cfg --work /scratch/li_runs [--only case1,case2] [--np-max 16]
python3 li.py --work /scratch/li_runs --compare-only     # re-compare existing runs
```

Exit status is 0 when every case meets its expectation. The report is printed and written to
`<work>/REPORT.txt`.

On Caltech, use `run_li.sh`. It loads gcc 13.2 and openmpi 5.0.1. Do not set
`OMPI_MCA_io=romio321`: openmpi 5 has no such component. The binaries are listed in
`bins.caltech.cfg`. Launch ON a compute node of a held allocation, detached:

```
ssh <node> "setsid nohup bash -l .../run_li.sh <work> --only ... > <work>.out 2>&1 < /dev/null &"
```

The dhj cases need the WASP-121b packages, the ck tables and the IC profile (env
`LI_DHJ_PKG`, default `/resnick/groups/carnegie_poc/jingze/w121r1_1101/pkg`). The whole suite
takes about 15 min on 16 cores; the cubed-sphere and Cartesian cases take seconds each.

## Coverage

- **Cartesian** (`built_in_pgens`): 3-D hydro and MHD linear waves (CT edges shared by 4 or 8
  blocks), 2-D implode with reflecting walls (physical-BC corners), Orszag-Tang, and an M1
  beam with user/vacuum boundaries.
- **Cubed sphere** (`cs_test`):
  - hydro and MHD blasts through a seam and a cube vertex, with clamp and with extend;
  - nghost 3 with ppm4, the HG/RSG halo depth;
  - the r1m seam keys (`cs_seam_flux = shared`, `cs_seam_mhd = shared_o`, `uct_hlld`);
  - Ohmic plus ambipolar diffusion with RKG STS;
  - scalar diffusion and a passive scalar;
  - probes at C8/C16/C32 that separate block size from blocks per panel.
- **Spherical polar** (`sp_test`):
  - hydro and MHD blasts across the polar boundary, with the phi direction split into 2, 4
    and 8 blocks, a theta split, and a radial split for hydro;
  - `polar_ring_avg` (auto group sizes), the ring Fourier filter, `polar_emf_diss` and
    `polar_hlle_rows`;
  - resistive decay with Ohmic, ambipolar and RKG.
- **dhj** (`deep_hot_jupiter_rt`), using the production r1m v6 MHD deck and the hydro planet
  deck on C8 x 32:
  - physics: well-balance cache, top sponge, Roche tide plus the `wb_tan` probe, `outflow1`
    and `hse` top BCs, bottom sponge, floors;
  - ck: the implicit ck two-stream (spherical beam) at cadence 4 and 1, plus a
    `ck_impl_jreuse = 0` control;
  - MHD: RKG Ohmic plus ambipolar diffusion and tracers;
  - layouts: 1 vs 2x2 blocks per panel and 1, 3, 4 or 6 ranks. dhj refuses a radial split.

## Not covered yet

- Restart starts from a permutation remap, i.e. production states. Fresh analytic ICs only.
- he_star_m1 (the HG/RSG pgen): M1 with vet_sc, vet_gd and vetlat.
- ry_per_accretor (the sp accretor deck).
- Refined meshes: `libin.assemble` handles single-level meshes only.
- GPU builds: run the same cases with GPU binaries in `bins.cfg`. The suite is
  backend-agnostic, but it uses mpirun on one node.
