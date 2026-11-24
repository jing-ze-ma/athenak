# 3-D star tests: spherical polar and cubed sphere

Two small, complete 3-D star decks that use the production scheme keys. Each one runs
in about 15 minutes on 4 A100 GPUs, and it starts on one CPU node, where a few cycles show
that the setup works. Use them to check a new build or machine before you start a production
run, and as starting points for your own star decks.

| deck | mesh | star | cells | MeshBlocks | run |
|---|---|---|---|---|---|
| `inputs/rad_vet_clean/star3d_sp.athinput` | spherical-polar wedge: r 20-80 Rsun, theta pi/3..2pi/3, phi 0..pi/3 (periodic lateral walls) | blue supergiant envelope (the BSG production deck at 1/4 of its lateral resolution) | 256 x 64 x 64 = 1.0 M | 16 of 256 x 16 x 16 | 4000 cycles, 9.8 d |
| `inputs/rad_vet_clean/star3d_cs.athinput` | cubed sphere, full sphere: r 25-2530 Rsun (stretched; photosphere at 834 Rsun), 6 faces x 32 x 32 | red supergiant (25 Msun, MESA profile; the RSG production deck at 1/4 of its lateral resolution) | 548 x 32 x 32 x 6 = 3.4 M | 24 of 548 x 16 x 16 | 1200 cycles, 13.9 d |

Both decks start fresh from a 1-D profile. On the first cycle the code relaxes the
radiation to its steady state in the initial gas (`clean_init = steady`). The cs deck sets
`clean_init_tol = 2.0e-7`, see section 3. Both decks use tile layout mode with
16 x 16 lateral tiles (`layout_tile = 16`): any number of ranks that divides the MeshBlock
count gives the same answer, bit for bit.

## 1. What you need

1. An AthenaK binary built with `-D PROBLEM=he_star_vet` ([BUILD.md](BUILD.md)): a CPU build
   for the start-up check, and a CUDA (or HIP) build for the GPU run.
2. The input data in `inputs/rad_vet_clean/star3d_data/`. The decks read them by the relative
   path `star3d_data/...`, so link that directory into the directory you start the run from.

| file | size | used by | what |
|---|---|---|---|
| `bsg_ic.txt` | 0.9 MB | sp | 1-D initial state (r, rho, e, F_r, E; he_ic_cols = 5): the Ma et al. (2026) BSG profile, radiative-only outward integration to the photosphere, isothermal atmosphere, F_r = 0 above the photosphere |
| `bsg_profile.txt` | 0.3 MB | sp | the MESA profile (r, T, rho, M, kappa, F_r) that `he_ic_profile` interpolates (the BSG reference initial state) |
| `bsg_rosseland.txt` | 0.8 MB | sp | Rosseland-mean opacity incl. electron scattering, OPAL via MESA `kap` (GS98, X 0.70, Z 0.008), on log T x log rho |
| `bsg_planck.txt` | 0.8 MB | sp | Planck-mean absorption opacity, TOPS above log T 4.2 and Ferguson et al. (2005) below (blend 4.0-4.2), same grid and composition |
| `rsg_ic.txt` | 2.6 MB | cs | 1-D initial state of the 25 Msun RSG (MESA model, alpha_MLT 3), built for the code's EOS and opacities |
| `rsg_rosseland.txt` | 0.8 MB | cs | Rosseland mean (OPAL / Ferguson via MESA `kap`, X 0.654, Z 0.0198) |
| `rsg_planck.txt` | 0.8 MB | cs | Planck mean (TOPS / Ferguson, same composition) |

The two BSG profiles are thinned copies of the production files (28 MB and 49 MB). Rows
are kept every 3e9 cm and wherever a column changes by more than 0.002 dex. Linear
interpolation of the copies reproduces the full files to better than 3e-4. Each file's
header records its source. The cs deck builds its gas EOS table (H/He ionisation) at start-up
from the `<hydro>` `eos_*` keys, so it needs no EOS file. md5 sums:

```
ee9bf2ba1faef43372c268e34ccd0db6  bsg_ic.txt
8c9460946b4580ba4a5a8b7173b26220  bsg_profile.txt
c116f402987a6ce03e2d8f872e0ff94c  bsg_rosseland.txt
512e97c29c4e8785d4ebcda33e667516  bsg_planck.txt
125db6efdff0a72289042e6a50aae3cf  rsg_ic.txt
a9d3440e1aa5a4718e081b03714bb8fe  rsg_rosseland.txt
d8382173d238bf2a3d9fe1afb4129ce6  rsg_planck.txt
```

## 2. Run

Set up a run directory once (paths relative to the AthenaK checkout `$ATHENAK`):

```bash
mkdir -p star3d && cd star3d
ln -s $ATHENAK/inputs/rad_vet_clean/star3d_data .
cp $ATHENAK/inputs/rad_vet_clean/star3d_sp.athinput $ATHENAK/inputs/rad_vet_clean/star3d_cs.athinput .
```

### One CPU node: does it start? (a few cycles)

Use one MPI rank per MeshBlock (16 for sp, 24 for cs) and one OpenMP thread per rank:

```bash
export OMP_NUM_THREADS=1
mpirun -np 16 $ATHENAK/build_cpu/src/athena -i star3d_sp.athinput -d run_sp_cpu time/nlim=3 > run_sp_cpu.log 2>&1
mpirun -np 24 $ATHENAK/build_cpu/src/athena -i star3d_cs.athinput -d run_cs_cpu time/nlim=2 > run_cs_cpu.log 2>&1
```

Measured on one 2 x 56-core Intel Xeon Platinum 8480+ node (gcc 13, Open MPI 4.1, ranks pinned to cores): sp
3 cycles in 128 s; cs 2 cycles in 125 s (685 s when the node is shared with other jobs).
Most of that is the start-up (the profile, the ray tables and the steady radiation start).
A cycle then costs about 25 s (sp and cs). Fewer ranks also work (any divisor of the MeshBlock
count), only slower.

### 4 GPUs: the test run

Use one MPI rank per GPU. With Slurm:

```bash
srun -N 1 -n 4 --gpus-per-node=4 env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 \
  $ATHENAK/build_cuda/src/athena -i star3d_sp.athinput -d run_sp > run_sp.log 2>&1
srun -N 1 -n 4 --gpus-per-node=4 env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 \
  $ATHENAK/build_cuda/src/athena -i star3d_cs.athinput -d run_cs > run_cs.log 2>&1
```

With `mpirun`, use `mpirun -np 4 -x KOKKOS_MAP_DEVICE_ID_BY=mpi_rank -x OMP_NUM_THREADS=1 ...`.
For a run over several nodes of A100s with InfiniBand (e.g. MPCDF Raven), also set
`UCX_IB_GPU_DIRECT_RDMA=0`. One node does not need it. Over several H200 nodes (e.g. Caltech HPC), set
`UCX_TLS=^xpmem UCX_IB_GPU_DIRECT_RDMA=no UCX_RNDV_SCHEME=put_ppln`: without them the first message
between nodes fails with `ib_mlx5 ... Local protection error`. With `put_ppln`, two nodes were bitwise
from run to run.

Measured on one node of 4 A100-40GB (MPCDF Raven, `gcc/13 cuda/12.6 openmpi_gpu/5.0`):

| deck | wall time | per cycle | simulated time | output |
|---|---|---|---|---|
| sp | 13.9 min (832 s) | 0.21 s | 4000 cycles, 9.8 d (dt about 260 s) | 6 dumps of hydro_w + vet, 0.6 GB |
| cs | 13.3 min (797 s) | 0.65 s | 1200 cycles, 13.9 d (dt 1000 s) | 6 dumps of hydro_w + vet, 1.9 GB |

To make a run shorter or longer, change `time/nlim` (or give `time/nlim=N` on the command
line). Dumps are written every 1000 cycles (sp) or 300 cycles (cs), set by `output2/dcycle`
and `output3/dcycle`, and the history every cycle.

## 3. What to look at

Outputs (in the `-d` directory; prefix `bsg3d` for sp, `rsg` for cs):

- `PFX.user.hst`: the energy budget and luminosities, one row per cycle. The column names are
  in the header line `#  [1]=time  [2]=dt  [3]=L_bot ...`; read columns by name, because
  the two decks have different columns (the RSG adds `w_mlt` and `L_MLT`).
- `bin/PFX.hydro_w.NNNNN.bin` (rho, v, e) and `bin/PFX.vet.NNNNN.bin` (E, F, and the
  flux factors): the dumps.
- The run log (stdout): every `CLEAN WARNING` line and the end-of-run solver summary
  ([OUTPUTS.md](OUTPUTS.md)).

Checks:

1. The return code is 0, and the log ends with the end-of-run summary.
2. `grep -c 'CLEAN WARNING' run_sp.log` gives 0 (the same for cs).
3. The steady start converged (log line `clean_init steady: converged after N passes`).
   The line `clean_init steady pass 1: ... = ...` is a fingerprint of the binary, and a
   new machine or compiler should reproduce it to a few digits: sp `9.025e-09`, cs
   `1.218e-07` (the same on CPU and A100). On the RSG the steady-start residual stalls
   near 1.1e-7 of max E. The default tolerance 1e-8 is then never reached, and the start
   ends after 40 passes with a `CLEAN WARNING`, so the cs deck accepts 2e-7. Do not use a
   large pseudo time step (`clean_init_dtfac = 1e6`, as in the BSG deck) on the RSG. It
   does converge, but it moves the RSG surface layers so far from equilibrium that `L_top`
   falls to 0.28 `L_in` within a day and stays there (measured; with the default it stays
   near 1).
4. Luminosity. `L_in` is the luminosity fed in at the inner boundary, `L_top` the
   luminosity leaving through the outer boundary, `L_bot` the radiative luminosity through
   the innermost faces (erg/s; on the sp wedge these are integrated over the wedge only).
   Expected at the end of the 4-GPU runs:

| deck | L_in | L_top at the end | L_top / L_in at the end | L_top / L_in over the second half |
|---|---|---|---|---|
| sp | 4.651e37 (wedge) | 4.83e37 | 1.039 | 0.958 - 1.066 |
| cs | 6.837e38 | 6.74e38 | 0.986 | 0.986 - 0.997 |

What the curves look like (the notebook plots them):

- sp: at t = 0, `L_top` is 1600 `L_in`, because the radiation of the initial isothermal
  atmosphere (E = a T^4, F = 0 above the photosphere) streams out. It falls to about
  1.4 `L_in` within 10 cycles and to within 5 % of `L_in` after about 0.3 d. After that it
  oscillates around `L_in` by up to 10 % (0.96 - 1.07 over the second half). `L_bot`
  swings between -1.4 and +2.2 `L_in` over the first 100 cycles while the inner wall
  adjusts, then stays within 10 % of `L_in` (0.94 - 1.05 over the second half).
- cs: `L_top / L_in` starts at 0.987, peaks at 1.10 after about 0.1 d and relaxes back to
  0.986 - 0.997 over the second week. `L_bot / L_in` stays within 3 % of 1, and within
  2e-4 of 1 after the first day.

A run that departs from these values by much more than the oscillation amplitude, or
that prints `CLEAN WARNING` lines, is worth a closer look ([OUTPUTS.md](OUTPUTS.md)).

## 4. Notebook

`docs/vet_clean/notebooks/star3d.ipynb` (numpy and matplotlib only) reads the two runs and
plots:

- `L_top`, `L_bot`, `L_in` and `L_top / L_in` against time, from `user.hst` by column name;
- the density contrast `rho / <rho>_shell - 1` and `log10 E` of the last dump in a
  meridional and an equatorial plane (sp: the middle phi column and theta row; cs: the planes
  y = 0 and z = 0 through the cube faces);
- shell-averaged rho and E of the first and the last dump.

Set the two run directories in the first cell (`RUNS`), or set the environment variables
`STAR3D_SP` and `STAR3D_CS`, then run all cells:

```bash
STAR3D_SP=run_sp STAR3D_CS=run_cs jupyter nbconvert --to notebook --execute \
  $ATHENAK/docs/vet_clean/notebooks/star3d.ipynb --output star3d_run.ipynb --output-dir .
```

The notebook in the repository was executed on the 4-A100 runs above.

## 5. From the test to a production run

The production templates (`star_bsg.athinput` and the others, [RUN.md](RUN.md)) differ from
these decks only in resolution, MeshBlock size (`layout_tile` must equal the lateral
MeshBlock size), run length and outputs. Both sp decks use the Gauss-Legendre 7 x 18 ray
set (`vet_gd_gl_nmu = 7`, `vet_gd_gl_nphi = 18`, `vet_gd_gl_stagger = true`). On a
phi-periodic wedge the width must be a multiple of the ray set's azimuthal period
(360 deg / nphi = 20 deg here, the wedge is 60 deg; LQ10: 90 deg; the default icosahedral set:
180 deg), else the run stops with a FATAL. A narrow wedge that cannot meet this (e.g. 0.064 rad)
names `vet_gd_wall_sym = true` with `vet_gd_wall_interp = true` ([KEYS.md](KEYS.md)).
