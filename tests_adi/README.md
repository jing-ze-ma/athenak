# Tests of the implicit transverse radiative operator on the cubed sphere

Branch `cs-implicit-transverse` (off `rt-integration` 7f9f68a7).  Design note:
`docs/dev/cs_implicit_transverse.md`.  Regression test:
`tst/test_suite/rad/test_rad_cs_implicit_ang_cpu.py`.

Everything below is reproducible from this directory.  `measure_cs.py` is the analysis
(the same fit `test_rad_cs_raddiff_cpu.py` uses), `amp_hist.py` the amplitude history,
`cmpbin.py` the bitwise comparison of two dump directories (it compares the DATA, not
the files: the binary dump embeds the effective input, and this branch adds parameters
to it, so `cmp` on the file always differs), `ladder.sh` the convergence ladder,
`build_gpu.sh` a GPU build, `run_gpu_cs.sh` / `run_gpu_bitwise.sh` the apudev jobs.

## Builds

```bash
# CPU, for the cubed-sphere problem and for the Cartesian bitwise arm
cmake -S .. -B ../build_cpu_chk -D PROBLEM=cs_test -D CMAKE_BUILD_TYPE=Release
make -C ../build_cpu_chk -j
cmake -S .. -B ../build_cpu_def -D CMAKE_BUILD_TYPE=Release
make -C ../build_cpu_def -j

# GPU (MI300A APU): modules gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0, hipcc,
# Kokkos_ENABLE_HIP, Kokkos_ARCH_AMD_GFX942_APU, Athena_ENABLE_MPI=ON
./build_gpu.sh cs_test          # -> ../build_gpu_cs_test
./build_gpu.sh default          # -> ../build_gpu_default
./build_gpu.sh default <head-worktree>   # the unmodified 7f9f68a7 binary
```

The reference ("head") binaries are built from a detached worktree at 7f9f68a7.

## (a) Cartesian bitwise regression

```bash
# head and new, 2D Gaussian transverse diffusion (pgen rad_diff2d), both solvers
for v in head new; do for s in sts adi; do
  $BIN -i rad_transverse_gauss_bin.athinput time/nlim=30 \
       hydro/rad_ang_solver=$s hydro/rad_adi_scheme=lod2
done; done
python3 cmpbin.py bb_head_sts/bin bb_new_sts/bin
python3 cmpbin.py bb_head_adi/bin bb_new_adi/bin
```

| arm | result |
| --- | --- |
| `rad_ang_solver = sts`, 14 `hydro_u` dumps | **BITWISE IDENTICAL** |
| `rad_ang_solver = adi` (`lod2`), 14 dumps | **BITWISE IDENTICAL** |
| `lap2d-errs.dat` | byte-identical |
| stdout | identical apart from wall-clock lines |

The operator's own conservation report on the Cartesian problem is unchanged at
`|sum V de|/sum V|de| = 1e-15 .. 1e-17`, i.e. round-off.

### the PRODUCTION Cartesian configuration: the He-star FeCZ box

The Gaussian test above is a thin gate -- ideal gas, one MeshBlock, no taper, no radial
solve.  The real one is `box_w8`, the running He-star production: tabulated general EOS,
`rad_implicit_x1`, the two-stream/tau blend and the transverse taper, the mode-3 implicit
RT column, `rad_ang_solver = adi` with `rad_adi_scheme = lod2`, and **2 x 4 MeshBlocks per
ADI line** on 2 ranks -- so the partitioned line solve, its reduced interface system and
its MPI ring gather are all in the comparison.  `he_box_w8_smoke.athinput` is that input
with the `<output*>` blocks replaced by a single `hydro_u` dump every step.

```bash
sbatch run_gpu_hebox.sh     # head vs new, nlim = 30, 2 x MI300A
python3 cmpbin.py hebox_head/bin hebox_new/bin
```

**BITWISE IDENTICAL over 32 dumps.**

It was not, at first, and the failure is worth recording because it is the one real trap
in this change.  The first version of the ADI tridiagonal row read

```c++
const Real vi = vcell(m,k,j,i);          // exactly 1.0 on a Cartesian mesh
Real dj = 1.0 + thtau*ai*(cl + cr)/vi;
```

and a second version hoisted the Cartesian and curvilinear forms into an `if (curv)
... else ...`.  Both are algebraically no-ops off a cubed sphere, and both moved the last
bit of the solution: dividing by an exact 1.0, or putting the expression in a branch,
stops `hipcc` contracting `1.0 + thtau*ai*(cl+cr)` into a single FMA, so `dj` changes by
one ULP.  It showed up as `mom2` differing by 3.55e-15 at dump 2 and growing
exponentially, while `max abs(de)` and `max z_i` stayed bit-identical for ten cycles --
i.e. the operator's largest increment was unchanged and only its tiniest cells moved.

The fix is to leave the three Cartesian statements exactly where they were and append the
curvilinear form as an overwrite.  The diagnosis was a five-step bisection on the GPU:

| binary | He box vs head |
| --- | --- |
| head rebuilt with a no-op extra View capture in the ADI RHS kernel | identical (so it is *not* register pressure / reduction trees) |
| mine, `BuildAngularCoeffs` edits reverted | differs |
| mine, `conduction_transverse.cpp` reverted to head | identical |
| mine, all ADI-section hunks reverted | identical |
| mine, RHS + factorisation hunks reverted | identical |
| mine, **factorisation hunk only** reverted | identical -> that hunk is the culprit |

Every other hunk of this branch -- the geometry locals, the boundary predicates, the
Gershgorin restructure, the RHS hunk, the cross-term iteration wrapper, the seam
sub-step -- is bitwise clean on this input.

## (b) Cubed-sphere diffusion

`cs_test` `iprob = 15`, `inputs/tests/cubed_sphere_raddiff.athinput`: uniform density at
rest, `T = T0 (1 + amp f)` with `f` a harmonic of degree `l` and no radial dependence, so
the radial fluxes vanish and the angular operator is the only thing that moves energy.
`amp` below is the fitted ratio of the discrete angular source to the exact
`-l(l+1) kappa (T-T0)/r^2`, i.e. **1 means the harmonic decays at exactly
`exp(-l(l+1) D t/r^2)`** and `|amp - 1|` is the rate error.

The ladder scales `dt` as `dx^2` at a fixed final time (`cfl = 0.03 / 0.0075 / 0.001875`,
`nlim = 20 / 80 / 320`), so the operator-splitting error is the same at every resolution
and what is compared is the space discretisation.

```bash
./ladder.sh
```

### l = 2, rate error |amp - 1|

| n per panel edge | explicit (reference) | `sts` (RKL1) | `adi` (`lod2`) |
| --- | --- | --- | --- |
| 16 | 0.362 % | 0.497 % | 0.511 % |
| 32 | 0.034 % | **0.328 %** | **0.335 %** |
| 64 | 0.032 % | 0.285 % | 0.288 % |

Both implicit solvers are inside the 2 % target at 32 cells per panel edge.  The rate
error **floors** at ~0.29 %, and that floor is *not* spatial: this ladder holds the
per-stage stiffness `z = tau lambda` fixed by construction, so the first-order splitting
error of running the operator after the RK update is the same at every resolution.
Measured separately with `z -> 0` (`nlim = 20`, so the field barely decays and what is
read off is the operator itself), the SPATIAL error converges better than second order:

| n | `sts`, z ~ 0.05 |
| --- | --- |
| 16 | 0.497 % |
| 32 | 0.064 %  (7.8x) |

(64 cannot be measured this way: the binary dump is `float32` and the per-step energy
change is below its resolution.)

### residual shape, seam signature and conservation, l = 2

| n | solver | L1 | max residual, panel-edge cells | ... panel interior | edge/interior |
| --- | --- | --- | --- | --- | --- |
| 16 | explicit | 0.0251 | 0.0475 | 0.0300 | 1.58 |
| 16 | `sts` | 0.0244 | 0.0376 | 0.0305 | 1.23 |
| 16 | `adi` | 0.0246 | 0.0374 | 0.0316 | 1.18 |
| 32 | `sts` | 0.0245 | 0.0380 | 0.0341 | 1.12 |
| 32 | `adi` | 0.0246 | 0.0395 | 0.0343 | 1.15 |
| 64 | `sts` | 0.0244 | 0.0348 | 0.0354 | 0.98 |
| 64 | `adi` | 0.0245 | 0.0357 | 0.0356 | 1.00 |

The panel seam is **not** a special place for either solver: the edge residual is within
15 % of the interior one at 32 and indistinguishable at 64.  (The explicit reference is
*worse* at the edge, 1.58x at n = 16.)

Shell-integrated conservation, the operator's own
`|sum_i V_i de_i| / sum_i V_i |de_i|` (worst over the run):

| n | `sts` | `adi` |
| --- | --- | --- |
| 16 | 1.83e-5 | 1.87e-5 |
| 32 | 2.14e-6 | 1.32e-6 |
| 64 | 7.83e-7 | 1.23e-6 |

against 1e-15 .. 1e-17 on a Cartesian mesh.  The difference is the **along-seam
resample**: the two panels meeting at a seam do not share cell centres, so each forms
the seam flux from its own resampled halo and the two are not bitwise opposite.  It
converges with resolution (24x from 16 to 64 for `sts`), and it is the error the
EXPLICIT cubed-sphere operator already makes at a seam -- nothing in this operator adds
to it.

### l = 6 (`problem/lharm = 6`, the zonal P_6, Laplacian_S f = -42 f/r^2)

| n | explicit | `sts` | `adi` |
| --- | --- | --- | --- |
| 16 | 2.65 % | 2.96 % | 3.01 % |
| 32 | 0.60 % | **1.25 %** | **1.28 %** |

Conservation 3.7e-4 (16) and 3.1e-4 (32): the seam resample error is larger because the
field varies four times faster along the seam, and at these resolutions it does not yet
converge.

### the seam sub-step weight (`hydro/rad_adi_seam_w`)

`adi`, n = 32, `cfl = 0.03`, `nlim = 20`; max residual at the panel edge / in the
interior:

| `rad_adi_seam_w` | edge | interior | edge/interior |
| --- | --- | --- | --- |
| 0 (seam flux at the pre-sweep state) | 0.1596 | 0.0598 | 2.67 |
| **0.5 (trapezoidal, the default)** | **0.0692** | 0.0371 | **1.87** |
| 1 (post-sweep state) | 0.1086 | 0.0341 | 3.18 |

The outliers are the **cube-vertex** cells, the only cells with two seam faces.

### the cross-term outer iteration (`hydro/rad_adi_cross_iter`)

`adi`, n = 32, `cfl = 0.0075`, `nlim = 80` (the ladder's 32 arm):

| `rad_adi_cross_iter` | amp | L1 | `|sum V de|/sum V|de|` |
| --- | --- | --- | --- |
| 1 | 0.99665 | 0.02462 | 1.3208e-6 |
| 2 | 0.99666 | 0.02462 | 1.3196e-6 |
| 3 | 0.99666 | 0.02462 | 1.3196e-6 |

It converges in two iterations and changes nothing here, which is the expected result:
the lag error of the explicit cross term scales with the per-stage stiffness `z`, and
this ladder runs at `z ~ 0.05`.  The switch exists for the production regime.

## (c) Stability

`hydro/rad_implicit_x1 = true` as well, or the radial conduction step shrinks `dt` with
the conductivity and the transverse diffusion number never gets large.

```bash
# z ~ 2e3
$BIN -i ../inputs/tests/cubed_sphere_raddiff.athinput mesh/nx2=32 mesh/nx3=32 \
     meshblock/nx2=32 meshblock/nx3=32 time/nlim=40 hydro/rad_kappa_fac=3.0e-3 \
     hydro/rad_implicit_x1=true hydro/rad_implicit_ang=true hydro/rad_ang_solver=adi
```

| `rad_kappa_fac` | max `z_i` | solver | result |
| --- | --- | --- | --- |
| 3e-3 | 2.24e3 | `adi` | bounded, **monotone**, no sign change, `A/A0 = 2.08e-3` after 40 cycles |
| 3e-3 | 2.24e3 | `sts` (50 substages) | bounded, `A/A0 = -4.97e-5` (one sign change at the 1e-4 level, i.e. fully relaxed) |
| 1e-5 | 6.71e5 | `adi` | bounded, `max abs(de)` falls monotonically 1.03e-2 -> 5.27e-3 -> 9.91e-4 -> 5.85e-4 over the first cycles |
| 1e-5 | 6.71e5 | `sts`, `rad_ang_maxit = 2000` | 859 substages, relaxes to `A/A0 = 3.5e-6` |
| 1e-5 | 6.71e5 | `sts`, `rad_ang_maxit = 200` (default) | **CLAMPED and blows up**: `max abs(de) = 2.3e303`, `dt -> 0` |

The last line is a pre-existing property of the RKL1 clamp, not of the cubed sphere: at
`z = 6.7e5` RKL1 needs 859 substages and the default ceiling is 200, at which point the
super-step is not covered and the scheme is an over-long explicit one.  The ADI takes
ONE step however stiff the row is, which is the reason the production He/B-star boxes
use it, and it is the solver the cubed sphere should use at production stiffness.

## GPU (2 x MI300A APU, 2 MPI ranks, `apudev`)

```bash
sbatch run_gpu_cs.sh          # the cubed-sphere ladder and the stability arms
sbatch run_gpu_bitwise.sh <head-gpu-binary>   # the Cartesian bitwise arm
```

Two ranks means the panel seams are also RANK boundaries, so the module's own one-
variable exchange and (for `adi`) the seam sub-step run over MPI.  Every number the GPU
run produces agrees with the serial CPU run **to all printed digits**:

| arm | GPU amp | CPU amp | GPU `|sum V de|/sum V|de|` | CPU |
| --- | --- | --- | --- | --- |
| `sts` n=16 | 0.99503 | 0.99503 | 1.831665e-5 | 1.831665e-5 |
| `adi` n=16 | 0.99489 | 0.99489 | 1.868048e-5 | 1.868048e-5 |
| `sts` n=32 | 0.99672 | 0.99672 | 2.139308e-6 | 2.139308e-6 |
| `adi` n=32 | 0.99665 | 0.99665 | 1.320816e-6 | 1.320816e-6 |

Stability arms on the GPU:

| `rad_kappa_fac` | max `z_i` | solver | `A/A0` after 40 cycles | monotone | sign changes |
| --- | --- | --- | --- | --- | --- |
| 3e-3 | 2.24e3 | `adi` | 2.08e-3 | yes | 0 |
| 3e-3 | 2.24e3 | `sts` (50 substages) | -4.97e-5 | no | 1 |
| 1e-5 | 6.71e5 | `adi` | 1.86e-3 | yes | 0 |
| 1e-5 | 6.71e5 | `sts` (859 substages) | 3.44e-6 | no | 2 |

## What is refused

* spherical polar (unchanged fatal);
* SMR/AMR (unchanged fatal);
* on the cubed sphere: `rad_sts_all`, `rad_sts_split`, `rad_tr_halo_faces_only`, and for
  `rad_ang_solver = adi` more than one MeshBlock per panel in a transverse direction
  (a tridiagonal line cannot cross a panel seam; `sts` has no such restriction).
