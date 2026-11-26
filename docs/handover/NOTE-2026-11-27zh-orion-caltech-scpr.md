# NOTE 2026-11-27zh orion -> Caltech: AthenaK PR #798 (short-characteristics RT) vs our rays, H200 timing + sphere accuracy

The USER asks how our short-characteristics rays compare with the open upstream AthenaK pull
request IAS-Astrophysics/athenak #798, "Gray LTE short-characteristics radiative transfer
(nr_radiation)", by apbailey, tip **c5ac58778**. The coordinator moved this study off Orion CPU,
so the timing runs on one idle H200 node. GPU only, no CPU runs needed. Matched problem only.

## Build (in a job, not on the login node)
- PR: `git fetch https://github.com/IAS-Astrophysics/athenak.git pull/798/head:upstream-pr798`
  (check `git rev-parse --short=9 upstream-pr798` = c5ac58778), in its own worktree. Its kokkos
  submodule is 4.7.02; your Kokkos 4.6.2 of the rel6 build is fine too (say which one you used).
  Use **CUDA 12.9 + g++13** (not 12.2, see the miscompile note), Release, `-D Kokkos_ENABLE_CUDA=On
  -D Kokkos_ARCH_HOPPER90=On -D Athena_ENABLE_MPI=ON -D PROBLEM=built_in_pgens`, the same flags
  as your rel6 suite binary otherwise. Post the md5.
- Ours: your rel6 suite binary (md5 8cc98fb3, src = release aa2646b92), the one 27ze/27zg use.
  No build needed.

## Files (this branch, `docs/handover/scpr-1127/`)
- `pr_sph.athinput`: PR deck, homogeneous emitting sphere (R 0.6, sigma_a 5, S 1) in vacuum on
  [-1,1]^3, 64^3, one MeshBlock, Carlson nmu = 6 (168 rays), writes `sc_rad_J` (bin).
- `our_sph.athinput`: ours, the same sphere (hot_blob with a sharp edge, chi = 5, a T^4 = 1,
  background rho 1e-12, gas frozen, vacuum faces), Cartesian rays `vet_sc_set = gd`, GL 8 x 20
  unstaggered (160 rays), writes `vet_vsc` (the rays' J is `vet_vj`).
- `sph_err.py`: J against exact transport over the interior cells r < 0.36
  (`BINPY=<athenak>/vis/python python3 sph_err.py FILE VAR [gl 8 20]`).
- `kp_ktime.cpp`: the SCSPEED kernel timer (`g++ -O2 -std=c++17 -shared -fPIC -o libkp_ktime.so
  kp_ktime.cpp`; you built it for 27ze, so reuse that one).

Let `U = problem/radius=10.0 problem/sigma_a=16.0 problem/r_mask=0.0 problem/tol=1e30
output1/dt=1e30 time/nlim=10 nr_radiation/iter_min=1 nr_radiation/iter_max=1` (the whole box is
then a uniform absorber, cell dtau 0.25 at 128^3; one sweep per cycle) and `N128 = mesh/nx1=128
mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128`.

## Runs (one rank per GPU, your usual binding; GPU numbers first)
0. Smoke, 1 GPU: `pr_sph.athinput` as is (`time/nlim=1`). If it FATALs, post that and stop.
1. **G1 timing, 1 GPU, 128^3, one MeshBlock, kernel timer on**
   (`KOKKOS_TOOLS_LIBS=.../libkp_ktime.so KP_OUT=<rundir>/kp`), order P1 P2 P1 P2:
   - P1 (PR, plane kernel, the default): `-i pr_sph.athinput $N128 $U`
   - P2 (PR, tiled kernel, its GPU option): the same plus `nr_radiation/sweep_kernel=tiled nr_radiation/tile_size=16`
   - Ours: the 27zg G1 A (Cartesian) and B (sp) runs. If they ran on this node in the same
     session, reuse them. Otherwise run G1 A and B once each here with the timer, as in 27ze.
2. **G4 timing, 4 GPUs, 128^3, 4 MeshBlocks split laterally** (`meshblock/nx1=128 meshblock/nx2=64
   meshblock/nx3=64`): P1 with `$U` and the timer, once. Ours = 27zg G4 A/B.
3. **G4 iterations (PR only), 4 GPUs**: `-i pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128
   meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 time/nlim=5 output1/dt=1e30
   nr_radiation/iter_tol=1e-6`. The PR iterates block sweeps + ghost exchange until J converges;
   post every `Iteration stops at niter: ...` line and the s/cycle.
4. **Accuracy, 1 GPU, 64^3**: `-i pr_sph.athinput` (post the `SC sphere:` line) and
   `-i our_sph.athinput` with our binary. Then
   `sph_err.py <PR last bin> sc_rad_J` and `sph_err.py <our last bin> vet_vj gl 8 20`. Post the
   printed lines. If time allows, repeat both at 32^3 (`mesh/nx1..3=32 meshblock/nx1..3=32`).

Per run please post: rc, the s/cycle over cycles 2..end, `cpu time used`. For runs with the timer,
the top 25 lines of `kp.0` (PR kernel names `sc_sweep_plane`, `sc_sweep_tiled`, `sc_carry_j`,
`sc_residual_j`, the moment kernels and the bvals pack/unpack). For ours, the `vet_sc: formal
solutions=...` / `vet_gd: ... fenced seconds ...` lines as in 27ze. Also any FATAL or WARNING
line. Delete the rst/bin files after the analysis. Post a short NOTE back (caltech -> orion,
"scpr"). Do not touch vet-clean-release or any production.

-- orion (SCPR, for the coordinator)
