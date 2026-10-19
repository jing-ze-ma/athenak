# TASK 2026-10-19 orion -> raven: Arepo2 + IDORT on A100 -- feasibility build and first measurements

**Test only, report back by NOTE** (docs/handover/NOTE-2026-10-2x-raven-arepo2-idort-gpu.md on bsg-files-1009).
No production runs. Do not push to the Arepo repo (Bitbucket); local build branches only.

## Why
Arepo2 IDORT (the IDORT radiation module in Arepo2) is validated on CPU for the RSG g25 star (Orion, 560 ranks). Per g25 top step
(late in the run) ~42.6 s: radiation ~30 s (IDORT, CPU only, pure MPI, **no OpenMP, no GPU code**), mesh + gravity + hydro
~12.7 s. On this branch only the mesh (VORONOI_CELL_GPU) and gravity (FORCETREE_NODELOCAL_GPU) have GPU code.
Questions this task answers:
1. Does Arepo2-IDORT-dev build and run with the developer's GPU recipe plus IDORT?
2. **How many MPI ranks per node do CPU work (IDORT) in GPU mode?** IDORT is one rank per core. If GPU mode leaves only a
   few ranks per node for CPU work, IDORT slows down by that factor and the GPU run loses overall. This is the deciding number.
3. Are the results the same as a CPU build (within noise; bitwise is not expected)?
4. The time split (mesh / gravity / hydro / IDORT) GPU vs CPU on the same node count.

## Code
- Repo: bitbucket `volkerspringel/arepo` (ssh host `ssh.bitbucket.org`), branch **Arepo2-IDORT-dev**, head **b0ed6a23db**
  (2026-10-19). IDORT: `src/modules/idort`; tests `examples/idort/` (README.md there: `./test.py idort/<dir>` from
  `examples/`, each dir has Config.sh, param.txt, create.py, check.py).
- Base Arepo2_WD is under heavy development by its developer, only tested on NVIDIA A100; feedback welcome.
- IDORT_DSA needs HYPRE (Orion: own static build of HYPRE 2.32.0, MPI, 32-bit ints). Pass
  `HYPRE_INCL=... HYPRE_LIB="-L... -lHYPRE"` to make; modules-lib.make blanks the HYPRE variables unless the Config switch is
  listed there (add IDORT_DSA if it is not). If HYPRE is a problem, drop IDORT_DSA for the first build and say so.

## GPU recipe (from the Arepo2_WD developer, relayed by the user 2026-10-12)
Config (add):
    VORONOI_CELL
    FORCETREE_NODELOCAL
    GPU
    VORONOI_CELL_GPU
    FORCETREE_NODELOCAL_GPU
    GPU_CUDA                     # (GPU_HIP for AMD, GPU_SYCL for Intel)
    TIMESTEPS_WAKEUP
    SHAREDNODEMEMORY
    HIERARCHICAL_DOMAIN_DECOMPOSITION
Config (remove): TREE_BASED_TIMESTEPS.
Params:
    NumIslandsPerNode    4
    CPUsUsedPerGPU       1
    MeshCandidateSource  1
    MeshCellRetirement   1
    GpuMaxMemSize        ~80 % of the GPU memory
    MaxMemSizeNode       similar or slightly smaller
    MaxMemSize           for the rest
CPU-only builds: drop the *_GPU / GPU* options; VORONOI_CELL + FORCETREE_NODELOCAL still need SHAREDNODEMEMORY and
TIMESTEPS_WAKEUP (or FORCE_EQUAL_TIMESTEPS).
Orion gotchas that may also apply: VORONOI_CELL runs > ~250k cells need `export AREPO_NGBNL_GRAPH_POOL=160`; parallel `make`
from one clone races on the git submodule lock (pre-touch `<builddir>/git_submodules.check`); batch mpirun
`--bind-to core --map-by ppr:N:node` is ~15 % faster than `--bind-to none`; link may need `HWLOC_LIB="-lhwloc -lopen-pal"`.

**Differences from the g25 production Config** (examples/idort/G25_setup/Config.sh) that the recipe forces -- report what
breaks, do not fix the physics:
- g25 uses the classic gravity tree (SELFGRAVITY, HIERARCHICAL_GRAVITY, EXACT_GRAVITY_FOR_PARTICLE_TYPE=1, EXACT_GRAVITY_REACTION,
  GRAVITY_PM_ENERGY_CONSERVING, ALLOW_DIRECT_SUMMATION, ...), not FORCETREE_NODELOCAL. Whether FORCETREE_NODELOCAL(_GPU) accepts
  these point-mass options is unknown.
- g25 uses TREE_BASED_TIMESTEPS (the IDORT radiation signal speed enters through it). Removing it changes the time steps.

## Steps
A. **Plumbing, small** (no input files needed): `examples/idort/I4_star_lite` (3-D star, VORONOI_CELL, EXTERNALGRAVITY,
   IDORT + DSA, 8 ranks; `check.py` checks the budget and the star). Build (i) its Config as is (CPU) and (ii) + the GPU recipe.
   Run both on 1 Raven node; check.py must pass for both. Report: which recipe options compile with IDORT, run.log
   "IDORT: timing" lines, the rank layout (MPI ranks per node, ranks driving a GPU, ranks doing CPU work; how
   NumIslandsPerNode / CPUsUsedPerGPU map to this), and cpu.txt totals for both.
B. **Rank layout scan** on I4_star_lite or A2_sphere_3d (235k cells): GPU build with 4, 8, 18, 36, 72 ranks per node (if the
   code allows more ranks than GPUs). Time per step split (cpu.txt: mesh, gravity, hydro, IDORT). This answers question 2.
C. **g25 reduced, only if A works**: IC `g25_red.hdf5` (1.19M cells, 143 MB, md5 999c75f387dfc36d3db1040b3836d2c3) is on
   Orion at `/orion/ptmp/jinma/idort_dev/runs_g25/ics/g25_red.hdf5`; **the user copies it to Raven** (Orion cannot ssh to
   Raven). Config: examples/idort/G25_setup/Config.sh (+ GPU recipe, with the gravity changes above), param.txt there
   (InitCondFile g25_red; smoke-test TimeMax). Same run on CPU-only on the same node count. Compare: the first full sync's
   IDORT solver line (idort_solver.txt), energy budget (idort_budget.txt, res_E), L(r) in idort_shell_flux.txt, wall per
   step split. If the gravity options do not combine with FORCETREE_NODELOCAL, run VORONOI_CELL_GPU only (classic gravity) and say so.

## Report (NOTE)
- Build outcome per option set (errors verbatim, short).
- Rank layout in GPU mode and the per-node CPU rank count left for IDORT.
- Tables: per-step wall split CPU vs GPU (A, B, C as done); check.py results; first-full-sync IDORT line CPU vs GPU.
- Anything to send to the Arepo2_WD developer (bugs, missing options), as a separate list.
Delete restart dumps and large snapshots after each run.
