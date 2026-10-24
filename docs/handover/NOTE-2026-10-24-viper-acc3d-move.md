# NOTE 2026-10-24 viper -> caltech: 3-D accretor run12 move bundle (answer to TASK-2026-10-24-caltech-acc3d-move.md)

## Status: FRESH START
run12 (deck acc3d_prod_D_B1v6_extrap2) has NOT started anywhere. The viper chain (12327914, 12327916-18, apu 4x2) and the Raven
chain (31291344..) are both still PENDING, and /viper/ptmp2/jinma/run12_viper_1022/run12 is empty. So the bundle holds no restart
and no hst: Caltech starts from t=0 with the static phase.

## Bundle
`/viper/ptmp2/jinma/acc3d_run12_move/` (README.txt inside, MD5SUMS covers every file)

| file | md5 |
| --- | --- |
| acc3d_prod_D_B1v6_extrap2_caltech_abs.athinput (ready to use, 3 paths -> /resnick/groups/carnegie_poc/jingze/acc3d_run12_move/in/) | a93b8f057521fc4be4a6c862817b7afd |
| acc3d_prod_D_B1v6_extrap2_caltech.athinput (same, relative `in/...` placeholders) | 31a122b746543302ce5ed85f2932f90b |
| deck_original/acc3d_prod_D_B1v6_extrap2_viper.athinput (as run on viper) | 2713294183c961d4ed34148b8b571bfc |
| deck_original/acc3d_prod_D_B1v6_extrap2_raven.athinput (as run on Raven) | c333a04c749ae2286df645f6bec0457f |
| in/column_mesakap_z0.02_p381.txt (`<problem>/env_ic_file`) | 5b813886cf9130a34546b201f6ba3e85 |
| in/rosseland_mesakap_gs98_x0.70_z0.02.txt (`<problem>/env_opac_table`) | d81b7a468b0d9470a09fd43ef03273d5 |
| in/planck_ext2_gs98_x0.70_z0.02.txt (`<problem>/env_planck_table`) | 99d668f69b16a502f50f1ec555fb0cd1 |
| scripts/link_viper.sbatch, scripts/prod_chain_raven.sbatch, scripts/run_raven.cfg | see MD5SUMS |
| scripts/CMakeCache_{viper_mi300a,raven_a100}_excerpt.txt | see MD5SUMS |
| ref/ : smoke logs + run.log + hst (viper 12355796, Raven 31453560) | see MD5SUMS |

The only file inputs are those 3 `<problem>` keys. The viper and Raven decks differ only in them. The general EOS is the in-code
Saha table, so no EOS file. The relative deck resolves `in/` against the process cwd, and both link scripts cd into
PD/static or PD/stream before srun. So use the `_abs` deck, or point env_ic_file / env_opac_table / env_planck_table at
/resnick/groups/carnegie_poc/jingze/acc3d_run12_move/in/.

## rsync (USER)
Direct, from a machine that reaches both:

    rsync -av viper:/viper/ptmp2/jinma/acc3d_run12_move/ <caltech-host>:/resnick/groups/carnegie_poc/jingze/acc3d_run12_move/

Two hops via the laptop:

    rsync -av viper:/viper/ptmp2/jinma/acc3d_run12_move/ ~/acc3d_run12_move/
    rsync -av ~/acc3d_run12_move/ <caltech-host>:/resnick/groups/carnegie_poc/jingze/acc3d_run12_move/

Then on Caltech: `cd /resnick/groups/carnegie_poc/jingze/acc3d_run12_move && md5sum -c MD5SUMS`. The bundle is about 6 MB.

## Phases (as in scripts/link_viper.sbatch and prod_chain_raven.sbatch)
- static, from t=0: `athena -i DECK problem/stream=false time/tlim=0.05 -d PD/static -t WALL`
- stream: `athena -r <newest PD/static/rst/ryper.*.rst> -i DECK -d PD/stream -t WALL`. Later links restart from the newest
  PD/stream rst, with no extra keys.
- XK (extra keys) is empty. WT = 85200 s of simulation wall per 24 h link.
- stop_check stops the chain on `NON-CONVERGED=[1-9]`, `FATAL`, `Picard solve DIVERGED`, or nan in the last 50 log lines.
  Raven keeps only the 3 newest stream rst (3.6 GB each).

## Binary / configure
Commit: fin51-r-fix **958cac780** (958cac78050834d414b6765304f57616cac7c8cd, on fork/fin51-r-fix and fork/csfix-1023).
- viper MI300A binary athena_apu_ry3d_958cac78_r71, md5 a282a9be93358724c2008a83c61f43bc:
  hipcc (rocm 7.1), `-D PROBLEM=ry_per_accretor -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off
  -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On`, Release.
- Raven A100 binary athena_a100_ry_per_accretor_958cac78_nofma, md5 1d93f5498af8c8e1b05dfaba8b414d0e:
  gcc/13 cuda/12.6, `-D PROBLEM=ry_per_accretor -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off
  -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_AMPERE80=On -D Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF
  -D CMAKE_CXX_COMPILER=<src>/kokkos/bin/nvcc_wrapper`, Release (Kokkos_ENABLE_CUDA_LAMBDA=ON, set by default).
- **Caltech H200 equivalent** (= your _nofma recipe):

      cmake -B build_acc -D CMAKE_BUILD_TYPE=Release -D PROBLEM=ry_per_accretor -D Athena_ENABLE_MPI=ON \
        -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF \
        -D CMAKE_CXX_FLAGS=-ffp-contract=off -D CMAKE_CXX_COMPILER=$PWD/kokkos/bin/nvcc_wrapper

## Rank layout
mesh 420x58x512, meshblock 420x58x64 -> **8 MeshBlocks -> 8 ranks, 1 per H200** (2 nodes x 4). This is the same as Raven
(2 nodes x 4 A100). The viper twin runs 4 nodes x 2 MI300A.

## Reference numbers (existing smokes, binary 958cac780, this deck)
The hst output cadence is coarse, so each smoke hst has a single row (static: the t=0 row, stream: the t=2.124e-4 row).

| | viper MI300A 12355796 (1 node, 2 ranks, 4 blocks/GPU) | Raven A100 31453560 (1 node, 4 ranks, 2 blocks/GPU) |
| --- | --- | --- |
| static t=0 hst: dt, mass, tot-E | 1.416096577201122e-05, 1.142252572723357e+06, -6.402933515344138e+11 | 1.416096577200850e-05, 1.142252572723363e+06, -6.402933515344166e+11 |
| static after 10 cycles (nlim 10) | time=1.416097e-04, Picard mean 3.379 (29 solves, max 4), NC 0 | same: time 1.416097e-04, Picard 3.379 / max 4, NC 0 |
| static speed | 4.57e6 zone-cycles/cpu_s | 7.77e6 zone-cycles/cpu_s; 0.83 s/cycle (cycles 1-9) |
| stream 20 cycles from smoke2 static rst ryper.00001 | hst t=2.124e-4: dt 1.415731279868494e-05, mass 1.142253008500556e+06, tot-E -6.402935117791125e+11; Picard mean 27.35, max 88, NC 0 | dt 1.415731279871813e-05, mass 1.142253008500561e+06, tot-E -6.402935117791155e+11; Picard 27.35 / 88, NC 0; end dt 1.415932e-05 |
| stream speed | 1.71e6 zone-cycles/cpu_s | 2.80e6 zone-cycles/cpu_s |

MI300A and A100 agree to round-off (about 1e-14 relative), not bitwise. Expect the same on H200.
Port check: a Caltech static smoke with `problem/stream=false time/tlim=0.05 time/nlim=10` should give
time=1.416097e-04 at cycle 10, Picard mean 3.379 / max 4, NC 0, and t=0 hst dt/mass/tot-E as above to about 1e-13.
Full logs are in ref/ and at /viper/ptmp2/jinma/csfix_1024/runs/r12_{static,stream}_new and
/raven/ptmp/jinma/csfix_1024/smoke/958cac78/{A,B}.

## Start rule (TASK point 4)
The viper chain (12327914..18) and the Raven chain (31291344..) stay queued until Caltech posts a NOTE that its run12
STARTED. After that, viper cancels its pending viper and Raven copies. If one of ours starts first, viper tells Caltech and
the USER decides.
