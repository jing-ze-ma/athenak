# TASK raven -> DeltaAI: production-scheme cost breakdown, AG Car A + reduced BSG (1 short ghx4-interactive job)

From Raven worker M (bsgaccel-1011). The Raven gpudev slot is ~45 min out (1 node, 3 jobs ahead), so we ask you too.
**First result wins**: post a STARTED line in the NOTE when the job runs; Raven cancels its copy (31084675) then, and posts
WITHDRAWN here if Raven finishes first. Short test: <= 40 min wall, 1 node x 4 GH200, 1 rank per GPU, plain srun,
`export MPICH_GPU_SUPPORT_ENABLED=1`. If ghx4-interactive will not start within ~2 h, reply SKIP.

## 1. Build (report md5)
Fork branch `prodfinal-1011` @ **df483770** (code commit; build the tip if it moved and report the hash).
`build_inc_deltaai.sh hes_gpu` (PROBLEM=he_star_m1, GH200), as your 530c38e2 bsgrec smoke build.
Also build Kokkos Tools simple-kernel-timer (you have it: bsg_hrdet_1009/kt/profiling/simple-kernel-timer); set
`KOKKOS_TOOLS_TIMER_BINARY=1` so it writes `<host>-<pid>.dat` files, read them with kp_reader.

## 2. Files (docs/handover/bsgaccel-1011/ on bsg-files-1009)
- `agc.athinput.in` = the AG Car A production FRESH deck (prodfinal_1011/prod/agcar_A_fresh.athinput, Raven md5 of the
  resolved file f892a381) + one line `implicit_timers = 0` in `<rad_m1>` (read only when named; 0 = off).
  `sed "s#@FILES@#<your agcarA_files copy>#g"` (docs/handover/agcarA_files/: IC + Rosseland + Planck tables).
- `bsg.athinput.in` = the bsgrec R deck (bsgR.athinput, production scheme on the BSG truerepro2 input, original IC) +
  `implicit_timers = 0`. `sed "s#@BSGFILES@#<your 10-09 BSG files dir>#g"`.
- `ana.py run.log K0 10` = per-cycle wall from the `elapsed=... cycle=...` lines (steady / rotation / avg).

## 3. Runs (one job; each run its own -d dir)
Common keys:
- agc: `time/ndiag=1 output1/dt=1e-6 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 output5/dt=1e30`
- bsg (reduced mesh 4 blocks 256x64x64, 1 per GPU): `mesh/nx2=128 mesh/nx3=128 mesh/x2min=1.3089969389957472
  mesh/x2max=1.832595714594046 mesh/x3max=0.5235987755982988 meshblock/nx2=64 meshblock/nx3=64 time/ndiag=1
  output1/dt=1e-6 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30`

| run | case | extra | nlim |
|---|---|---|---|
| agc_plain | agc | - | 60 |
| bsg_plain | bsg | - | 40 |
| agc_tmr | agc | `rad_m1/implicit_timers=5` | 60 |
| bsg_tmr | bsg | `rad_m1/implicit_timers=5` | 40 |
| agc_kt10 / agc_kt60 | agc | KOKKOS_TOOLS_LIBS=libkp_kernel_timer.so | 10 / 60 |
| bsg_kt10 / bsg_kt40 | bsg | KOKKOS_TOOLS_LIBS=libkp_kernel_timer.so | 10 / 40 |

## 4. Report: NOTE-2026-10-11-deltaai-bsgaccel-prof.md
- Per run: rc, FATAL, NON-CONVERGED, Picard mean/max line, `ana.py run.log 5 10` output, nvidia-smi peak MiB + mean util.
- The two `<rad_m1> timers` lines of the tmr runs (verbatim).
- Per case: the kp_reader top-40 kernels of (long - short) per rank, mean over ranks (name, s, calls), plus "Total
  Execution Time" and "Total Time in Kokkos kernels" of both runs. Attach the kt_*.txt files.
- GPU-h used.
