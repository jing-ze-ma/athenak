# TASK raven -> Caltech (2026-10-18): He giant 14p5M -- production-layout SMOKE + one floor-variant screen (tests only, NO production chain)

Raven queues are hours out; USER: "caltech is up" -> run these two GPU tests on Caltech H200 now. Production itself stays HELD
(separate TASK after the USER's go).

## 1. Build
Branch **prodmerge-1018-hg @ 4c494982** (pushed) = fin51-r-fix bb360d8f + mltp-1017 (opt-in problem/mlt_closure) + ic-mesascaf
tools. Same cmake recipe as your fin51-r-fix build (PROBLEM=he_star_m1, CUDA H200). Report md5.

## 2. Files: docs/handover/hg145-prod-1018/ (md5)
- `hg145_prod.athinput.in` d70305a7 and `hg145_mkfl15_nobh.athinput.in` a0dceb7b -> `sed "s#@FILES@#<files dir>#g"`; change no key.
- `ic_hegiant_14p5M_mk.txt` f9678c2a, `ic_hegiant_14p5M_mk_fl15.txt` 13c23d32, `rosseland_mesakap_afe2dd8da2bd.txt` e2868b75,
  `planck_ext2_x0.0009496_z0.0215304_416a61196c64.txt` 23ee543f, `cmplast145.py` (restart compare: `python cmplast145.py <C> <R>`).

## 3. Test A -- production smoke (1 node x 4 H200, deck hg145_prod)
4 ranks, 1 MeshBlock 704x64x64 per GPU (mesh 704x128x128), your CUDA_VISIBLE_DEVICES wrapper, --cpu-bind=cores, UCX_TLS=^xpmem.
Keep vet_gd_twin_lowmem = true.
- C (fresh): `time/nlim=300 time/ndiag=1 output1/dt=1e-6 output2/dt=3600 output3/dt=3600 output4/dt=1e30 output6/dt=1e30 output5/dt=21500`
- R (restart): `-r <C>/rst/hegiant.00001.rst time/nlim=300` -> hst tail rows + final bins BITWISE vs C (cmp + cmplast145.py).
Report: rc, FATAL/NaN/DIVERGED/summary NON-CONVERGED, DT_RETRY count, Picard mean/max, dt, s/cycle, GPU peak memory, the startup
lines `mlt_closure = adaptive`, sponge and floor radius lines; R vs C bitwise yes/no.

## 4. Test B -- floor/sponge variant screen (1 H200, deck hg145_mkfl15_nobh = dfloor 1e-19, pfloor 1e-12, he_sponge_dmax 2e-15, fl15 IC)
1 rank, keys: `mesh/nx1=704 meshblock/nx1=704 mesh/nx2=16 mesh/nx3=16 meshblock/nx2=16 meshblock/nx3=16 problem/he_seed=0.0
problem/mlt_closure=adaptive time/tlim=2.592e5 output1/dt=300 output2/dt=4320 output3/dt=4320 output5/dt=1.728e5 output6/dt=1e30`
(3 d; ~3 h on an A100). Report rc/FATAL/NC as above and hst columns vs time (L_top/L_in, KE); keep run dir with bins for Raven
to analyse (the interface wave at 1.015-1.04 R_ph). rsync-able path please.

## 5. Report: NOTE-2026-10-18-caltech-hg145-smoke.md on bsg-files-1009 (A as soon as done, B when done).
STOP rule: only end-of-run summary NON-CONVERGED=[1-9], FATAL, "Picard solve DIVERGED", NaN, rc != 0.
