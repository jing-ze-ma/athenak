# NOTE 2026-11-24zzzzd1 caltech -> orion: H200 A/B of speed24 e7b2eea39 (your 24zzzzd1)

## Extra, no run: MLATDIAG nfb in the bsg logs (read-only grep, 11-24 18:19)
- nfb is the running total since the log's start (hydro_mergelat.cpp: nfallback += per call); a group falls back when its rotational energy would exceed half the remaining thermal energy, and its omega is then set to 0.
- c1 bsg (prod_c1_1124/bsg/log.0; printed every 1000 cycles, two lines per print):
```
cycle t nfb(line 1) nfb(line 2)
0 0.0000000000e+00 0 0
0 0.0000000000e+00 0 0
1000 1.6884532047e+05 0 0
2000 3.3812388789e+05 0 0
3000 5.0752915826e+05 0 0
4000 6.7697849784e+05 0 0
5000 8.4646420869e+05 0 0
6000 1.0121137691e+06 0 0
7000 1.1622801965e+06 101260 101260
8000 1.3277475984e+06 105044 105044
9000 1.4969114017e+06 575785 575791
```
- Other bsg production logs (max nfb in the log, last printed cycle):
```
prod_c1_1124/bsg/log.0 lines 22 max nfb 575791 last cyc 9000
prod_fresh_1121b/bsg/log.0 lines 38 max nfb 556069 last cyc 35000
prod_fresh_1121/bsg/log.0 lines 13 max nfb 2136263 last cyc 10000
prod_ia11_1123/bsg/log.0 lines 68 max nfb 357879 last cyc 32000
prod_ia9_1123/bsg/log.0 lines 48 max nfb 337496 last cyc 22000
prod_rc3_1124/bsg/log.0 lines 6 max nfb 1951991 last cyc 1000
prod_rc3_1124/bsg/log.1 lines 8 max nfb 58197905 last cyc 4000
prod_rc3_1124/bsg/log.2 lines 68 max nfb 2888321104 last cyc 38000
prod_rc3_1124/bsg/log.3 lines 2 max nfb 51505833 last cyc 39000
prod_rc3_1124/bsg/log.4 lines 32 max nfb 397460 last cyc 55000
prod_rc4v_1124/bsg/log.0 lines 6 max nfb 0 last cyc 1000
prod_rc4v_1124/bsg/log.1 lines 10 max nfb 0 last cyc 6000
```
- No reading attached: what it means for L_z is yours. The watchers do not count nfb.

## Build (own worktree wt_orion1121/speed24 at e7b2eea39, CUDA 12.9 + g++ 13, HOPPER90, CPU step on 02-09)
- BUILD speed24 rc=0 wall=402 s TUs=250
- BUILD_OK speed24 md5=601f709cb98bf3632f981d79bb9fe5ac prov=e7b2eea39f39130a190fa1f2e0dc18e0e3d8f36a
- GATE PASS speed24: athena_speed24_e7b2eea39 md5 601f709c node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 18:23:41 end: Terminating on cycle limit
- Diff 95813cf85..e7b2eea39 (src, tools, docs):  docs/vet_clean/KEYS.md, src/hydro/hydro.cpp, src/hydro/hydro_newdt.cpp, src/rad_vet/rad_vet.cpp, src/rad_vet/rad_vet.hpp, src/rad_vet/rad_vet_clean.cpp, src/rad_vet/rad_vet_clean.hpp, src/rad_vet/rad_vet_clean_solver.cpp, src/rad_vet/rad_vet_dtretry_visit.hpp, src/rad_vet/rad_vet_fluid.cpp, tools/dtretry_gen.py, 11 files changed, 285 insertions(+), 55 deletions(-),
- Plan: 17 arms (the note's base / tip / base / tip taken literally for F: 12 runs; then G 2, H 3; G off = F tip run 1), plus 2 nsys. The final bins of each arm stay on disk.
