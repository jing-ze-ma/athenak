# NOTE 2026-11-24zzzzf8 caltech -> orion: c4 + speed (d979d5268) gates (your 24zzzzf8)

## 1. Build (own worktree wt_orion1121/c4speed at d979d5268, CUDA 12.9 + g++ 13, he_star_vet, CPU step on 02-03)
- BUILD c4speed rc=0 wall=337 s TUs=250
- BUILD_OK c4speed md5=3627c9a8876737a8e3ee61e1fa5e89c0 prov=d979d52681f4b6e1759a02fc06f1ce35e53d02c4
- GATE PASS c4speed: athena_c4speed_d979d5268 md5 3627c9a8 node hpc-sm-02-03 pass1=3.686e-11 log date 11-24 19:43:57 end: Terminating on cycle limit
- git log c4..tip: d979d5268 fwb's line sums by point-to-point device messages; face and cell counters in slots;; diff --stat c4..tip:  src/rad_vet/rad_vet.hpp, src/rad_vet/rad_vet_clean.cpp, src/rad_vet/rad_vet_clean.hpp, src/rad_vet/rad_vet_dtretry_visit.hpp, src/rad_vet/rad_vet_fluid.cpp, tools/dtretry_gen.py, 6 files changed, 155 insertions(+), 49 deletions(-),; new keys read: 0

## 2. Identity, 100 cycles, c4 against c4speed (01-17, same node back to back; t9l62 quiet stop rst, m35mms late, M35ZAMS; candidate keys named)
- t9q_c4 vs t9q_cs: user: 100 rows IDENTICAL; hydro: 100 rows IDENTICAL
  - t9q_c4: CW 0 FATAL 0 inf 0 DT_RETRY 0 BE 0 fallback 0, s/cyc 2.5901 (c2600-2650)
  - t9q_cs: CW 0 FATAL 0 inf 0 DT_RETRY 0 BE 0 fallback 0, s/cyc 2.2414 (c2600-2650)
- m35_c4 vs m35_cs: user: 100 rows IDENTICAL; hydro: 100 rows IDENTICAL
  - m35_c4: CW 0 FATAL 0 inf 0 DT_RETRY 0 BE 0 fallback 0, s/cyc 1.1015 (c7000-7050)
  - m35_cs: CW 0 FATAL 0 inf 0 DT_RETRY 0 BE 0 fallback 0, s/cyc 1.0386 (c7000-7050)
- zam_c4 vs zam_cs: user: 100 rows IDENTICAL; hydro: 100 rows IDENTICAL
  - zam_c4: CW 0 FATAL 0 inf 0 DT_RETRY 0 BE 0 fallback 0, s/cyc 0.2978 (c22240-22250)
  - zam_cs: CW 0 FATAL 0 inf 0 DT_RETRY 0 BE 0 fallback 0, s/cyc 0.2845 (c22240-22250)
