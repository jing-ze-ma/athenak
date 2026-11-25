# NOTE 2026-11-24zzzzf8 caltech -> orion: c4 + speed (d979d5268) gates (your 24zzzzf8)

## 1. Build (own worktree wt_orion1121/c4speed at d979d5268, CUDA 12.9 + g++ 13, he_star_vet, CPU step on 02-03)
- BUILD c4speed rc=0 wall=337 s TUs=250
- BUILD_OK c4speed md5=3627c9a8876737a8e3ee61e1fa5e89c0 prov=d979d52681f4b6e1759a02fc06f1ce35e53d02c4
- GATE PASS c4speed: athena_c4speed_d979d5268 md5 3627c9a8 node hpc-sm-02-03 pass1=3.686e-11 log date 11-24 19:43:57 end: Terminating on cycle limit
- git log c4..tip: d979d5268 fwb's line sums by point-to-point device messages; face and cell counters in slots;; diff --stat c4..tip:  src/rad_vet/rad_vet.hpp, src/rad_vet/rad_vet_clean.cpp, src/rad_vet/rad_vet_clean.hpp, src/rad_vet/rad_vet_dtretry_visit.hpp, src/rad_vet/rad_vet_fluid.cpp, tools/dtretry_gen.py, 6 files changed, 155 insertions(+), 49 deletions(-),; new keys read: 0
