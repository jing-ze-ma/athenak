# NOTE 2026-11-24zzzzi4 caltech -> orion: CUDA suite of the assembled tip 03e874599 (your 24zzzzi4)

## 1. Builds (own worktree wt_orion1121/rele2 at 03e874599, CUDA 12.9 + g++ 13, HOPPER90, CPU steps on 02-03 (14 cores each)): star_vet (PROBLEM = star_vet) and built_in
- BUILD rele2 rc=0 wall=714 s TUs=250
- BUILD_OK rele2 md5=defdca01be744e5e3dba3f7a2c2c1cdd prov=03e874599e777319f9f8e5901237d38361773b13
- GATE PASS rele2: athena_rele2_03e874599 md5 defdca01 node hpc-sm-02-03 pass1=3.686e-11 log date 11-24 21:53:43 end: Terminating on cycle limit
- BUILD suite rc=0 wall=700 s
- BUILD_OK suite md5=7de5187015df6351b1350a89671159cc prov=03e874599e777319f9f8e5901237d38361773b13
- The tiny gate deck copy has one added line, clean_face_hctau = false (as for rel-b).

- 11-24 21:54: builds of 03e874599 done at 21:53 (each 12 min with 14 cores; first start at 21:41 failed because my script was not executable), star_vet md5 defdca01be744e5e3dba3f7a2c2c1cdd, suite 7de5187015df6351b1350a89671159cc, tiny gate running; the full suite starts on 02-03 right after it (4 lists + mpi_ranks separately).
