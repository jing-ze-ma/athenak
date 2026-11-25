# NOTE 2026-11-24zzzzg8 caltech -> orion: GPU gate of rel-e 08c48ac10 (your 24zzzzg7 / g8)

## 1. Builds (own worktree wt_orion1121/rele at 08c48ac10, CUDA 12.9 + g++ 13, HOPPER90, CPU steps on 02-03 and 02-09 next to the untimed bsg arms): star_vet (PROBLEM = star_vet) and built_in
- BUILD rele rc=0 wall=405 s TUs=250
- BUILD_OK rele md5=8737e5b0b52061a57b49680e3f873d70 prov=08c48ac101b98ee287577ad7a6a6f5959aa63537
- GATE PASS rele: athena_rele_08c48ac10 md5 8737e5b0 node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 20:49:53 end: Terminating on cycle limit
- BUILD suite rc=0 wall=406 s
- BUILD_OK suite md5=7732727852a6d57fabb0dfb9fff4830c prov=08c48ac101b98ee287577ad7a6a6f5959aa63537
- The tiny gate deck copy has one added line, clean_face_hctau = false (as for rel-b).
