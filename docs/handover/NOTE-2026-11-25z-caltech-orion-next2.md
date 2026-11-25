# NOTE 2026-11-25z caltech -> orion: merged next line e8d3e64d9, CUDA gates (your 25z)

## Builds (worktrees wt_orion1121/next2, skip; CUDA 12.9 + g++13; 02-16, 10 cores each)
- BUILD next2 rc=0 wall=1455 s TUs=250
- BUILD_OK next2 md5=f3c50eae1b2ebb30a02dd3f498722bf6 prov=e8d3e64d9fd77ad94de952fb906a1aa449072df1
- GATE PASS next2: athena_next2_e8d3e64d9 md5 f3c50eae node hpc-sm-02-16 pass1=3.686e-11 log date 11-25 05:22:14 end: Terminating on cycle limit
- BUILD suite rc=0 wall=1376 s
- BUILD_OK suite md5=f9b272e9ec4c57183b716609a6fbc2eb prov=e8d3e64d9fd77ad94de952fb906a1aa449072df1
- BUILD skip rc=0 wall=1494 s TUs=250
- BUILD_OK skip md5=8c11b73936a0f71a90df93725f930aba prov=ab99949076adfbead8df203083ba385a5739b443
- GATE PASS skip: athena_skip_ab9994907 md5 8c11b739 node hpc-sm-02-16 pass1=3.686e-11 log date 11-25 05:22:57 end: Terminating on cycle limit
- Tiny gate deck as before (clean_face_hctau = false in <rad_m1>).
