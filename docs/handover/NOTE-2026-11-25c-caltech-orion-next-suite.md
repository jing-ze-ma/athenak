# NOTE 2026-11-25c caltech -> orion: CUDA suite of the next line 79c1cabef (your 25c)

## Builds (own worktree wt_orion1121/next at 79c1cabef; CUDA 12.9 + g++ 13, HOPPER90; CPU steps on 02-03, 12 cores each; suite on the first free test node after it)
- BUILD next rc=0 wall=817 s TUs=250
- BUILD_OK next md5=9bebf6a3972983aedf3ba060d8a5c9cb prov=79c1cabef142b9f399c09c1a9370d3072558f9e2
- GATE PASS next: athena_next_79c1cabef md5 9bebf6a3 node hpc-sm-02-03 pass1=3.686e-11 log date 11-25 00:47:13 end: Terminating on cycle limit
- BUILD suite rc=0 wall=791 s
- BUILD_OK suite md5=710b7ae350c2404058790c5ba5d14127 prov=79c1cabef142b9f399c09c1a9370d3072558f9e2
- Tiny gate deck as for rel-e (copy with clean_face_hctau = false in <rad_m1>).
