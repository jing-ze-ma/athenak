# TASK raven -> DeltaAI (2026-10-18): He giant 14p5M -- production-layout SMOKE + floor-variant screen (tests only)

USER: Raven queue too long, Caltech crowded, DeltaAI is up -> run here. Same tests as TASK-2026-10-18-caltech-hg145-smoke.md
(read it for the full spec: build, files, keys, report items). Race: whichever site starts first wins; post a one-line NOTE
as soon as a job STARTS so the other copies can be cancelled.

Differences for DeltaAI:
- Build prodmerge-1018-hg @ 4c494982 for GH200 (your usual recipe, PROBLEM=he_star_m1); export MPICH_GPU_SUPPORT_ENABLED=1.
- Test A (smoke): 1 node x 4 GH200, 4 ranks, 1 MeshBlock 704x64x64 per GPU; walltime <= 1 h. Keys and report as in the Caltech TASK sect. 3.
- Test B (fl15 screen): 1 GH200, keys as in the Caltech TASK sect. 4, but split into <= 1 h links: link 1 fresh with `-t 00:50:00`,
  links 2.. restart from the newest rst (`-r <rst> -t 00:50:00`) until t = 2.592e5 (at most 4 links). Keep the run dir; give an rsync path.
- Report: NOTE-2026-10-18-deltaai-hg145-smoke.md on bsg-files-1009.
