# NOTE DeltaAI: vgdrot-1011 BSG timing + bitwise (TASK-2026-10-11-deltaai-vgdrot-bsg-timing)

**Picked up** 2026-10-10 22:40 CDT. BASE d28dacd8 built (athena_hes_gpu_d28dacd8aaa6, md5 315c33f7); NEW a749711a building.
Job 3357303 (ghx4-interactive, 1 node x 4 GH200, 0:30) submitted 22:49 CDT: base1, new1, base2, new2, BSG bsg_hr_dc5 reduced to
4 blocks 256x64x64 as in the TASK, 32 cycles, hst every cycle, bins only at cycle 0 + final, no rst. Results here when done.

**WITHDRAWN, acknowledged** 22:53 CDT (Raven 1b296899). Job 3357303 had started 22:51; only base1 had finished (rc 0, BASE only,
nothing to compare), so no GH200 cross-check to post. Cancelled after 46 s (~0.05 GPU-h); NEW a749711a build stopped unfinished.
Job script + analysis kept in /work/nvme/bivj/jma20/vgdrot_1011 (job.sub, analyse.py) if a GH200 cross-check is wanted later
(note vgdrot-1011 has moved on to 0d8eef6e). DeltaAI idle.
