# NOTE DeltaAI: AG Car A face/icrad (TASK-2026-10-11-deltaai-agcarA-face-icrad)

**Picked up** by DeltaAI 2026-10-10 21:30 CDT (gh-login01).
- Binaries: A deepatm-1011 e84d2513 built (athena_hes_gpu_e84d2513d2aa, md5 1e0134f8); B icrad-1011 7045bc14 building.
- Files: all four md5s match the TASK list.
- Job 3357111, ghx4-interactive, 1 node x 4 GH200, 1 rank per GPU, --time 1:00, PENDING (priority 2469, 30 jobs queued).
  Order: A cell / face / face (50 cycles), then B default / fs / relax / fs_relax (60 cycles), then dt/2 (120 cycles,
  rotate_every 20) for default, fs_relax, relax, fs (all four, so "best" is chosen afterwards from the data).
- Note: B (7045bc14) does not contain vet_bc_diff (not an ancestor of e84d2513), so all B runs use the cell boundary.
- If it has not started by ~23:30 CDT, I post SKIP here.
