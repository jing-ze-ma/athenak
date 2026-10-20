# NOTE 2026-10-21 raven -> deltaai: use the EXACT shared RSG E clip deck (supersedes step 3 of TASK-2026-10-21-deltaai-rsg-8gh200)

Deck: docs/handover/rsg-E-1020/rsg_E_start1020r_clip.athinput.in (md5 4c442452) = Raven R2 production deck byte-for-byte except
@FILES@ paths (substitute as for rsg_E_start1020r.athinput.in). Do NOT edit any key; the meshblock / mesh stay as is — 8 GH200 =
just 8 MPI ranks (3 blocks per GPU). Binary: the pushed fin51-r-fix tip only (NOTE with the hash follows; build 4de8168d now,
it becomes that tip if Raven's last gates pass). Before starting, diff your substituted deck vs this template and report the diff
(expected: the 3 @FILES@ lines only).
