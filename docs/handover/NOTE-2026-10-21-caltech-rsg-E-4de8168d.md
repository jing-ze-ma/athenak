# NOTE 2026-10-21 caltech -> raven: RSG E 5126654 switched to 4de8168d + clip deck (reply to 93577914)

Done in place, no resubmit (link.sh and rdiff563.sh source run.cfg AT JOB START; all jobs still PENDING, none started):
- run.cfg BIN = builds_fin51/bin/athena_gpu_he_star_m1_4de8168d_nofma, BIN_MD5 9a56b4a109c575b7aa87a9a94d0124eb (verified)
- run.cfg IN = rsg_E_1020/rsg_E_1020_clip.athinput, md5 2ee06daba59f0c1395d129f4e2d16cd3
  = template rsg_E_start1020r_clip.athinput.in (md5 4c442452, verified) with ONLY the 3 @FILES@ lines substituted
  (-> /resnick/groups/carnegie_poc/jingze/rsg_E_1020/files); diff vs our old deck = the one `implicit_cs2_clip = rad` line.
- run.cfg md5 8bf41a58e289160495a920aa25f21b29 (old 3896d529 cfg kept as run.cfg.3896d529).
- Layout unchanged: 2 nodes x 4 H200, fresh t=0, tlim 1e9, same STOP rules; links 5126654 -> 55 -> 56. No new check.
- rdiff triplet 5140304 now runs on 4de8168d (det decks unchanged, no clip key).
Queue estimate for 5126654 still ~days (GPU partition ~3 d for 1 node). Will NOTE when it starts.
Accretor: nothing queued; 4de8168d accretor build (0ff5fc54) kept.
