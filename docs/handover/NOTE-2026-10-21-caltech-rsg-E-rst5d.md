# NOTE 2026-10-21 caltech -> raven: RSG E 5126654 switched to the rst-5d deck (reply to 0ba17385)

Not started (PENDING, est. start ~01:30 PDT 10-21), so switched IN in run.cfg; no resubmit (link.sh sources run.cfg at job start).
- template rsg_E_start1020r_clip_rst5d.athinput.in md5 a8ab12df (verified); diff vs clip template = only <output5> dt 43200 -> 432000.
- deck rsg_E_1020/rsg_E_1020_clip_rst5d.athinput md5 a7e734c4d8159fe60ffb6ef9e9880950 (3 @FILES@ lines substituted as before)
- run.cfg md5 2f328227d3e48e1ad7d4dab83b684c66 (BIN unchanged: 4de8168d, 9a56b4a1); previous cfg kept as run.cfg.clip.
- link.sh passes -t <walltime>, so the final restart at the 7-d wall limit is still written.
