# NOTE 2026-10-21 raven -> caltech + deltaai: RSG E writes restarts every 5 d on ALL sites (USER: write less output)

New shared template docs/handover/rsg-E-1020/rsg_E_start1020r_clip_rst5d.athinput.in (md5 a8ab12df) = clip template 4c442452 with ONLY
<output5> (rst) dt = 432000.0 (5 d, was 0.5 d). Output-only; physics identical. Raven R2 and Viper V2 already switched (not started).
- CALTECH (RSG E 5126654, not started): switch IN to this template (3 @FILES@ lines substituted as before) before it starts; report the
  run.cfg + deck md5. If it has already started: add output5/dt=432000.0 on the command line at the next link instead.
- DELTAAI: already on 5 d (your deck 44b2ef21 = same change). Dumps stay daily (hydro_w + m1, the page source); with your 1 TB quota,
  thin dumps OLDER than the page's last reduction to every 5 d (keep every 5th day + the newest 10 d), after the page has reduced them.
- HG global unchanged (already rst 5 d). HG wedge unchanged.
