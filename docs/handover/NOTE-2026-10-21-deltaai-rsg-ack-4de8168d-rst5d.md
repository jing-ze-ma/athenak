# NOTE 2026-10-21 deltaai -> raven: ack 0ba17385 (rst 5 d) + 14db112d (RSG E stays on 4de8168d)

- Binary: DeltaAI RSG E never switched; BINARY = athena_hesnf_gpu_4de8168d2b23 md5 c56693cb for all links. cdf71b98's RSG part ignored.
- Deck: rsg_E_clip_rst5d.athinput (44b2ef21) with @FILES@ re-substituted == shared template a8ab12df except the comment text on the
  <output5> dt line (value identical, 432000.0). Restarts also get output5/dt=432000.0 on the command line.
- Status: link 1 done (0 -> 24.5 d, NOTE d5730ffd). Link 2 = interactive 3426005 pending, Slurm est 10-21 ~17:16 CDT (low priority);
  batch 3425614 -> 3425615 pending (~10-26).
- Dump thinning (every 5th day + newest 10 d, only dumps the page already reduced): will do at the next session (pending USER OK
  on deleting data); quota is not tight before then (/projects ~120 GB used of 1 TB, run waiting in the queue).
