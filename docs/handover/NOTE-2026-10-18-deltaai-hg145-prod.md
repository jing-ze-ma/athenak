# NOTE DeltaAI -> Raven (2026-10-18 05:02 UTC): hg145 PRODUCTION link 1 QUEUED (TASK 70c0836b)

- **Job 3403162**, ghx4, 1 node x 4 GH200, walltime 48 h (ghx4 MaxTime), athena `-t 47:45:00`. Slurm start estimate ~11:26 UTC
  (~6.5 h; ghx4 estimates have been unreliable both ways). Run dir /work/nvme/bivj/jma20/hg145_prod_1018/run, FRESH t = 0.
- Binary athena_hes_gpu_4c494982ffb3 (md5 385860369a7d3d3801df2b4649c0ced7, checked at job start); deck hg145_prod.athinput.in d70305a7
  with @FILES@ only; layout and env exactly as Test A. tlim 1e9; restart links use `-r <newest rst> time/tlim=1e9`; rst keep newest 3.
  STOP rule exactly as the TASK (summary NON-CONVERGED=[1-9], FATAL, "Picard solve DIVERGED", NaN in the log tail, rc != 0).
- **Deviation (USER 10-18): ONE link only for now, no afterany links yet.** Reason: SU and storage.
  - SU: 882 left on bivj-dtai-gh; a 48 h link = 192 SU (3 links = 576). Links 2+ after the USER sees link 1.
  - **Storage: /work/nvme/bivj quota has ~300 GB free (699 of 1000 GB used).** Deck outputs = ~600 MB per 0.25 d (hydro_w 231 MB +
    m1 185 + m1_vet 185) = ~2.4 GB per sim day; at Test A speed (dt 100 s, 0.52 s/cycle, ~7.5 min wall per sim day) that is
    ~20 GB per wall hour, enough to fill nvme in ~15 h. Outputs are NOT reduced. Instead an in-job mover copies each finished bin dump
    (unmodified for 5 min) to **/work/hdd/bivj/jma20/hg145_prod_1018/bin** and removes the nvme copy only after `cmp` confirms the copy.
    No dump is deleted. hst and rst stay in run/. Log: hg145_prod_1018/mover.log.
- Next: a NOTE when link 1 starts (s/cycle, dt, Picard, memory) and a one-line NOTE at link end.
