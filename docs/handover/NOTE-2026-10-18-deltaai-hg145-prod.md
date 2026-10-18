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

---
## Link 1 STARTED (2026-10-18 10:18 UTC, job 3403162, gh126); status at 10:28 UTC (10 min in)
- FRESH t = 0, clean so far: FATAL 0, Picard DIVERGED 0, DT_RETRY 0 retries. Picard 3.0 per solve (hst). L_top/L_in 0.9999 at t 1.23e5.
- **0.49 s/cycle**, dt 100 s flat -> t 1.20e5 (1.39 d) after 10 min wall = **~8.3 sim days per wall hour**. If dt holds, link 1
  (athena -t 47:45) reaches ~**390-400 d**. GPU memory 21.8 GB per GPU (of 96), the same as Test A.
- **Dumps: ~2.4 GB per sim day** (hydro_w 231 MB + m1 185 + m1_vet 185 per 0.25 d) = ~20 GB per wall hour. The mover is working:
  finished dumps go to /work/hdd/bivj/jma20/hg145_prod_1018/bin with cmp verification (mover.log).
  **Storage warning:** at this rate link 1 alone writes ~950 GB of dumps. That about fills the /work/hdd quota (1000 GB), with nvme
  (~280 GB free) as overflow, because the mover keeps any file it cannot copy on nvme. A link 2 at the same output cadence would NOT fit:
  dumps would need moving off DeltaAI, or a lower bin cadence (USER/Raven decision).
- One-line NOTE at link end (~47.8 h from start).

---
## Link 1 (3403162) ENDED 14:04 UTC: FATAL on rst write, nvme QUOTA FULL (DeltaAI's fault, not the physics); resubmitted 3405272
- At cycle 24150 (t 2.415e6, 28.0 d) writing hegiant.00056.rst: `Disk quota exceeded` -> restart.cpp:1254 FATAL
  ("vet_rotate_avg history cache not written correctly"); job hung after the FATAL and was cancelled (3 h 46 min, ~15 SU).
- Cause: the deck writes a **~7.5 GB rst every 0.5 d (~15 GB per sim day)** and my link script pruned rst only at link end.
  I had budgeted the bin dumps only. 56 rst files (~420 GB) filled /work/nvme/bivj.
- Physics up to the crash was clean: FATAL 0 before the I/O error, DIVERGED 0, NC 0. dt 100 s flat. Pace dropped from 0.49 to ~0.68 s/cycle
  after ~20 d as Picard rose 3 -> 5 per solve. **L_top/L_in drifted 1.00 -> 0.95 by ~27 d** (worth a look on Raven).
- Recovery: rst_pick.py chose hegiant.00055.rst (t 2.376e6 = 27.5 d, cycle 23760, valid); the truncated 00056 is quarantined in rst/invalid/.
  Pruned rst 00001-00052 (TASK rule keep newest 3; kept 00000 = IC + 00053-55). All bin dumps kept (0-2.4e6, 61 GB on /work/hdd).
  hst has rows for 2.376e6-2.415e6 that the restart will write again (duplicate times; drop the earlier rows when analysing).
- Fix: an in-job rst pruner (keeps IC + 4 newest, the newest possibly still being written) runs every 5 min; the restart uses
  rst_pick.py (validity check). Storage now stays bounded: rst <= ~40 GB on nvme, bins ~2.4 GB/sim day to hdd.
- **Resubmitted as job 3405272** (ghx4, 4 GH200, `-r hegiant.00055.rst time/tlim=1e9`, **44 h** wall so the total stays ~192 SU = the one
  approved link). Slurm start estimate ~17:42 UTC 10-19 (~27.5 h).

---
## 15:25 UTC: also chaining 2 h ghx4-interactive links (USER: "maybe you can try interactive node")
- ghx4 queue is ~2200 jobs deep; new 1-node jobs estimated ~4 days out. ghx4 link resubmitted as **3405590** (replaces 3405272; it needed the
  new run lock). Interactive chain started: first link **3405591** (2 h, est. start ~03:21 UTC 10-19). Interactive costs x2 = 16 SU
  per link, ~1 h 50 min of run each (~9 sim days at the late pace).
- Both share run/ through a lock: an interactive link exits if the run is locked; the ghx4 link waits (<= 130 min) for the lock. The chain
  stops submitting once the ghx4 link runs, and the ghx4 link trims its own TimeLimit so the WHOLE run stays within **192 SU**
  (= the one approved 48 h link; 15 SU used so far). Each link restarts from the newest VALID rst (rst_pick.py).

---
## 20:08 UTC: CANCELLED per USER (Raven NOTE 4820331f)
Chain stopped (STOP_chain + killed), ghx4 3405590 and ghx4-interactive 3405591 cancelled while PENDING (0 s run). No movers or pruners are
alive (they only ran inside jobs). **Total SU for the hg145 production on DeltaAI: 15.1** (link 3403162, 3 h 46 min on 4 GH200). Run kept:
/work/nvme/bivj/jma20/hg145_prod_1018/run (rst 00000 + 00053-55, last valid 27.5 d, hst, logs; 32 GB) and dumps 0-28 d at
/work/hdd/bivj/jma20/hg145_prod_1018/bin (64 GB, 113 dumps). Page bundle work: /work/hdd/bivj/jma20/hg145_prod_1018/page. Nothing deleted.
