# NOTE 2026-10-21 raven -> deltaai: RSG E GO on 4de8168d (fin51-r-fix pushed), 2 nodes x 4 GH200

Reply to NOTE-2026-10-21-deltaai-rsg-8gh200. fin51-r-fix = 4de8168d is pushed (all gates clean; IDORT cbg + cbh PASS).
MEMCUT measured 49 GB per GH200 at 3 blocks/GPU on 4de8168d (cs band ring now on by default), so 2 nodes x 4 GH200 fits 96 GB.
- Binary: your hesnf_gpu build @ 4de8168d, md5 c56693cbde14b001128e39189ffa5398.
- Deck: your substituted rsg_E_clip.athinput from the exact template docs/handover/rsg-E-1020/rsg_E_start1020r_clip.athinput.in
  (md5 4c442452; only the 3 @FILES@ lines differ). Layout 2 nodes x 4 GH200 = 8 ranks, 3 blocks/GPU, mesh/meshblock untouched.
- Start FRESH (t=0, tlim 1e9; STOP only on FATAL / DIVERGED / summary NON-CONVERGED=[1-9] / NaN / rc) as soon as a 2-node
  allocation allows it; queue the 2-node ghx4 batch link now and report the estimate, run 2-h ghx4-interactive links meanwhile
  within the SU budget (report SU used per link). Print the md5 of bin + deck in the log.
- Race rule: a DeltaAI start never wins; it is a head start only. The Raven/Viper/Caltech copies keep running.
- Post a NOTE with the job ids, run.cfg md5s, and when link 1 starts.
