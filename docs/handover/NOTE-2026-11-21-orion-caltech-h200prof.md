# NOTE 2026-11-21 orion -> caltech (USER): H200 profile arms PRF_rsg / PRF_hg (follows NOTE ed)

The USER asked whether any significant speed-up is specific to the H200. These two arms profile one production step of each
star. Files: docs/handover/caltech_runs_1121/h200prof/ (profc.sh + wrap_nsys.sh, wrap_ncu.sh, picktop.py, gaps.py, ncusum.py,
common.py; it reuses ../batch4/kname.py and ../rsg_v203.keys / ../hg_v203.keys).

**When and where:** production has priority. Use one whole node that production does not use (02-16 after R_LF/RA finish, or any
node the fresh production leaves idle). Never preempt or share a production node. Run the two arms one after the other on that node.
Each takes about 15-25 min: nsys ~2 min, then ~10 short ncu runs of ~1 min each.

**Commands** (same env convention as armc.sh; BIN = the B_v203 binary athena_v203_e688cb7bc, md5 7df2b648...):
```
cd docs/handover/caltech_runs_1121
HOLDER=<holder> NODE=<node> J=<J> OUT=<orion_runs_1121> h200prof/profc.sh PRF_rsg <B_v203 bin> <J>/rsg_E_1020/prodH/rst/rsg.01295.rst 5 rsg_v203.keys
HOLDER=<holder> NODE=<node> J=<J> OUT=<orion_runs_1121> h200prof/profc.sh PRF_hg  <B_v203 bin> <J>/hgglob_1021/prodG/rst/hegiant.00139.rst 5 hg_v203.keys time/dt_retry_host=m1
```
(= the V203L_rsg / V203L_hg recipe. No restart dumps are written. NR defaults to 4.)

**ncu:** profc.sh takes ncu from PATH, or else from the directory nsys is in, or else /usr/local/cuda*/bin. If ncu lives somewhere else,
set NCU=<path>. If the first ncu run fails (for example ERR_NVGPUCTRPERM: performance counters are blocked), profc.sh writes this in
c/rc.txt and skips the rest of phase C; the nsys part still counts. To skip ncu entirely, set NCU_ON=0.

**Please send back** caltech_results/PRF_{rsg,hg}/: res.txt, info.txt, n/{rc.txt,kname_r*.txt,gaps_r*.txt,stats_r0.txt,mem.log,run.log.gz},
c/{rc.txt,top.txt,ncusum.txt,k*.txt,k*.csv}. Please keep the .nsys-rep, .sqlite and .ncu-rep files on Caltech until 11-23; we may ask for one.
If a script fails, send the error and do not retry or rerun more than once.

-- orion (for the USER)
