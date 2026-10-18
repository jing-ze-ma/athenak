# NOTE DeltaAI -> Raven, Caltech (2026-10-18 04:40 UTC): hg145 smoke STARTED on DeltaAI (TASK 9537ef9d)

**Test A RUNNING** (job 3403091, 1 node x 4 GH200, started 04:36 UTC). Test B link 1 (job 3403094, 1 GH200) starts in a few minutes.
Built prodmerge-1018-hg @ 4c494982 -> athena_hes_gpu_4c494982ffb3, md5 385860369a7d3d3801df2b4649c0ced7. All task file md5s match.
**Caltech: please cancel 4925620-22** (hg145 build/A/B). Results follow in this NOTE (A within ~1 h, B in ~2-4 h).
Planned deviation for A: if C's hegiant.00001.rst is the end-of-C file (cycle 300, so nlim=300 would run 0 cycles),
R restarts from hegiant.00000.rst (t=0) and runs 300 cycles instead (logged).

---
## RESULTS (2026-10-18 04:50 UTC): Test A and Test B both DONE, both CLEAN. Total cost ~0.4 GPU-h.

### Test A: production smoke (job 3403091, gh046, 4 GH200, deck hg145_prod, keys exactly as TASK sect. 3)
| | C (fresh, nlim 300) | R (`-r C/rst/hegiant.00001.rst time/nlim=300`) |
|---|---|---|
| rc / FATAL / NaN / Picard DIVERGED | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |
| summary NON-CONVERGED | 0 | 0 |
| DT_RETRY | 300 checkpoints, **0 retries** | 85 checkpoints, 0 retries |
| Picard mean / max | 3.06 / 4 | 3.00 / 3 |
| dt | 100 s flat (cycle 300 = t 3.0e4) | 100 s |
| s/cycle | 0.52 (cycles 75-300, incl. hst every cycle + bins) | 0.70 (short run, rst read) |
| GPU peak memory | 21.9 GB / GPU (of 96) | 21.7 GB |
- No deviation needed: hegiant.00001.rst was written mid-run at t 21500 (cycle 215); 00002 = end of C (t 3e4, cycle 300).
- **R vs C BITWISE: yes.** All 85 overlapping hst rows (hydro + user) byte-identical; `cmplast145.py C R` -> CMPLAST BITWISE
  (hydro_w, m1, m1_vet final bins at t 3e4, every variable). Raw `cmp` of bin files differs only in the PAR_DUMP header length
  (different command lines); data identical.
- Startup lines: `he_star_m1: mlt_closure = adaptive, relax time 1500 s, every 1 cycle(s), taper 0.01, update clock 0, max L_sub = 1.16688e+37`;
  `... r_top = 2.47391e+13, sponge rate = 0.001 above r = 2.22652e+13`; floor: `floor from r = 8.27659e+12 (1.04154 R_ph)`
  (atm from cell 294, R_ph 7.94683e+12, 39 cells above floor). fofc fired 0 cell-stages.
- L_top/L_in: 0.847 (t 0) -> 1.004 by 0.17 d -> 0.9996 at 0.35 d. KE_int 9.1e37 at t 3e4.

### Test B: fl15 screen (job 3403094, 1 GH200, deck hg145_mkfl15_nobh, keys exactly as TASK sect. 4)
- Finished in ONE link: 1815 cycles to t = 2.592e5 in 4.7 min (~0.15 s/cycle; dt 114-143 s). Links 2-4 exited unused (DONE).
- rc 0, FATAL 0, NaN 0, DIVERGED 0, summary NON-CONVERGED 0, DT_RETRY 0 retries, Picard mean 3.04 max 5. GPU peak 2.9 GB.
- hst (time, L_top/L_in, KE_int):
  0 d 0.872 0 | 0.25 d 1.001 2.8e37 | 0.5 d 0.999 4.7e37 | 1.0 d 1.000 5.9e37 | 1.5 d 1.000 7.5e37 | 2.0 d 1.001 8.8e37 |
  2.25 d 1.004 1.1e38 | 2.5 d 0.998 2.2e38 | 2.75 d 0.999 4.2e38 | 3.0 d 0.998 6.3e38.
  L_top/L_in min 0.755 (startup transient) max 1.019. Mdot_top 8.7e14 -> 8.1e14. M_tot constant to 7 digits.
  **Flag for Raven:** KE_int roughly flat-linear to 2 d, then grows ~7x between 2.0 and 3.0 d (8.8e37 -> 6.3e38) while L_top/L_in stays 1.00.
  Not analysed here; the bins are kept so you can look at the interface wave at 1.015-1.04 R_ph.
- Run dir (rsync): `jma20@dtai-login.delta.ncsa.illinois.edu:/work/nvme/bivj/jma20/hg145_1018/B/` (791 MB: bin/ every 4320 s, hst, rst 00000-00002).
  Test A: `/work/nvme/bivj/jma20/hg145_1018/A/{C,R}` (rst 5-8 GB each; ask before syncing).

Binary athena_hes_gpu_4c494982ffb3 (prodmerge-1018-hg @ 4c494982), md5 385860369a7d3d3801df2b4649c0ced7. Nothing left running on DeltaAI.
