# NOTE 2026-11-24zzzzc2 orion -> Caltech: profiles of the star steps and the t9l62 merged-base arms (USER: 24 h per star)

**Order:** production first; then c4 (24zzzzb8 items 2 and 3), the dt limiters (24zzzzb9), CSHERO D (24zzzzb6), sdorder5 (24zzzzb5); then section 1 priority 1, section 2, the rest. Run an item earlier only on a node nothing above can use. No key of any production run changes. If a rst named below no longer exists, say so and skip that arm. MLQ1 switches the merge on at a restart, so its first cycles may show a start-up signal: report cycles 0-10 apart from the rest. Text by SPEED24:


USER goal: t9l62, m35mms, M35ZAMS within 24 h each on one node. Idle nodes only, production has priority, copies only (nothing in a production dir is touched), no rst / bin kept.
Binary for all arms: **95813cf85 (md5 30d01765)**, decks / keys = candidate_1124 (every candidate key named on the command line of a restart, as in the CAND_1124 arms; cfl 0.3, conv 1e-5 named).

## 1. PRF: nsys profile (profc.sh of docs/handover/caltech_runs_1121/h200prof, NCU_ON=0, NR=4, NC=12; `time/nlim`, `job/final_restart=false` as the script sets)
| arm | start state | priority |
|---|---|---|
| PRF_t9q | t9l62 quiet stop rst (the PQ state, c2553) | 1 |
| PRF_m35 | m35mms late: the CAND_1124 m35mms_R rst (c6950); mg 4 | 1 |
| PRF_zam | M35ZAMS: the m35zams_R rst (c22155), mg 2 | 2 |
| PRF_t9o | t9l62 00015 (PO state) | 2 |
| PRF_m35c | m35mms newest c1 rst, a COPY, if t > 6e4 | 3 |
- Name `rad_vet/clean_verbose=2` in each (host printf; only for matching passes).
- Return per arm: res.txt, n/kname_r0..3.txt, n/gaps_r0..3.txt, n/stats_r0.txt, the `^CLEAN (PASS|cycle)` lines, and the plain s/cyc of the same state from the CAND arm (nsys slows the step: the kernel SHARES are what we use).
- Questions: per family (fam.py), what is the t9l62 fixed 1.79 s (rays are 0.53)? how much is line solve + apply + reductions on m35mms / M35ZAMS? GPU idle share (MPI / sync) on M35ZAMS at 1.6 M cells per GPU?

## 2. ML: t9l62 with the innermost cells merged 1 x 2 x 2 (existing key, bsg uses it; no build)
Extra keys (restart: "added from the command line" NOTEs expected, grep them):
`hydro/merge_lat_rc=2.022e12 hydro/merge_lat_f=2 hydro/merge_lat_blockend=exchange hydro/merge_lat_diag=10`
(merge_lat_blockend must be named: a restart without it takes flat.) rc = 1.2 x1min = innermost 27 of 512 shells. IC prediction (the code's formula; it reproduces the measured 96.85): dt 96.85 -> **152.5** (radial bound at the base cell; first unmerged shell 166).
| arm | start | cycles | keys |
|---|---|---|---|
| MLQ0 | t9l62 quiet stop rst (PQ) | 100 | candidate, dt_report=10 |
| MLQ1 | same | 100 | + the merge keys, dt_report=10 |
| MLQ1r | same | 50 + rst + 50 | + the merge keys: hst bitwise vs MLQ1 (delete the rst) |
- Read: dt, DTREP lines (cell, dir, speeds), passes and its per cycle, s/cyc (the larger dt may raise passes), CW / BE / FATAL / DT_RETRY, dt_retry gain, MLATDIAG M / E / P / L drift per cycle vs MLQ0, L_top / L_in and radial E, T profiles vs MLQ0 at equal model time (100 x 96.85 s; scored on hst / existing bins), floor_pay, clean_wb counters.
- A pass is not a claim that the answer is unchanged: the lateral resolution of the 27 innermost shells changes. USER decision.
**Late (E/F): no arm is asked yet.** The merge at 1.2 x1min cannot raise the late dt beyond what the next binder allows. With the binary of 11-23 on 00037 (cfl 0.6, its DTREP), besides the wall cell at 96.4, dir-3 cells at r 2.57e12 (1.52 x1min, theta 47, v -2.3e7) sat at 94.5-96 and a dir-1 cell at r 4.4e12 (v -9.4e7) at 93.7 (both x0.3/0.6). To move those rc would need > 1.55 x1min for the first and the second is radial. So the late gain is claimed as 0.
## 3. Late state on the CURRENT scheme (needed for the E/F projection)
t9l62 restarts of the candidate scheme beyond 4.6e5 s do not exist. Two options; ONE of them when a node is free, after 1 and 2:
- LW: rst_p1_seg 00038 (t 1.598e6, dt 127.5, the newest of the old-scheme files with 00039) with candidate keys + `rad_vet/clean_conv_max=1e-5` named, warm 300 cycles (about 20 min) with clean_verbose=2, rst at the end; then 100 cycles from it. State "settled" = FCT rounds per cycle and passes per cycle flat over the last 100 (if 49-285 rounds are still falling at 300, say so and stop). Report s/cyc, passes, FCT rounds, dt, limiting cell over the last 100; keep the LW rst until T9LATE says otherwise.
- or the first 5e5-7e5 rst of the T2 / c1 t9l62 run when it exists (zero cost).
-- orion (coordinator), drafted by SPEED24
