# TASK raven -> Caltech: WV wavefix-1013, the FULL IDORT suite on ONE new binary (359 run dirs)

From Raven worker WV, via the coordinator. USER: "GO on all three fixes" (wave tests rswtv / rlatt), full IDORT suite,
compared against FIN3. Bubbles and 2-source are excluded (as in the mbband bundle). Not the same as OD's mbband or OT's
enth_be TASKs: this one is **wavefix**.

## 1. Binaries (report both md5)
- OLD = plain CPU build of **9972a56d** (rotavg-rst2-1013, the production hash; NOT 281dd686: its rst files differ in
  layout, and the BIT step compares rst data).
- NEW = branch **`wavefix-1013` @ `41c8b210`** (fork; = 9972a56d + 3 commits:
  e84d7b9f implicit_bc_advect_wall, d458ff31 vet_gd_band_dir, 41c8b210 implicit_hr_q_face (default OFF, not used here)).
- Build both with EXACTLY the same flags (compiler, OpenMP, MPI, Release, PROBLEM = built_in_pgens, and the same
  Athena_M1_VGD_NMAX if you set one for the p2 decks).

## 2. Run
`bash run_wavefix.sh <OLD> <NEW> 16 10` on one node (debug-QoS resubmits are fine; finished steps are skipped).
Bundle: docs/handover/open-1013/caltech_wavefix/ (decks = the mbband / CVALL decks: 45 cart + 194 sp + 4 ss + 10 p2 + 106 lad).

Per dir:
1. BIT, 10 cycles each, same threads: OLD; NEW (a) = deck as is (fixes 1 + 3 on by default); NEW (z) = keys off
   (implicit_bc_advect_wall=false vet_gd_band_dir=false). bit.txt = "<OLD vs a> <OLD vs z> <dir> <s OLD> <s a> <s z>".
   (z) MUST be BITWISE in every dir: report any (z) DIFF at once.
2. FULL (a): NEW to the end in <dir>/full_a, only for dirs whose (a) is DIFF (two-stage dirs: stage A, then -r, as in CVALL).

## 3. Report: docs/handover/NOTE-2026-10-13-caltech-wv-wavefix.md
- md5, job ids; bit_all.txt with the per-arm counts; every FAIL / non-zero rc line; the FULL_a lines.
- `tar czf wavefix_runs.tgz runs --exclude='*.rst' --exclude='bitO' --exclude='bitA' --exclude='bitZ'`, pushed under
  docs/handover/caltech_wavefix/ (split if > 90 MB). Post interim bit_all counts when the BIT step is done.

## Coordinator notes
- METHOD fixes, user GO. Expected (a) DIFF: decks with a reflecting gas x1 wall under a Marshak radiation end (rsw*, rlat*,
  rswt*, ...), and vet_gd wedges with nx2 != nx3 whose band was capped (rlatt*, rlat*). Production is unchanged (user gas
  BCs; multi-rank keeps the old band cap); WV checks reduced AG Car / bsgR bitwise on Raven.
- Priority: after the running od / ot tasks unless the coordinator says otherwise.
