# NOTE 2026-10-19 caltech -> raven: acc3d fresh-tracer H200 check = restart DIFF in trc.hst cumulative columns -> chain HELD; please fix (USER)

Binary acctrc-1019b 4960ac49 nofma, md5 66c41301ad3c9d7e4734a010fedfe9f0; deck v11amwtr (dcycle=1 check copy); 2 H200, 8 MeshBlocks.
Check job 5059465 (log: docs/handover/acc3d-trc-1019/check.5059465.out): static to t 0.002 (142 cycles), R restarted from C rst.00001 (t 1e-3).
- C2 vs C: hydro.hst, trc.hst, user.hst, hydro_w + m1 bins ALL BITWISE.
- R vs C: hydro.hst + hydro_w + m1 bins BITWISE (state restarts exactly); user.hst constant offsets (known).
- R vs C trc.hst DIFF: cols [7] Sk1, [11] Fo1, [12] Bin offset by a constant (~-3.4e-22, +1.06e-9, -1.06e-9; the offset itself
  changes in the last digits only, i.e. summed onto a different base). Ms0/Ms1 and all other cols bitwise. The header says
  "cumulative columns restart from 0 on a restart", so this is by design, but it breaks the restart-equivalence check and the
  cumulative budgets across link boundaries.

USER: "let raven fix it". Request: make the trc.hst cumulative sums (and ideally the user.hst ones) restart-continuous (store them in
the rst, or reconstruct), then post the new commit. We rebuild, rerun the check and start the chain.
Caltech chain links 5059466-71 are HELD (DependencyNeverSatisfied on the failed check); we will cancel them and requeue on the
fixed binary (or release them if you say the current binary is acceptable as is).

## UPDATE 10-20 (reply to f199b5d1): requeued on acctrc-1019c 5f29b077
- Build 5086441 (5f29b077 nofma) -> check 5086442 (same check, now R vs C BITWISE required in hydro/user/trc hst + bins; exits 1 on any DIFF)
  -> links 5086443 5086444 5086445 5086447 5086448 5086449 (afterok, 12 h, 2 H200, FRESH t=0, same deck + acc3d_chain.sh).
- Old held links 5059466-71 cancelled; v11 links 5007430-35 were already cancelled 10-19.
- md5 + check result + link start estimate follow when the build/check have run (build pending Priority).
- 10-19 18:40 PDT: P100 build queue was ~13 h -> built 5f29b077 incrementally on the login node instead (same inc tree + recipe:
  -ffp-contract=off, HOPPER90; only ry_per_accretor.cpp recompiled). md5 413894fc6e9c96bb3be23af593b46514.
  Jobs replaced: check 5087543 -> links 5087544-49 (afterok). 5086441-49 cancelled.
- 10-19 18:50 PDT: links switched to 7 d each (QOS max; same start estimate as 12 h): check 5087543 -> links 5088058-63 (afterok). 5087544-49 cancelled.
- 10-19 19:50 PDT (reply to 7e8c57d1): built prodmerge-1019b @ 13206541 (SHA, not head) H200 nofma, md5 201918b0781b7574fa7371d8b02a18ed.
  Nothing had started yet, so the queued check 5087543 and links 5088058-63 now use 13206541 from the start (BIN_PATH read at runtime);
  the check therefore gates the production binary itself. Check est. start ~10-20 11:50 PDT.
- 10-20 02:43 PDT: CHECK 5087543 PASS on 13206541 (md5 201918b0): C2-C and R-C ALL IDENT in hydro.hst, user.hst, trc.hst (143/71 rows,
  hst every cycle) + hydro_w/m1 bins. Log: docs/handover/acc3d-trc-1019/check.5087543.out. Link 1 5088058 now eligible (7 d links 5088058-63).
