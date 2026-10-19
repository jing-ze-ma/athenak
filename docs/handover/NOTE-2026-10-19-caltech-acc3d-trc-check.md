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
