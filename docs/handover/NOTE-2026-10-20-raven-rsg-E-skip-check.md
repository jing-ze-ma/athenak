# NOTE 2026-10-20 raven -> caltech: RSG E third copy — SKIP the check job (USER)

USER (10-20): skip the pre-start check for the Caltech RSG E copy; start the chain directly.
- scancel 5126653 (check), and remove link 5126654's dependency on it:
  `scontrol update jobid=5126654 dependency=""` (links 5126655/56 keep afterok on the previous link).
- Reason: 563ab474 is fully gated (fin51-r-fix PUSHED = 563ab474; A100 + CPU + MI300A strict gbit clean; RSG E deck
  C==C2 + restart + forced retry BITWISE on A100), and the USER skips tests that aren't necessary.
- Deck md5: 393d443a (committed template) vs ae8d2e36 (Raven path-substituted) is expected — only the @FILES@ lines differ.
- Race rule unchanged: first copy to START wins; we tell you when another copy started, then you cancel yours.
- Heads-up: a BE fallback inside time2_reject=retry is coming (a dt-independent retry FATAL in near-vacuum cells on IDORT
  tests). When it is gated, a NOTE will name the new binary; switch at the next link (newest-binary rule).
Reply with a short NOTE: new state of 5126654 (pending reason / start estimate).
