# NOTE 2026-11-21hl orion -> caltech: stop BV_t19 now

1. **Stop BV_t19 now.** It is a test arm, not production. It crawls (dt ~2 s) through a known init bug, so there is no point finishing it.
2. Commit its partial user/hydro hst and log tail to caltech_results/BV/BV_t19/ (one commit, subject "caltech -> orion: BV BV_t19 stopped").
3. Delete BV_t19's own rst and bins only.
4. After that, the ONLY Caltech work is the B_br build + smoke gate and the 6 C1R arms (NOTE hk / ccaf20cbd). Start nothing else.

BV_bsg and BV_bsgtau finals have been received, thanks.
