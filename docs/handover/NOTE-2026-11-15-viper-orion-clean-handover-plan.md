# NOTE 2026-11-15 viper -> Orion: plan — the clean radiation scheme moves to Orion

The user decided that the Orion session takes over the whole clean-scheme development (code AND CPU tests) on
branch `fork/clean-cc`. All clean-scheme work is CPU work (1-D waves, single-block IDORT decks, Python
prototypes), and Orion can run the edit/build/test loop directly, which viper cannot (viper is GPU-only for jobs,
its login node is capped at ~12 cores).

What this means for the two TASKs pushed earlier today:
- `TASK-2026-11-15-orion-cc_idort_head.md` (24 IDORT runs, code dc67dba7b = current clean-cc head): still valid.
  Run it as your first job if you like — it is the baseline of the scheme you take over.
- `TASK-2026-11-15-orion-cc_fixedD.md` (156 wave runs, code 87e045577): built on an older commit. Rerun it on the
  current head dc67dba7b instead, and use `clean_tol_nl >= 1e-12` for the STAGGERED layout (its Newton loop stalls
  at 1e-13). No need to report it as a separate NOTE; it becomes part of your own work.

Coming in the next commit (within the hour): `HANDOVER-2026-11-15-clean-scheme-viper-to-orion.md` plus a package
`docs/handover/handover_1115_clean/` (progress logs, 1-D prototype, scorers with exact references, wave drivers,
run lists, the user's rules and the open problems in priority order). Please wait for it before changing code.

Ownership from now on:
- Orion owns `fork/clean-cc` (and may create branches off it). viper makes NO further commits to clean-cc.
- viper keeps: productions and dhj runs (GPU), the restyle branch `restyle-1114` (merged into `clean-1114` per
  module; its CPU bitwise checks come to you as separate A/B TASKs), GPU validation of the clean scheme on request.
- Report results as NOTE-*-orion-viper-*.md on bsg-files-1009 as usual; viper polls.
