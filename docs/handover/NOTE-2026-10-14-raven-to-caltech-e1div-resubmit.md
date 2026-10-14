# NOTE Raven -> Caltech: please RESUBMIT the E1 divergence arms (USER, 10-14)

The user asks Caltech to run the E1 restart arms again (R0/R1 were removed when the queue was cleared for bsgLIc).
Priority order (all from rst 00022, t 1.1e7 s, dirs e1div_rst/ as before; each 1 node x 2 H200):
1. **D1** (TASK-2026-10-14-caltech-e1div-dtretry.md, part 3): H200 build of fin4-dtretry @ 5b7e1bc5, restart from rst 00022 with
   `time/dt_retry=true`, tlim 1.30e7. THE key run: does dt-retry get E1 through the infall? (Run its short gate F22/F88 first only
   if it does not delay D1 much; Raven already passed the full-size gate on A100.)
2. **R1**: dad00933 binary, `time/cfl_number=0.15`, tlim 1.30e7 (does half dt alone survive?).
3. R0 (baseline divergence repro) is optional -- skip if GPUs are short.
These may run alongside bsgLIc and the He giant (user decision); report start + results (t reached, DT_RETRY lines, NC/DIVERGED,
Picard, profiles via e1div_prof.py at 1.15e7/1.2e7/1.25e7/1.3e7) in NOTE-2026-10-14-caltech-e1div-rst.md.
