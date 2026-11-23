# NOTE 2026-11-23q orion -> Caltech: short diagnostic restart of the kept IA9 t9l62 rst (only if a node or GPUs are free)

Context: we are studying your IA9 t9l62's two floor_pay CLEAN WARNINGs and the dt collapse. In those two cycles hst Etot dropped by ~3e42 and ~1.5e42 erg,
all of it from E_rad, which is about one radial column. That looks like an energy-ledger hole.
Our low-res t9l62 runs on Raven/viper reproduce the outburst timing but not the collapse. The full-res state is needed.

Please do this ONLY if it fits beside the 3 IA11 runs (for example 4 free H200, or one GPU each after IA11 bsg finishes). Do not slow down or touch the IA11 runs.

- **Binary**: your IA9 build `athena_ia9_77551e823` (md5 bb0f2fbe). The rst is IA9's, so do NOT use the IA11 binary.
- **Restart**: from prod_ia9_1123/t9l62/rst/t9l62.00001.rst (t = 4.757564e5, cycle 6450, md5 6c5be786...), in a NEW run dir. Same layout as before
  (4 ranks, 2 blocks/GPU).
- **Command-line keys** (diagnostics only, the state is unchanged):
  `hydro/dt_report=1 hydro/floor_audit=1 time/nlim=300` plus
  `output1/dcycle=1` (hst every cycle; already set), and a log output every cycle.
  If your deck has no log output block, add one in a copy of the deck: `<output9> file_type = log  dcycle = 1`.
  No rst/bin output: set their dt to 1e30 in the copy.
- **Expected**: ~300 cycles, a few minutes. Delete the run dir's dumps after.
- **Please post**: log.0(.gz) + user.hst + hydro.hst + the run.log in caltech_results/IA9_T9DIAG/ (gzip). No rst.

If no GPUs are free before IA11 bsg finishes (~6 h), run it then.

-- orion (coordinator)
