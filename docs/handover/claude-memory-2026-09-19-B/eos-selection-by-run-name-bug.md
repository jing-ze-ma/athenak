---
name: eos-selection-by-run-name-bug
description: "The soleos analysis scripts picked the EOS by guessing from the run name (or by argument position), which silently analysed tabulated runs as an ideal gas; fixed, and how to spot it"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 6d2fc30e-834e-473c-b800-b2c8b589b054
  modified: 2026-08-15T07:02:17.712Z
---

Found and fixed 2026-08-15 in `athenak/tools/solar_convection/` (commits `ee11563f`,
`c4de025b`). The analysis scripts chose between `IdealEOS` and `TableEOS` by GUESSING:

- `acoustic.py`: `'table' in run or 'lid' in run`
- `slices.py` / `tau_surface.py`: `run.startswith('table')`
- `compare.py`: by ARGUMENT POSITION — first run ideal, second run table
- `lidtest.py`: hardcoded `TableEOS` for every run

A run directory named `sponge08` matches none of the name patterns, so a tabulated-EOS
run gets analysed with ideal-gas formulas and NO error. Damage found: the tau=2/3 maps
for the three `lidtest/lid*` runs came out at **8700-10600 K** instead of ~5000-5600 K
— the mu 0.602 vs 1.26 ratio of ~2.1 in `T = p/(rho R)` — and because `get_kapr(rho,T)`
then places the tau surface at the wrong height, the v_z panels were wrong too. Nothing
else was affected, and no CONCLUSION changed (the lid-test result came from `lidtest.py`,
which used the tabulated EOS throughout). Re-rendered correctly.

**Why:** an EOS is a property of the RUN, recorded in its input file; deriving it from a
directory name or an argument slot is a guess that fails silently and produces
plausible-looking numbers rather than an error.

**How to apply:** `eoslib.eos_for_run(run_dir, grid)` parses `<hydro>/eos` and
`general_eos` from the run's own `sun.athinput` and RAISES if the file is missing or names
an EOS it cannot evaluate. Use it everywhere; never add another name-based branch.

**The tell I missed:** the colour limits were 4600-6000 K, so every temperature map in
those sets was uniformly saturated, and the per-frame log lines printed `T[8835, 10425] K`
right in front of me. When a script prints a physical range, read it against what the run
should give — a photosphere is ~5800 K, not ~10000 K. See [[solar-convection-general-eos]].

**The audit that found the full scope**, worth repeating after any similar fix: enumerate
every output directory, map it back to its run, and evaluate the OLD selection rule against
`eos_for_run` for each — do not reason from memory about which scripts were run on what.
