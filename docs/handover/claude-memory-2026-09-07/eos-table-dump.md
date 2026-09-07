---
name: eos-table-dump
description: New <hydro>/eos_table_dump option writes the generated general-EOS table to ASCII so post-processing can recover p and T
metadata:
  type: project
---

**Committed as f882159f** (2026-08-25, eighth session) on `polar-average-perf`, on top of
3429f59f; documented in `docs/general_eos.md`. Pushed to the fork 2026-08-26. Adds `<block>/eos_table_dump = <file>`: after `BuildEOSTable` samples the grid,
rank 0 writes an ASCII dump of the two surfaces at every node --
`log10(e/rho)` and `log10(p/rho)`, both cgs -- with the grid and composition as `#` comment
lines so the body loads with a bare `np.loadtxt`.

**Why it exists:** binary `hydro_w` dumps carry only `dens, velx, vely, velz, eint`. Under
the tabulated EOS there is no constant `gamma-1`, so pressure and temperature cannot be
recovered from a dump without the EOS. There is no pressure/temperature derived output
variable in `src/utils/derived_vars.cpp` either. With the dump, post-processing inverts
`e(rho,T)` for `T` then evaluates `p(rho,T)` -- the run's own EOS, not a fit.

**How to use it:** the option must be added to the athinput; `ModifyFromCmdline` can only
override parameters that already exist (`parameter_input.cpp:390` fatals otherwise). A
serial CPU build suffices (the table is host-side) and takes ~1 min:
`cmake -S . -B <dir> -D PROBLEM=deep_hot_jupiter_rt -D CMAKE_BUILD_TYPE=Release`.
Run with `time/tlim=0` and a shrunk mesh -- but `nx3/meshblock_nx3` must be EVEN with
polar boundaries (`mesh.cpp:334`), so nx3=32/mb=16 works and nx3=16/mb=16 does not.
The reader + inverter live in the session scratchpad as `eos.py`.

**POST-PROCESSING TRAP: see [[eos-inversion-nan-trap]].** The NaN rows below wreck a naive
root find, silently and plausibly. Validate any inversion against `problem/ck_dump_file`.

**Finding that fell out of it:** the table is **NaN for log10 T < 1.85 (T < 71 K)**, over
the whole density range. That is the mechanism behind the already-measured
"`tfloor_kelvin = 50` crashes the run" in [[dhj-floors-for-1e-6-bar]] -- 50 K sits inside
the NaN region. The production floor of 200 K is clear of it, and zero cells of the
[[ck-hydro-long-run]] snapshot touch a NaN node.

**On the style gate:** `tst/test_suite/style/check_athena_cpp_style.sh` reports ~2363
cpplint errors on the UNTOUCHED tree, so a non-zero exit from it says nothing about your
change. Diff the per-file counts against `git show HEAD:<file>` instead. The three it flags
here are cpplint misreading `for (int j=0; j<ny; ++j)` and `h_tbl(j,i,ITE)` -- the file's
own existing idiom, which it already flags 17 times.
