# cs collapse: analysis scripts and results

The scripts that closed the cubed-sphere vertex dt collapse
(see `docs/handover/HANDOVER-2026-09-11.md` and the memory note
`cs-vertex-dt-collapse-0907-defaults.md`), copied out of the working directories
so they survive the machine.

**All `bench/...` paths below are VIPER-ONLY** -- they live under
`/viper/u2/jinma/ATHENAK/bench/` and are not in git and not visible from any other
machine. The scripts hard-code those absolute paths; on another machine they are a
record of the method, not runnable as-is. Several also `sys.path.insert` on
`athenak/vis/python` (for `bin_convert`), on `athenak/docs/handover/scripts`
(for `read_cyclediag.py`, which IS in git, one directory up) and on
`bench/cs_ens/analysis` (for `dhjcs.py` and `csgeom.py`, which are NOT).

## Cold-map geography -- from `bench/cs_ens/analysis/coldmap/`

| script | what it does | reads |
| --- | --- | --- |
| `coldmap.py` | The main geography analysis: where are the cold columns (>30% of cells below 1500 K) and the drained columns (rho(i=52) < 1e-9) -- at the cube vertices, at the panel edges, or in the panel interiors? Compares ctl2/s01 (explicit RT, dies rot 5.71) against si/s01 (semi-implicit, survives) at rot 1-5, and writes the per-panel cold fractions to `coldfrac.npz`. Produced the vertex/edge/interior enrichment numbers and the 418/322/102/770-vs-59/0/0/0 column counts. | `bench/cs_ens/{ctl2,si}/s01/bin/dhj.hydro_w.*.bin`, the EOS table `bench/cs_ens/analysis/eos_table.txt` |
| `coldfig.py` | Plots the cold-fraction maps written by `coldmap.py`: the 6 cubed-sphere panels side by side with the cube vertices marked, for ctl2 and si. | `bench/cs_ens/analysis/coldmap/coldfrac.npz` |

## si/s05 post-mortem -- from `bench/cs_ens/analysis/si_s05/`

The one semi-implicit death (rot 5.34).

| script | what it does | reads |
| --- | --- | --- |
| `coldmap_s05.py` | `coldmap.py`'s method applied to si/s05 at rot 1-5, plus the meshblock-gid mapping needed to point `problem/diag_gid` at the abort column. Showed s05 had NO cold patch and NO cold bias but 65 drained columns in the antistellar panel interior. | `bench/cs_ens/si/s05/bin/*.bin`, `bench/cs_ens/si/s01/bin/*.bin` |
| `profiles_s05.py` | Radial profiles (rho, T, v_r) of the deepest drained columns of si/s05 at rot 1-5, with the same columns in si/s01 for reference. Gave the gid 0 (J,K)=(14,11) plateau 2-5e-10 over i 52-68 at 2340 K. | same |
| `root_s05.py` | Root-cell analysis for a si death, following `root2_task1..4`: finds the first cell above the deep interior to exceed 5000 K, then dumps its 200-cycle history, the drain history at i and i+-4, and the per-cycle block counts. Takes `<cyclediag_dir> <k> <j> <abort_cycle>` on the command line. | a `cyclediag/` directory written by `problem/diag_gid` |
| `root_fast_shim.py` | Local copy of `root2_fast.py` so `root_s05.py` can import it from its own directory. Identical content. | -- |
| `si_s05_root2_fast.py` | The si_s05 directory's own copy of `root2_fast.py` (same file, kept under a prefixed name because of the name collision with the twin2 copy). Byte-identical to `root2_fast.py` and `root_fast_shim.py`. | -- |

## Twin cyclediag analysis -- from `bench/cs_ens/ctl2_s01_twin2/analysis/`

All three twin arms restart from the same rot-5.700 seed
(`bench/cs_ens/ctl2_s01_rst570/rst/dhj.00006.rst`): `exp_g3` (explicit, `diag_gid=3`, the
block the explicit collapse lands in), `exp_g6` (explicit, `diag_gid=6`), `si_g6`
(semi-implicit, `diag_gid=6`). Every script here reads per-cycle
`<arm>/cyclediag/dhj.cyclediag.<cycle>.dat` files written by the instrument (1d988879);
`root2_task*` are run from inside `bench/cs_ens/ctl2_s01_twin2/`.

| script | what it does | reads |
| --- | --- | --- |
| `twin2_analysis.py` | The main twin analysis: the RT flux components (`rt_Ft`, `rt_Fb`, `rt_Qs`, `rt_Em`) in the collapsing vertex column across all three arms, purely descriptive (no modelling, no fitting). Established that `rt_Qs == 0` on the nightside, so the source is purely `-(Ft-Fb)/dx`. | `ctl2_s01_twin2/{exp_g3,exp_g6,si_g6}/cyclediag/` |
| `root2_fast.py` | Not an analysis: a fast partial reader for the cyclediag format that seeks to and reads only the selected variables' `(k,j,:)` rows (`read_rows`) or selected full arrays (`read_full`), instead of loading all 23 arrays. Imported by most of the `root2_task*` scripts. | -- |
| `root2_scan.py` | Scans an arm's cyclediag files for the first cycle at which any cell in column `(k,j)` exceeds a temperature threshold. Usage: `root2_scan.py <arm> <k> <j> <thr>`. | `<arm>/cyclediag/` |
| `root2_scan2.py` | Same scan, but ignores the cells that were already above the threshold in the first file, i.e. finds the first NEW crossing -- which is how the root cell was located. | `<arm>/cyclediag/` |
| `root2_task1.py` | The si_g6 root-cell history table (column k=17, j=2, i=68): dt, T, rho, v_r, e_int, delta_u, `rt_src`, `rt_de`, `rt_Fb`, `rt_Ft`, `rt_Em`, `rt_clip` and the neighbours' T, every 5 cycles over the drain and every cycle near the abort. Result: `results/root2_task1.txt`. | `si_g6/cyclediag/` |
| `root2_task1_scan.py` | Finds the first NEW cell above 5000 K in the si_g6 column (17,2) -- the si ignition cycle. | `si_g6/cyclediag/` |
| `root2_task1b.py` | Fills the 138855-138858 cycle gap in the si_g6 table and adds the earlier i=68 history every 20 cycles from 138354. | `si_g6/cyclediag/` |
| `root2_task2.py` | The exp_g3 column (k=16, j=17) cycle-by-cycle over i=48..58 for cycles 138769-138800: T, rho, v_r, `rt_de` and `delta_u` per cycle -- the table that shows delta_u > 0 while `rt_de` is a large sink, i.e. HYDRO supplies the energy at the root. Result: `results/root2_task2.txt`. | `exp_g3/cyclediag/` |
| `root2_task2b.py` | The same data transposed into one matrix per variable (cycle down, i across), easier to read. Result: `results/root2_task2_matrix.txt`. | `exp_g3/cyclediag/` |
| `root2_task3.py` | exp_g3 whole-block counters per cycle: number of cells at the T floor (200 K), number below rho 1e-9, number with `rt_clip` set, the summed clip magnitude and the block max T -- the cascade timeline. Result: `results/root2_task3.txt`. | `exp_g3/cyclediag/` |
| `root2_task4.py` | The same block counters for exp_g3 AND si_g6 every 20 cycles, for the explicit-vs-semi-implicit comparison. Result: `results/root2_task4.txt`. | `exp_g3/cyclediag/`, `si_g6/cyclediag/` |

## `results/`

Text output of the runs above, kept because the data directories they read are
viper-only and 4-31 GB each.

| file | from |
| --- | --- |
| `coldmap.txt` | `coldmap.py` -- the cold/drained column geography, ctl2 vs si, rot 1-5 |
| `panels.txt` | the per-panel summary that goes with it |
| `coldmap_s05.txt` | `coldmap_s05.py` -- si/s05 cold and drained maps + the gid mapping |
| `profiles_s05.txt` | `profiles_s05.py` -- the deepest drained columns of si/s05 |
| `comparison.txt`, `precursor.txt`, `log_extract.txt`, `replay_dt_trace.txt` | the rest of the si/s05 post-mortem: si/s05 vs si/s01, the hst precursor (radial KE x3.9, floor counters x3 in the last 0.1 rot), the log extract around the abort, and the dt trace of the non-reproducing replay |
| `root2_task1.txt`, `root2_task2.txt`, `root2_task2_matrix.txt`, `root2_task3.txt`, `root2_task4.txt` | the `root2_task*` scripts above |

`bench/cs_ens/ctl2_s01_twin2/analysis/run_out.txt` (the full 574 KB console log of
`twin2_analysis.py`) was **NOT copied**: it is over the 500 KB limit for this bundle.
The `.npz` and `.png` products in those directories were not copied either.

## flake8

Checked with `python3 -m flake8 --max-line-length=90` (the repo's Python limit, see
`setup.cfg`). Only whitespace, blank-line, continuation-indent, statement-splitting
(`import a, b` and `a=1; b=2` one-liners), line-length and unused-import fixes were
applied; no logic, no numbers, no name, no format-string content and no paths changed
(verified by comparing the token stream of every file against the pre-edit copy, with
adjacent string tokens merged so that implicit-concatenation wraps are proven
byte-identical). That took the count from **132 to 33**. `twin2_analysis.py` is fully
clean; the remainder is left deliberately, and these scripts are not part of the linted
tree anyway (`tst/run_test_suite.py --style` covers `src/`, `vis/` and `tst/`):

- **26x E402** "module level import not at top of file" -- every one of these scripts must
  `sys.path.insert` on `athenak/vis/python`, `athenak/docs/handover/scripts` and
  `bench/cs_ens/analysis` BEFORE importing `bin_convert`, `dhjcs`, `csgeom`,
  `read_cyclediag` or `root2_fast`. Not fixable without changing how the scripts find
  their modules.
- **3x E741** ambiguous variable name `I` (`root2_task1.py`, `root2_task1b.py`,
  `root_s05.py`) -- `I` is the root-cell radial index throughout the analysis; renaming it
  is a logic-adjacent edit and was out of scope.
- **3x F841** local `n1` assigned but never used, in the three identical copies of the
  cyclediag fast reader (`root2_fast.py`, `root_fast_shim.py`, `si_s05_root2_fast.py`).
  It comes from the tuple unpack `n1, n2, n3 = h['n1'], h['n2'], h['n3']` in `read_full`,
  where only `n2` and `n3` are used; breaking the unpack for one lint code was not worth
  the edit, and the three copies must stay byte-identical to each other.
- **1x E501** `coldmap.py:32` -- a 102-character one-line docstring. Wrapping it in
  parentheses would stop it being a docstring, so it stays.
