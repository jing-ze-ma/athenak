# T0 results — decision gate for `ck_implicit` on the tm sweep (2026-09-22)

Brief: `T0_BRIEF.md`; design: `DESIGN_tm.md` §4. No code changes. Run directory
`/viper/u2/jinma/ATHENAK/bench/impl_t0_0922/`.

## Verdict: **GO**, on the first criterion alone

At the production configuration (tm + spherical + beam, semi-implicit gas coupling) from the
settled prod3 restart, the sweep-to-gas gap `rt_desum` is **-4.7 % to -12.9 %**. That is
about 5-13x the 1 % threshold, and it holds steady after the restart transient. The second
criterion (T moving by more than 5 K between CFL 0.3 and 0.15) is **not decidable as
measured**, because trajectory divergence dominates it (see A). The GO does not depend on it.

## B. `rt_desum` and the clip count (job 11942457, apudev, 2 GPUs, 171 s wall)

* Binary: a copy of `bench/cs_mhd_prod4/athena`, md5 `f95130b22e6b7f9f57d86ca44da22058`
  (commit 395db5bc). Restart: `bench/cs_mhd_prod3/rst/dhj.00567.rst` (ncycle 4429831).
  Input: `bench/tm_prof_growth/prof_tm/deep_hot_jupiter.athinput` plus `rt_cell_report = true`,
  `rt_report_every = 20` and `rt_outer_verbose = true`, all added in the input's trailing
  `<problem>` block. The log confirms the pseudo-spherical beam, the spherical correlated-k
  form and the probe-free tm sweep (`ck_sweep_form = 1`).
* 2000 cycles (4429831 to 4431831). Hydro dt = **13.2-18.9 s** (the `dt=` lines in `run.log`).
  This is the restart state's CFL dt, not the ~5.3-5.9 s of prod4's cold start (the A run
  logs dt = 5.9 s at t = 5.2e4 s).
* **Not per MeshBlock, as the brief asked.** `rt_desum` (`two_stream_rt.hpp:6064-6069`) is a
  rank-wide sum over that rank's 12 MeshBlocks, printed once per RK stage by every rank with
  no rank label. Each report cycle therefore gives 4 lines, in 4 slots (2 ranks x 2 stages).
  Getting a per-MeshBlock figure would need a code change.
* `rel = sum(de dx)/sum(src dt dx) - 1`. The window is cycles 4430040-4431820, i.e. from
  200 cycles after the restart, 90 report cycles per slot:

| slot (line order within a cycle) | mean | min | max |
| --- | --- | --- | --- |
| 0 | -11.73 % | -12.87 % | -10.86 % |
| 1 | -9.84 % | -11.34 % | -9.00 % |
| 2 | -6.21 % | -6.83 % | -5.72 % |
| 3 | -5.16 % | -5.98 % | -4.71 % |
| all 360 lines | -8.23 % | | |

  The first report (cycle 4429840) gave -14.8 / -13.5 / -7.9 / -7.2 %. After the transient
  the gap settles at the values in the table and does not shrink.
* **`rt_de_max` clip:** `nclip = 0` at all 10 `rt_clip` census points (every 200 cycles,
  rank 0, both stages), with `efix_tot = resc_eq = resc_floor = 0`. The one-time
  `LimitRTSource` warning never printed. No rescued-cell `rt_cell_report` lines appeared.
  The gap therefore comes from the semi-implicit relaxation itself, not from the limiter.
* Caveat: dt here is ~2.5x the cold-start production dt. If the gap scales roughly with
  dt, that would put it near -2 to -5 % at 5.3 s. That figure is inferred, not measured, and
  it is still above 1 %.

## A. T(p) between CFL 0.3 and 0.15 at equal time

**Not possible from the existing dumps.** The last dumps of `g_cfl03` and `g_cfl015` are at
t = 81525.9 s and t = 52206.9 s (wall-clock exits), and the only common time is t = 0.
**What I did instead (job 11942456, apudev, 389 s):** reran `g_cfl03` from scratch with the
same binary (md5 `a6dd8e0b...`) and the same input, with `time/tlim = 52206.9`
(`impl_t0_0922/A_cfl03_eqt`), then compared its final dump with `g_cfl015` 00001 at the same t.
T comes from inverting the EOS table (`bench/cs_ens/analysis/eos_table.txt`). The table's
header says `mh 0`; whether that matches `eos_metal_condensation = true` was not checked.
Both arms go through the same inversion. Day and night are split by the sign of
cos(lat) cos(lon_ss). Script: `impl_t0_0922/analysis/tp_diff.py`; output:
`analysis/A_tp_diff.txt`.

| p [bar] | day d<T> | day max/rms \|dT\| | night d<T> | night max/rms \|dT\| | global max/rms |
| --- | --- | --- | --- | --- | --- |
| 1e-6 | -0.2 | 6.5 / 0.7 | +5.6 | 2434 / 222 | 2434 / 157 |
| 1e-5 | +0.7 | 81 / 6.4 | +8.1 | 1344 / 104 | 1344 / 74 |
| 1e-4 | +0.7 | 56 / 4.5 | -7.3 | 596 / 98 | 596 / 69 |
| 3e-4 | +1.1 | 107 / 8.6 | -14.3 | 1422 / 124 | 1422 / 88 |
| 1e-3 | +0.9 | 363 / 21 | -6.9 | 694 / 89 | 694 / 65 |
| 1e-2 | +1.2 | 292 / 25 | +2.5 | 557 / 102 | 557 / 74 |

(K, CFL 0.15 minus CFL 0.3; the intermediate levels are in the output file.) On the day
side the horizontal-mean shift is at most 1.2 K. The night side differs by rms ~100 K per
column and by up to 14 K in its mean. This is 0.17 rotations after a cold start, while the
night-side flow is still developing. Without a same-dt control (for example CFL 0.29
against 0.3) there is no way to separate dt error from chaotic divergence, so this
criterion is **inconclusive**, not passed.
