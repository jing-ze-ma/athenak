# `tests_gpu_ab_0922` — analysis of the 2026-09-22 apudev smoke A/B arms

Read-only analysis of `bench/ck_sph_ab/{prodbin,sph,sphbeam,offfix}`, `bench/cs_noang/noang`
(comparison 1) and `bench/cs_recon_ab/{plm_ctl,ppmx,wenoz}` (comparison 2). `offfix`
(`bench/ck_sph_ab/offfix`, binary `8e41808b`, `ck_spherical=false` on the repaired
plane-parallel path, same rot-283 restart as the other comparison-1 arms) was added
2026-09-22 to isolate the effect of `ck_spherical` from the effect of the binary; see
`bench/ck_sph_ab/README.md` §5 and the new section below. Nothing in `bench/` was written to.

Scripts (all `python3 <name>.py` from this directory):

| file | what it does | output |
| --- | --- | --- |
| `common.py` | readers (hst / `dhj.log` / stdout / `rt_desum`), the cubed-sphere metric reproduced from `coordinates.cpp`, and an EOS-table inverter that also returns `Gamma_1` | — |
| `logs.py` | Q1 (dt), Q3 (floors), Q4 (energetics), Q5 (cost) | `logs_tables.md`, `fig_dt.png`, `fig_floors.png`, `fig_floors_recon.png`, `fig_energy.png` |
| `cfl.py` | Q1: per-cell CFL limit rebuilt exactly as `mhd_newdt.cpp` does it; 10 most limiting cells; floor census | `cfl_tables.md` |
| `tp.py` | Q2: T(p) day / terminator / night, and T(r) | `tp_tables.md`, `fig_tp.png`, `fig_Tr.png` |
| `noang_diff.py` | Q6 | `noang_diff.md` |
| `vertex.py` | comparison 2: vertex band vs panel interior | `vertex_tables.md`, `fig_vertex.png` |

The EOS table used for T, p and `Gamma_1` is the run's own dump
`bench/cs_ens/analysis/eos_table.txt` (281x451, X=0.7381, Y=0.2485, H2 + ionization +
metal ionization on) — its header matches this run's startup print exactly.

---

## Caveat that governs all of comparison 1 — read first

`prodbin` is **not** the HEAD binary with the switch off. It is the production binary
(md5 `6d77fe19…`, identical to `cs_mhd_prod3/athena`, commit lineage `18b9c563`/`bc972cc8`),
while `sph`, `sphbeam` and `noang` are commit `916dc953`. They differ in far more than
`ck_spherical`: HEAD's `dhj.log` has three extra counter columns, and HEAD warns
*"restart file has no general-EOS temperature cache … this restart is not bitwise"*,
which `prodbin` does not.

The intended HEAD/plane-parallel control (`ck_sph_ab/off`, job 11931202) is invalid
(regressed hybrid), **but it is still the HEAD binary on the same restart**, and it is
decisive for one thing: its magnetic energy at the first history row is 4.0507e33, versus
4.0560e33 in `sph` and 4.9767e33 in `prodbin`. The ~19 % magnetic-energy step and the
ensuing decay therefore belong to the HEAD binary / restart, **not** to `ck_spherical`.
By the same token the dt result below cannot be attributed to `ck_spherical` alone —
`off`'s own dt collapsed to ~1 s for an unrelated reason, so there is no usable control.

Clean one-line comparisons in this set: `sph` vs `sphbeam`, `sph` vs `noang`, and now
`sph` vs `offfix` (see the new section below) — the last is the actual `ck_spherical`
on/off pair on (nearly) the same binary lineage, and `offfix` vs `prodbin` isolates the
364-commit binary gap on top of it.

---

## Q1 — what sets dt

`fig_dt.png` (left), `cfl_tables.md`.

Per-cell reconstruction of `mhd_newdt.cpp` (dt_i = dx_i/(|v_i|+cf_i), cf_2,3 divided by
`sin_cell`, dx from the polynomial radial stretch and the exact gnomonic arc lengths,
`Gamma_1` and p from the run's EOS table) gives, at each arm's last dump:

| arm | hst dt [s] | reconstructed MHD CFL [s] | radial / xi / eta |
| --- | --- | --- | --- |
| prodbin | 19.82 | **19.81** | 19.8 / 82.3 / 88.8 |
| sph | 12.53 | 17.36 | 17.4 / 80.2 / 83.5 |
| sphbeam | 12.85 | 18.00 | 18.0 / 84.5 / 89.2 |
| noang | 12.53 | 17.85 | 17.8 / 80.2 / 83.5 |
| offfix | 19.92 | **19.93** | 19.9 / 81.1 / 88.3 |

**In `prodbin` the MHD CFL sets dt** (19.81 reconstructed vs 19.82 reported — a 0.05 %
match, which also validates the metric and the EOS inversion). Its limiting cells are
real gas: r/Rp 1.25, equator, night-ish (lon +150…+165), p 1e-2…1e-6 bar, T 1900–2700 K,
|v| 6–12 km/s, none on `dfloor`. The limit is radial in every one of the top 10.

**`offfix` reproduces this mechanism, not the HEAD-arm one.** Reconstructed CFL 19.93 s
vs reported 19.92 s (0.05 % match, same as `prodbin`). Its top-10 cells sit at r/Rp
1.25–1.40, lon +150…+167 (night side by the same convention), p 3e-4…1.4e-2 bar,
T 1850–2725 K, none on `dfloor` — the same real-gas mechanism as `prodbin`, at slightly
larger radius. So the Ohmic-cap collapse to ~12.5–12.9 s seen in `sph`/`sphbeam`/`noang`
is **not** a generic property of "a HEAD-lineage binary on this restart" — `offfix` is
built from the current `rt-integration` tip and still runs at the `prodbin` dt. It is
specifically tied to `ck_spherical = true` cooling the r/Rp 1.23–1.25 shell enough to
saturate `max_eta` there (see the new section below).

**In the three HEAD arms the MHD CFL is *not* what limits dt.** Their MHD limit is 17–18 s,
only ~10 % below `prodbin`, yet they run at 12.53/12.85 s. The reported dt is also
*quantised*: over 45 history rows `sph` takes the value 12.8523 s 23 times and 12.5323 s
17 times, and `prodbin` never repeats a value.

Those numbers are the **Ohmic timestep with `eta` pinned at the `max_eta = 1e13` cap**:

```
dt = cfl * (1/6) dx1(i)^2 / max_eta          (resistivity.cpp, NewTimeStepGeneralResist)
  i = 42, dx1 = 5.06997e7  ->  12.8523 s     (reported: 12.8523)
  i = 41, dx1 = 5.00646e7  ->  12.5323 s     (reported: 12.5323)
  i = 40, dx1 = 4.94653e7  ->  12.2341 s     (reported: 12.2341)
```

Six-significant-figure agreement on three separate values. The `fac = 1/6 dx^2/coef`
form is unique to resistivity (conduction's is `dx^2 rho c_v / kappa`, and its radial
branch is skipped here because `rad_implicit_x1 = true`). So dt in the HEAD arms is set
by cells at **r/Rp = 1.230–1.241** (the H2-dissociation / recombination shell, the finest
part of the radial grid) whose resistivity has saturated the cap.

Why there: `fig_Tr.png` shows the night-side horizontal mean at r/Rp 1.236–1.246 is
**330–450 K colder** in the HEAD arms than in `prodbin`. `ohmic_resistivity = eos` makes
eta rise as the electron fraction falls, so a colder shell saturates `max_eta` over a
much larger volume; in `prodbin` the saturated region must lie above r/Rp ≈ 1.35
(dx1 ≥ 6.3e7) or the Ohmic limit would already have beaten 19.8 s. **This last step is an
inference**: `eta` and x_e are not in the dumps, so the cap could not be mapped directly.

Not a restart transient: dt falls to the plateau inside the first ~0.03 rot and then sits
flat for the remaining 0.44 rot with no drift and no recovery (`fig_dt.png`).

Answer to the literal question ("night-side top slab on dfloor? day side? a vertex?"):
**none of those** — neither in `prodbin` (where it is a real night-side equatorial cell at
r/Rp 1.25, not floored, not at a vertex) nor in the HEAD arms (where a diffusive limit at
r/Rp 1.23 beats the CFL entirely).

## Q2 — T(p) by region

`fig_tp.png`, `tp_tables.md`. Horizontal means, day = |lon_ss| < 30°, terminator
75–105°, night > 150°. Each arm at its own last dump (times differ by up to 0.05 rot).

Night side, T [K]:

| p [bar] | prodbin | sph | sphbeam | noang |
| --- | --- | --- | --- | --- |
| 1e-6 | 1951 | 1427 | 1906 | 1434 |
| 1e-5 | 2109 | 1561 | 2066 | 1556 |
| 1e-4 | 2111 | 1629 | 1949 | 1624 |
| 1e-3 | 2028 | 1600 | 1616 | 1593 |
| 1e-2 | 2127 | 1768 | 1789 | 1767 |
| 1e-1 | 2661 | 2627 | 2631 | 2628 |
| 1e+0 | 3334 | 3329 | 3331 | 3329 |

Day-night contrast [K]: 2370 (prodbin) / 2819 (sph) / 2373 (sphbeam) / 2812 (noang) at
1e-6 bar, falling to 9–13 K at 1 bar. Differences > 50 K (full list in `tp_tables.md`):

- **`sph` vs `prodbin`**: night −428…−553 K at 1e-6…1e-3 bar, terminator −104…−337 K,
  day −76…−305 K. Everything at and below 1e-1 bar agrees to < 35 K in every region.
- **`sphbeam` vs `sph`** (the clean one-line comparison): the pseudo-spherical beam
  **warms the terminator by +721 K at 1e-6 bar**, +414 K at 1e-5, +270 K at 1e-4, and
  brings the night side back to within 45 K of `prodbin` at 1e-6/1e-5 bar. This is the
  predicted twilight illumination, and it is the largest single effect in the whole set.
  Per `916dc953` the new beam deposits 1.106x the exact absorbed power, so a warmer
  `sphbeam` is not by itself validation.
- **`noang` vs `sph`**: ≤ 7 K everywhere.

## Q3 — floors

`fig_floors.png`, `logs_tables.md`, `cfl_tables.md`. `fofc` is identically 0 (not set),
`eos_fail` is 0 in every arm, `eos_vceil` 0.

| arm | dfloor/cycle | efloor/cycle | tfloor/cycle |
| --- | --- | --- | --- |
| prodbin | 1.41e5 | 1.58e2 | 1.25e2 |
| sph | 1.59e5 | 2.03e4 | 6.16e3 |
| sphbeam | 1.05e5 | 3.62e3 | 7.93e2 |
| noang | 1.59e5 | 2.02e4 | 6.16e3 |
| offfix | 1.44e5 | 1.50e3 | 5.39e1 |

Density flooring is comparable (`sphbeam` is 25 % *lower* than `prodbin`), but energy and
temperature flooring are **130x and 49x higher** in `sph`/`noang` than in `prodbin`, and
`sphbeam` takes about 5/8 of that back. This is the same cold night-side top seen in Q2 —
the spherical thermal form cools the thin gas onto `tfloor_kelvin = 200` far more often.
The claim that the spherical form floors *less* in the upper atmosphere is **not
supported here**; the opposite is seen (with the binary caveat above).

Location, last dump: floored cells (`rho <= dfloor` to 1e-4) are 5.00 % of the mesh in
`prodbin`, 5.79 % in `sph`, 5.82 % in `noang`, 3.23 % in `sphbeam`, **5.00 % in `offfix`**
— matching `prodbin` almost exactly. They are almost purely night-side (day-side fraction
0.003–0.014, 0.002 in `offfix`), spread over all longitudes with mean |lat| 28–31° (30.2°
in `offfix`), and the **inner edge of the floored slab moves down** from r/Rp 1.508
(prodbin) to 1.404 (noang) / 1.381 (sph), while `sphbeam` keeps it at 1.499 and `offfix`
sits at **1.588** — the deepest of any arm, and if anything slightly less floored than
`prodbin` in depth even though the two per-cycle rates agree. `p <= pfloor` occurs in 242
(sph) / 230 (noang) / 1 (sphbeam) / 0 (prodbin) / **0 (offfix)** cells. On both energy and
density-floor location, `offfix` tracks `prodbin`, not the `ck_spherical=true` arms —
consistent with Q1: the switch, not the binary, is what floors the upper atmosphere more.

## Q4 — energetics and `rt_desum`

`fig_energy.png`. Over the run: `prodbin` |dE/E| and |dM/M| are below the 6-digit history
precision (< 2e-6) over a full rotation; `sph`/`noang` lose 2.56e-4 of total energy and
2.8e-5 of mass over 0.47 rot, `sphbeam` 2.03e-4 over 0.35 rot. `offfix` loses **8.85e-5**
of total energy and 1.13e-5 of mass over 0.84 rot — between `prodbin` and the
`ck_spherical=true` arms, and much closer to `prodbin` per rotation (see below). Magnetic
energy: 4.98→4.71e33 (−5.3 % over 0.97 rot, i.e. ≈ −5.5 %/rot) in `prodbin`, 4.06→2.92e33
(−28.0 % over 0.44 rot, ≈ −64 %/rot) in `sph`/`noang`, 4.06→2.95e33 in `sphbeam`, and
**4.06→2.85e33 in `offfix` (−29.8 % over 0.84 rot, ≈ −35 %/rot)**. Kinetic energy rises
~8 % in the `ck_spherical=true` arms, falls ~6 % in `prodbin`, and falls ~2 % in `offfix`
(1.881→1.838e34).

**The expectation that the spherical form is the more conservative one is now
contradicted, not just unconfirmed.** `offfix` (switch off, repaired binary) decays its
magnetic energy at ≈ −35 %/rot — faster than `prodbin`'s −5.5 %/rot but distinctly slower
than `sph`'s ≈ −64 %/rot on the same restart. So turning `ck_spherical` on makes the
secular ME drift *worse*, not better, on this binary; part of `prodbin`'s superior
conservation is a genuine binary/algorithm difference over the 364 commits (offfix still
loses 6x more ME per rotation than prodbin), but the further loss from `offfix` to `sph`
is attributable to the switch itself.

`rt_desum`: `prodbin`'s binary has no `rt_cell_report` parameter (`strings` count 0), so
the plane-parallel side of the 36–40 % vs 2–5 % comparison could not be measured on the
production binary. **`offfix` now supplies that plane-parallel-off baseline** (28 lines):
mean **−27.9 %**, range **−94.6 % … +47.5 %** — an order of magnitude larger in magnitude
and, unlike the "on" arms, sign-changing from cell-report to cell-report. The "on" arms
(`sph`/`sphbeam`/`noang`) still give a tight mean −7.8 % (range −12.4 %…−4.4 %). This
matches the *direction* of the 1-D column result ("off" swallows/mismatches far more than
"on") and puts a first real number on it in 3-D, though the "off" magnitude here (28 %,
highly variable) is not the clean 36–40 % of the 1-D column, and the sign flips the 1-D
result never showed — so this closes the direction of the open question in
`tests_ck_sph/README.md` §0.3 but not the magnitude or the sign behaviour.

## Q5 — cost

| arm | cycles | wall s | cycles/s | sim s per wall s | relative |
| --- | --- | --- | --- | --- | --- |
| prodbin | 15200 | 601 | 25.3 | 505 | 1.00 |
| sph | 10900 | 768 | 14.2 | 182 | 0.36 |
| sphbeam | 8700 | 764 | 11.4 | 147 | 0.29 |
| noang | 10800 | 763 | 14.2 | 181 | 0.36 |
| offfix | 13300 | 769 | 17.3 | 345 | 0.68 |

Per cycle the `ck_spherical=true` arms are 1.78x more expensive (`sphbeam` 2.2x); the
smaller dt costs a further 1.58x, so the throughput penalty is **2.8x** (`sphbeam` 3.4x).
`offfix` (switch off, current binary) is **1.46x** more expensive per cycle than
`prodbin` — this is the pure 364-commit binary overhead, with `offfix`'s dt essentially
unchanged from `prodbin`'s (19.98 vs 19.93 median) so throughput tracks per-cycle cost
1:1 (0.68x). Subtracting: of the 1.78x per-cycle cost in `sph`, **1.46x is the binary and
only ~1.22x is `ck_spherical` itself**; the remaining 2.8x/1.46x ≈ 1.9x extra throughput
loss on top of the binary overhead is essentially all the Ohmic-cap dt collapse (Q1), not
switch compute cost. `noang` is **not** cheaper than `sph` (14.16 vs 14.20 cycles/s,
0.3 %), so the transverse conduction planes are not a measurable cost here — contrary to
the expectation in `cs_noang/README.md`.

## Q6 — `noang` vs `sph`

`noang_diff.md`. There *is* a difference, and it appears immediately: 1-KE and 3-KE differ
by 2.1e-3 and 3.3e-4 relative at the very first history row (0.03 rot). But it stays at
that level — over 0.47 rot the maximum relative difference in any history column is
6.7e-3 (dt, i.e. one quantisation step), 2.1e-3 (1-KE), 1.6e-4 (1-ME); mass and total
energy are identical to printed precision at every row. Region-mean T agrees to ≤ 7 K at
every level, the dt plateau is the same two values, and the floor counters agree to 0.3 %.

Where the difference lives: the largest per-cell differences in the last dumps (which are
21 cycles apart, so this mixes physics with evolution) are all on the **night-side upper
atmosphere**, lat −7…+16°, lon −150…−176°, radial index 17–55 — exactly where the thin gas
is and where the transverse operator acts. rms relative differences: density 4.9e-5,
eint 3.7e-5, velocity components 6–18 %, field components 3.7–4.8 %.

**Conclusion: switching the horizontal radiative conduction off changes the night-side
upper atmosphere at the per-cell level but is dynamically and thermally negligible over
half a rotation, and costs nothing either way.** The expectation of a colder night side
and a sharper terminator is not seen at this duration.

---

## Isolation of `ck_spherical` (offfix)

`offfix` = `bench/ck_sph_ab/offfix`, HEAD binary `8e41808b` (post-`ed188f69`, the fix to
the `ck_spherical=false` plane-parallel path), `ck_spherical=false`, `ck_beam_sph=false`,
same rot-283 restart, same `sph/deep_hot_jupiter.athinput` with one line changed. It is
**not** the same literal binary as `sph`/`sphbeam`/`noang` (`916dc953`) or `prodbin`
(`27ca5b13`), but it is on the current `rt-integration` tip and is the first valid
plane-parallel control since the invalid `off` (`916dc953`) arm.

**(1) Does `offfix` reproduce `prodbin`'s dt plateau and night-side T? Yes, closely.**
dt: reconstructed MHD-CFL 19.93 s vs `prodbin`'s 19.81 s (+0.6 %); reported median 19.98 s
vs 19.93 s. Both are set by the same mechanism — real night-side/equatorial gas, none on
`dfloor` — not by the Ohmic cap. Night-side T at 1e-6 bar: **1960 K** (`offfix`) vs
**1951 K** (`prodbin`), +9 K (0.5 %); at 1e-5/1e-4/1e-3 bar the gap is +11/+12/+60 K (all
`tp_tables.md`). Floored-cell fraction (5.00 % vs 5.00 %) and per-cycle `eos_dfloor`
(1.44e5 vs 1.41e5) also match to a few percent.

**(2) `offfix` vs `sph` — the pure effect of the switch on one binary lineage.** Night
side warms by +533 K at 1e-6 bar (1960 vs 1427 K), +559 K at 1e-5, +494 K at 1e-4,
falling to +391/+34/+4 K at 1e-2/1e-1/1 bar — i.e. `ck_spherical=true` cools the whole
upper night side by several hundred K, exactly the pattern Q2 attributed to the switch
using the invalid `off` baseline, now confirmed with a valid one. dt: 19.9 s vs 12.5 s —
turning the switch on is what saturates the Ohmic cap at r/Rp 1.23–1.25 by cooling that
shell (Q1); `offfix` shows this shell never gets cold enough to cap when the switch is
off. Cost per cycle: `offfix` 0.0577 s/cycle vs `sph` 0.0704 s/cycle, i.e. **`sph` is only
~1.22x the per-cycle cost of `offfix`** once the binary term is removed — most of the
throughput penalty in `sph` is the dt collapse, not extra compute. Magnetic-energy decay:
≈ −35 %/rot (`offfix`) vs ≈ −64 %/rot (`sph`) — the switch roughly doubles the secular ME
loss rate on this binary. `rt_desum`: −27.9 % mean, sign-changing (`offfix`, switch off)
vs a tight −7.8 % (`sph`/`sphbeam`/`noang`, switch on) — the switch also tightens and
roughly quarters the sweep-to-gas mismatch, consistent with it being the intended fix.

**(3) `offfix` vs `prodbin` — the pure effect of the 364 commits, `ck_spherical` off in
both.** dt and night-side T agree to <1 %/60 K as above — no regression there. What does
differ: `offfix` loses energy 0.84 rot (dE/E = −8.85e-5, dM/M = −1.13e-5) where `prodbin`
loses none measurable in a full rotation, and its ME decays at ≈ −35 %/rot vs `prodbin`'s
≈ −5.5 %/rot — a real ~6x gap that has nothing to do with `ck_spherical`. Per-cycle floor
counters differ too: `eos_efloor`/cycle is 9.5x higher in `offfix` (1.50e3 vs 1.58e2)
though `eos_tfloor`/cycle is lower (53.9 vs 125) and `eos_dfloor`/cycle matches (1.44e5 vs
1.41e5). Per-cycle cost is 1.46x higher in `offfix`. **The first-row ME step noted as
undetermined in the caveat above is present in `offfix`** too: its first history row gives
ME = 4.056e33 erg, matching `sph`/`noang`/the old invalid `off` (4.051–4.056e33) and
**not** `prodbin` (4.977e33, +18.5 %). Since `offfix` and `sph` are different literal
binaries that share only the post-`916dc953` lineage and the same restart file/warning
("restart file has no general-EOS temperature cache … not bitwise"), this step is now
seen on two independent HEAD builds and is confirmed to be a **restart/binary-lineage
effect** (most likely the missing general-EOS temperature cache forcing a cold C2P start
on the first cycle), reproducible regardless of `ck_spherical`, and unrelated to the
switch.

---

## Comparison 2 — reconstruction (cold start, nghost = 3)

`fig_dt.png` (right), `fig_floors_recon.png`, `fig_vertex.png`, `vertex_tables.md`.
These arms also carry `ck_spherical = true` (the resubmitted jobs 11931395-7).

| arm | cycles | wall s | cycles/s | sim s per wall s | t reached | dt at the end |
| --- | --- | --- | --- | --- | --- | --- |
| plm_ctl | 10600 | 768 | 13.8 | 81.2 | 8.54e4 (0.28 rot) | 5.10 |
| ppmx | 10200 | 774 | 13.1 | 56.4 | 6.83e4 (0.22 rot) | 3.35 |
| wenoz | 10300 | 771 | 13.3 | 55.5 | 6.56e4 (0.22 rot) | 4.64 |

**Is the falling dt the normal spin-up transient? Yes.** The original production cold
start (`cs_mhd_prod3`, nghost = 2, plm, plane-parallel ck) went 28.98 s → 6.73 s by
t = 9.15e4 on exactly the same trajectory shape; `plm_ctl` is at 5.10 s at t = 8.54e4.
(The production history is written only every 3.05e4 s, so the dashed curve in
`fig_dt.png` has 4 points — the comparison is of endpoints, not of shape.)

Cost per cycle is within 5 % across the three arms; the whole difference in throughput is
dt. `ppmx` and `wenoz` are **1.44x more expensive per simulated second** than `plm`, and
neither is cheaper at any point after the first 0.1 rot.

Floors are a null: 6.20–6.28e4 `eos_dfloor`, 3.39–3.70e4 `eos_efloor`, 1.38–1.73e4
`eos_tfloor` per cycle, i.e. within 25 % across arms, `eos_fail` 0 everywhere, no NaN, no
`dt_min` collapse, no FATAL. (All three are ~5x above the `cs_mhd_prod3` cold start's
1.19e4/6.02e2/1.12e2 per cycle — but that run had nghost = 2 *and* plane-parallel ck, so
that gap is the ck switch, not the reconstruction.)

Energetics through the transient settle rather than run away: KE 2.44–2.68e33, ME
3.66–4.75e32, mass −0.55 %, total E −0.68…−0.69 % — the three arms agree to 10 %.

**Vertex band vs panel interior** (band = within 4 cells of a cube vertex in *both*
gnomonic angles, 384 cells; interior = the middle half of each panel, 1536 cells; arms
compared at the nearest available dumps, 0.215–0.250 rot, not the same instant):

| arm | rms residual wind, vertex / interior | grid-scale roughness, vertex / interior |
| --- | --- | --- |
| plm_ctl | 0.798 | 0.974 |
| ppmx | 0.810 | 0.944 |
| wenoz | 0.840 | 0.894 |

("residual" = departure from the zonal mean at the same latitude and radius; "roughness" =
the |1,−2,1| stencil on v_xi, v_eta inside each MeshBlock.)

**There is no vertex-localised spurious flow in the real atmosphere to reduce.** The
vertex bands are *quieter* than the panel interiors on both metrics in all three arms, at
every radius sampled except r/Rp ≈ 1.45 and 2.04 (`vertex_tables.md`). The higher-order
schemes do shift the vertex/interior roughness ratio down monotonically (0.97 → 0.94 →
0.89), which is the direction `tests_cs_regions` reported, but they *raise* the absolute
roughness in both bands (4.9e4 → 6.3e4 → 6.6e4 cm/s) and cost 1.44x. The same metrics on
the rot-283 restart arms show the same picture (ratios 0.56–1.02).

Note `cs_recon_ab/README.md` §1: on the cubed sphere the radial sweep hardcodes
position-aware PLM, so this experiment varies the **angular** reconstruction only. The
null result is therefore consistent with the radial scheme, not the angular one, being
what matters for this problem.

---

## Ohmic cap in `sphbeam`

`ohmic.py` + `xe_model.py`, tables in `ohmic_tables.md`, figure `fig_ohmic.png`.
Both per-cell limits are rebuilt from each arm's **last dump**, exactly as the code
computes them:

* CFL — `mhd_newdt.cpp`: `dt_i = cfl dx_i/(|v_i| + cf_i)`, `cf_{2,3}` divided by
  `sin_cell` (this is `cfl.py`, reused unchanged);
* Ohmic — `Resistivity::NewTimeStepGeneralResist`: `dt_i = cfl dx_i^2/(6 eta)` with
  `dx` the *physical* curvilinear length, and
  `eta = ResistivityEOS(x_e, T, max_eta) = 230 sqrt(T)/max(x_e, 230 sqrt(T)/max_eta)
  + 5.2e11 * 20 / T^1.5`, clipped at `max_eta`.

All four arms set `<mhd>/ohmic_resistivity = eos`, `max_eta = 1e13`, `use_rkg_sts =
false`, `cfl_number = 0.3`. `x_e` therefore comes from the EOS's own Saha table, which
the dumped `eos_table.txt` does **not** carry (it has only `log10 e/rho` and
`log10 p/rho`). `xe_model.py` is a vectorised python port of
`src/eos/eos_composition.hpp` — H2 dissociation, Saha for H/He/He+, the six metal
donors (Na, K, Ca, Al, Mg, Fe) and the two-pass per-species condensation rainout — built
on the *same* node grid (`log10 rho` −14…0 by 0.05, `log10 T` 1.5…6 by 0.01) the run
tabulates, then bilinearly interpolated (the code uses a bicubic Hermite patch on those
same nodes).

**Validation.** The two limits combined reproduce each arm's running dt to 0.1–2.5 %:

| arm | run dt (median, last 200 cycles) | reconstructed | CFL-only | Ohmic-only | binds | cells at `eta = max_eta` |
|---|---|---|---|---|---|---|
| prodbin | 19.92 | 19.81 | 19.81 | 22.45 | **CFL** | 0.00 % |
| sph     | 12.85 | 12.53 | 17.36 | 12.53 | **Ohmic** | 1.43 % |
| sphbeam | 12.85 | 12.85 | 18.00 | 12.85 | **Ohmic** | 0.30 % |
| offfix  | 19.96 | 19.93 | 19.93 | 322.30 | **CFL** | 0.00 % |

So the cold night side *is* still the whole story of the cap — it is simply not the
night *top* (which `ck_beam_sph` does fix) but the night **r/Rp 1.23–1.28 shell**, which
the beam barely touches.

### (1) Which cells hold `sphbeam` at 12.85 s

The 20 most limiting cells are all **Ohmic**, all on the **night side** (100 % of capped
cells in both `sph` and `sphbeam` are at |lon_ss| > 105°; not one day-side or
terminator cell is capped in any arm), all in one shell:

| arm | binding shell | r/Rp | state of the limiting cells |
|---|---|---|---|
| sph | i = 41 | 1.236 | T 1385–1464 K, rho 1.4–1.7e-7, p 7–9 mbar, x_e 3–9e-10, beta 1e2–2e4 |
| sphbeam | i = 42 | 1.241 | T 1406–1461 K, rho 1.1–1.5e-7, p 5.6–7.8 mbar, x_e 5–9e-10, beta 5e1–4e3 |

They are **the same night-side shell cells as in `sph`, capped for the same reason**, not
a new region: the beam does not create a twilight population at the cap. The only thing
that changed is how many there are and how deep they reach:

| i | r/Rp | dx_min [cm] | prodbin ncap / night T 1 % | sph | sphbeam | offfix |
|---|---|---|---|---|---|---|
| 40 | 1.230 | 4.95e7 | 0 / 2084 | 0 / 1716 | 0 / 1756 | 0 / 2074 |
| 41 | 1.236 | 5.01e7 | 0 / 1992 | **5** / 1566 | 0 / 1605 | 0 / 2001 |
| 42 | 1.241 | 5.07e7 | 0 / 1903 | 31 / 1446 | **18** / 1470 | 0 / 1948 |
| 43 | 1.246 | 5.14e7 | 0 / 1862 | 127 / 1371 | 67 / 1406 | 0 / 1898 |
| 45 | 1.257 | 5.28e7 | 0 / 1771 | 123 / 1329 | 33 / 1387 | 0 / 1847 |
| 48 | 1.275 | 5.52e7 | 0 / 1745 | 175 / 1240 | 25 / 1336 | 0 / 1808 |

Because every capped cell has the *same* eta (the ceiling), the binding cell is simply
the **innermost capped cell**, where `dx_min` is smallest. `sphbeam` recovered exactly
one shell — i = 41 is now clear — which is worth 12.53 → 12.85 s, 2.6 %. It did not
recover i = 42. The beam's warming of this shell is only 20–40 K at the cold 1st
percentile (1446 → 1470 K at i = 42), whereas production sits 430 K higher there
(1903 K). The "within 45 K of production" recovery is real at the night **top**; at
r/Rp 1.24 the night side is still ~350 K (median) to ~430 K (cold tail) below `prodbin`.

`fig_ohmic.png` maps `log10 eta` on each arm's binding shell (lat vs lon from
substellar, terminators marked): `prodbin`/`offfix` show a smooth day–night eta contrast
topping out near 1e12, while `sph`/`sphbeam` show a mottled night hemisphere saturated
at the 1e13 ceiling.

### (2) The temperature at which eta reaches `max_eta`

`eta = max_eta` exactly when `x_e <= 230 sqrt(T)/max_eta` (~8e-10 near 1400 K), i.e.
when the alkali donors are not yet ionised. Bisecting `eta(rho,T) = max_eta` on T:

| rho [g/cm^3] | 1e-10 | 1e-9 | 1e-8 | 1e-7 | 1e-6 | 1e-5 | 1e-4 |
|---|---|---|---|---|---|---|---|
| T_cap [K] | 1209 | 1278 | 1356 | **1444** | 1544 | 1658 | 1790 |

At the binding shell's density (~1.3e-7) the threshold is **T_cap ≈ 1450–1475 K**.

| arm | median T of capped cells | median T_cap there | median deficit | deficit at the binding shell |
|---|---|---|---|---|
| sph | 998 K | 1096 K | 89 K | **6 K** |
| sphbeam | 1000 K | 1105 K | 87 K | **14 K** |
| prodbin / offfix | — (no capped cells) | — | — | — |

The binding cells are **marginal** — 6 K (sph) and 14 K (sphbeam) below their own
threshold. The capped population as a whole (mostly higher up, r/Rp out to 2.17, where
`dx` is too large to bind) sits ~90 K below. `prodbin`'s coldest night cell anywhere in
i = 38…49 is 1732 K, ~300 K above threshold; `offfix`'s is 1755 K.

### (3) What would change the binding limit

**Raising `max_eta` makes it strictly worse** — the Ohmic dt goes as `1/eta`, and the
capped cells are pinned *at* the ceiling, so the ceiling is the timestep:

| arm | CFL-only | 1e12 | 3e12 | 5e12 | **1e13 (current)** | 3e13 | 1e14 |
|---|---|---|---|---|---|---|---|
| prodbin | 19.81 | 19.81 | 19.81 | 19.81 | 19.81 | 19.81 | 19.81 |
| sph | 17.36 | 17.36 | 17.36 | 17.36 | **12.53** | 4.28 | 1.29 |
| sphbeam | 18.00 | 18.00 | 18.00 | 18.00 | **12.85** | 4.40 | 1.57 |
| offfix | 19.93 | 19.93 | 19.93 | 19.93 | 19.93 | 19.93 | 19.93 |

**Lowering** `max_eta` to 5e12 hands both `sph` and `sphbeam` back to the CFL (17.4 and
18.0 s) at no cost to `prodbin` or `offfix`, which never touch the ceiling — a pure
input-side change, unrelated to RT. That is the cheapest lever, at the price of a
factor-2 weaker resistivity cap in the coldest night cells (which are unresolved dead
zones there in any case).

**Warming the twilight** works too, but needs a lot of it. Applying a uniform `dT` and
recomputing `x_e`, eta and the Ohmic limit (CFL held fixed):

| arm | dT = 0 | 25 | 50 | 100 | 200 | 400 K |
|---|---|---|---|---|---|---|
| sph | 12.53 | 12.53 | 12.53 | 12.85 | 14.77 | 17.36 |
| sphbeam | 12.85 | 12.85 | 13.19 | 13.19 | 16.18 | **18.00** |

`sphbeam` needs ~+400 K on the night r/Rp 1.24 shell to reach its CFL ceiling — which is
exactly the ~430 K deficit against `prodbin` measured above, and consistent: the thing
`ck_beam_sph` has not yet fixed is the deep night shell, not the night top.

**Ceiling if the cap were removed entirely** (CFL-only): `sph` 17.36 s, `sphbeam`
18.00 s, against `prodbin` 19.81 s. So even a perfect fix to the Ohmic cap leaves
`sphbeam` ~9 % short of production dt; the residual is a genuine MHD fast-mode
difference on the night side, not resistivity.

---

## What could not be determined from this data

1. **Whether `ck_spherical` alone causes the dt loss, the magnetic decay or the energy
   drift — now resolved in direction with `offfix` (§"Isolation of ck_spherical"), not
   with a bitwise-identical control.** `offfix` (HEAD binary `8e41808b`, repaired
   plane-parallel path, same restart) reproduces `prodbin`'s dt plateau and night-side T,
   and shows worse ME decay and a much larger/sign-changing `rt_desum` gap than `prodbin`
   even with the switch off — so a real binary/lineage gap to `prodbin` remains (item 3
   below is now measurable, not closed). But `offfix` vs `sph` isolates the switch itself:
   turning `ck_spherical` on costs +533 K of night-side cooling at 1e-6 bar, collapses dt
   from 19.9 to 12.5 s via the Ohmic cap, roughly doubles the ME decay rate, and only
   modestly (~1.22x) raises per-cycle compute cost. `offfix` and `sph` are still not the
   same literal binary, so this is not a bitwise A/B; a same-binary `ck_spherical=false`
   rerun (i.e. `sph`'s own binary with the flag flipped) would remove the last gap.
2. **Where `eta` saturates `max_eta`.** The dt arithmetic identifies the shell
   (r/Rp 1.230–1.241) to six digits, but resistivity and electron fraction are not in
   `mhd_w_bcc`, so the map of the capped region and its change between arms is inference
   from T(r), not measurement. `offfix` (which never saturates the cap) is consistent with
   this being caused by `ck_spherical`'s cooling of that shell, not an independent effect.
3. **The plane-parallel `rt_desum` baseline** — `prodbin`'s binary has no such counter, so
   the true production-binary value is still unmeasured. `offfix` gives a HEAD-binary
   plane-parallel value instead (mean −27.9 %, range −94.6 %…+47.5 %, 28 lines) — larger
   in magnitude than the −7.8 % "on" arms, matching the 1-D column's direction, but with a
   sign-changing spread the 1-D result did not show, so the exact 36–40 %/2–5 % figures
   from `tests_ck_sph/README.md` §0.3 are still not reproduced in 3-D.
4. **Same-instant comparisons.** Every arm stopped at its own wall clock, so all dump
   comparisons are at times differing by 0.02–0.06 rot (comparison 1) or 0.035 rot
   (comparison 2). `offfix`'s dump is at t = 8.6667e7 (rot ≈ 284.13), 0.02–0.06 rot from
   the other comparison-1 dumps, same caveat. All history comparisons are exact in time;
   all dump comparisons are not.
5. **Whether the `noang` difference is physics or chaotic divergence.** Both arms are the
   same binary with one parameter changed, but no bitwise-repeat control was run, so the
   2e-3 level differences cannot be separated from amplified round-off.
6. **Whether the 364-commit binary gap seen in `offfix` vs `prodbin` (6x faster ME decay,
   9.5x higher `eos_efloor`/cycle, 1.46x per-cycle cost, the first-row ME step) is one
   regression or several**, and whether any of it is itself `ck_spherical`-adjacent code
   that changed behavior even with the flag off. Not investigated here — would need
   bisection between `27ca5b13` and `8e41808b`.
