# `tests_gpu_ab_0922` — analysis of the 2026-09-22 apudev smoke A/B arms

Read-only analysis of `bench/ck_sph_ab/{prodbin,sph,sphbeam}`, `bench/cs_noang/noang`
(comparison 1) and `bench/cs_recon_ab/{plm_ctl,ppmx,wenoz}` (comparison 2).
Nothing in `bench/` was written to.

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

Clean one-line comparisons in this set: `sph` vs `sphbeam`, and `sph` vs `noang`.

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

**In `prodbin` the MHD CFL sets dt** (19.81 reconstructed vs 19.82 reported — a 0.05 %
match, which also validates the metric and the EOS inversion). Its limiting cells are
real gas: r/Rp 1.25, equator, night-ish (lon +150…+165), p 1e-2…1e-6 bar, T 1900–2700 K,
|v| 6–12 km/s, none on `dfloor`. The limit is radial in every one of the top 10.

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

Density flooring is comparable (`sphbeam` is 25 % *lower* than `prodbin`), but energy and
temperature flooring are **130x and 49x higher** in `sph`/`noang` than in `prodbin`, and
`sphbeam` takes about 5/8 of that back. This is the same cold night-side top seen in Q2 —
the spherical thermal form cools the thin gas onto `tfloor_kelvin = 200` far more often.
The claim that the spherical form floors *less* in the upper atmosphere is **not
supported here**; the opposite is seen (with the binary caveat above).

Location, last dump: floored cells (`rho <= dfloor` to 1e-4) are 5.00 % of the mesh in
`prodbin`, 5.79 % in `sph`, 5.82 % in `noang`, 3.23 % in `sphbeam`. They are almost purely
night-side (day-side fraction 0.003–0.014), spread over all longitudes with mean |lat| 28–31°,
and the **inner edge of the floored slab moves down** from r/Rp 1.508 (prodbin) to 1.404
(noang) / 1.381 (sph), while `sphbeam` keeps it at 1.499. `p <= pfloor` occurs in 242
(sph) / 230 (noang) / 1 (sphbeam) / 0 (prodbin) cells.

## Q4 — energetics and `rt_desum`

`fig_energy.png`. Over the run: `prodbin` |dE/E| and |dM/M| are below the 6-digit history
precision (< 2e-6) over a full rotation; `sph`/`noang` lose 2.56e-4 of total energy and
2.8e-5 of mass over 0.47 rot, `sphbeam` 2.03e-4 over 0.35 rot. Magnetic energy:
4.98→4.71e33 (−5 % / rotation) in `prodbin`, 4.06→2.92e33 (−28 % in 0.47 rot) in
`sph`/`noang`, 4.06→2.95e33 in `sphbeam`. Kinetic energy rises ~8 % in the HEAD arms and
falls ~6 % in `prodbin`.

**The expectation that the spherical form is the more conservative one is not confirmed —
but neither is it refuted**, because the invalid `off` arm shows the same initial magnetic
step and decay on the HEAD binary with the switch *off*.

`rt_desum`: `prodbin`'s binary has no `rt_cell_report` parameter (`strings` count 0), so
the plane-parallel side of the 36–40 % vs 2–5 % comparison **cannot be measured here at
all**. On the 3-D grid the spherical arms give mean −7.8 % (range −12.4 %…−4.4 %), sph and
sphbeam agreeing to 0.1 pp — i.e. consistent with the "on" end of the 1-D column result,
with no baseline.

## Q5 — cost

| arm | cycles | wall s | cycles/s | sim s per wall s | relative |
| --- | --- | --- | --- | --- | --- |
| prodbin | 15200 | 601 | 25.3 | 505 | 1.00 |
| sph | 10900 | 768 | 14.2 | 182 | 0.36 |
| sphbeam | 8700 | 764 | 11.4 | 147 | 0.29 |
| noang | 10800 | 763 | 14.2 | 181 | 0.36 |

Per cycle the HEAD arms are 1.78x more expensive (`sphbeam` 2.2x); the smaller dt costs a
further 1.58x, so the throughput penalty is **2.8x** (`sphbeam` 3.4x). Part of the 1.78x
is the binary difference, not the switch. `noang` is **not** cheaper than `sph`
(14.16 vs 14.20 cycles/s, 0.3 %), so the transverse conduction planes are not a measurable
cost here — contrary to the expectation in `cs_noang/README.md`.

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

## What could not be determined from this data

1. **Whether `ck_spherical` alone causes the dt loss, the magnetic decay or the energy
   drift.** There is no valid HEAD-binary plane-parallel control; `ck_sph_ab/off` is the
   regressed hybrid. The `off` arm does show that the first-output magnetic step is a
   binary/restart effect, not the switch. A 15-min `ck_spherical = false` run on the
   `916dc953` binary would settle all three at once.
2. **Where `eta` saturates `max_eta`.** The dt arithmetic identifies the shell
   (r/Rp 1.230–1.241) to six digits, but resistivity and electron fraction are not in
   `mhd_w_bcc`, so the map of the capped region and its change between arms is inference
   from T(r), not measurement.
3. **The plane-parallel `rt_desum` baseline** — `prodbin`'s binary has no such counter.
4. **Same-instant comparisons.** Every arm stopped at its own wall clock, so all dump
   comparisons are at times differing by 0.02–0.06 rot (comparison 1) or 0.035 rot
   (comparison 2). All history comparisons are exact in time; all dump comparisons are not.
5. **Whether the `noang` difference is physics or chaotic divergence.** Both arms are the
   same binary with one parameter changed, but no bitwise-repeat control was run, so the
   2e-3 level differences cannot be separated from amplified round-off.
