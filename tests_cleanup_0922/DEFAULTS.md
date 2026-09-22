# Default-value changes, 2026-09-22

Branch `rt-integration`, on top of `ac42334d`.  Six switches, one section each: where the
default lives, what it changes bitwise and for which production, and the gate result.
Nothing is committed.

## Harness

`tests_cleanup_0922/{build.sh,gates.sh,compare.sh,rstgate.sh}` reused unchanged.
`gates_old.sh` is new: it is `gates.sh` with every flipped switch written back into the
input file EXPLICITLY at its old value (only into inputs that were SILENT about it), and
is what gate (a) runs.

The run trees were deleted at the end; the `*.md5` manifests remain, and `ref/` (bitwise
equal to the regenerated `ref2`) is still the reference tree.

Reference binaries were REBUILT from the unmodified HEAD first, because the `ref/` tree
predates the clean-up.  `compare.sh ref ref2` -> `same=52 differ=0`, i.e. the clean-up
commits are still bitwise neutral and `ref2` is a valid reference.

| gate | what | result |
| --- | --- | --- |
| (a) `compare.sh ref2 oldexp2` | new binary, old values stated explicitly | **same=52 differ=0** |
| (b) `compare.sh ref2 newdef2` | new binary, inputs untouched | **same=46 differ=6** — the six files are the two dhj runs only |
| restart | `rstgate.sh newdef` | straight/chain rows match bit-for-bit; `rt_surface.bin`, `rt_profile.bin` IDENTICAL, and also identical to `rst_ref2` |

The six differing files under (b) are `c_dhj_{false,true}/{dhj.hydro.hst, photo.txt,
bin/dhj.hydro_w.00001.bin}`.  Box (a1, a2 m3/m0), red giant (b1, b2), cubed-sphere blast
(d1, d2) and spherical-polar MHD blast (e) are byte-identical, which is what switches
4 and 5 predict (see below).  Magnitude of the dhj change at cycle 20: mass and total
energy unchanged to 6 digits, 1-KE 8.86194e31 -> 8.86230e31 (4e-6), 1-mom 1e-6; the only
large *relative* change is 2-mom, which is a 1e14 residual of a cancelling sum against a
1e28 1-mom.  An O(dt) change in the radiative source, as expected.

Style: `git diff` adds no line over 90 columns, no tab, no `}}`, no trailing whitespace;
`src/` files stay mode 644.  cpplint run on the four touched C++ files before and after
produces the identical set of categories (only the pre-existing `whitespace/line_length`
and `readability/fn_size` counts move, both filtered by `CPPLINT.cfg`).  `flake8` clean on
the touched test.

Tests: `tst/test_suite/rad/test_rad_dhj_ck_cpu.py` and `..._mpicpu.py` both pass (the
former after the one edit recorded under switch 1/2).  `tst/test_suite/nr/
test_nr_geneos_repro_gpu.py` is GPU-only and is a same-binary self-comparison, so it is
unaffected; its grid has nx1 = 64 <= RT_NNC.

---

## 1. `problem/ck_spherical`

**Default line.** `src/pgen/deep_hot_jupiter_rt.cpp:461` (was line 438,
`GetOrAddBoolean("problem","ck_spherical",false)`).  The flag itself lives at
`src/utils/two_stream_rt.hpp:865`.  It is read by **no other problem generator**.

**New default.** `(pmesh->use_spherical_polar || pmesh->use_cubed_sphere) && rt_ck &&
!rt_layer_legacy`.  The mesh flags tested are exactly `pmy_mesh_->use_spherical_polar`
and `pmy_mesh_->use_cubed_sphere`, captured at the top of `UserProblem` as
`use_spherical_polar` / `use_cubed_sphere_`.

The two extra conditions are not decoration: the flag's own guards in
`two_stream_rt.hpp:1524-1551` FATAL if it is true with `rt_ck` false ("converts the
CORRELATED-K kernel") or with `rt_layer_legacy` true ("the staggered whole-cell ck layers
were NOT converted").  Defaulting it true unconditionally on a curvilinear mesh would
therefore turn every grey / picket-fence spherical dhj run and every legacy-reproduction
run into a startup fatal.  `rt_layer_legacy` is read further down, so the default peeks at
it with a `GetOrAddBoolean` that the later read then returns unchanged.

**Cartesian refusal.** Already present at `two_stream_rt.hpp:1525-1532`, but only at the
first RT call.  A startup refusal was added at `deep_hot_jupiter_rt.cpp:465-476`, covering
both this switch and `ck_beam_sph`.  VERIFIED by running the dhj binary on a Cartesian
copy of `inputs/tests/dhj_ck_spherical.athinput` (`use_spherical_polar=false`,
`use_polar_boundary=false`, both stretches off, real x2/x3 extents, periodic x2):

* with `problem/ck_spherical=true` it dies with
  `### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_spherical and problem/ck_beam_sph
  are the SPHERICAL forms of the correlated-k thermal two-stream and of the direct
  stellar beam.  They have no meaning on a Cartesian mesh ...`;
* with the flags absent (the Cartesian default, false) it passes the RT setup, runs the
  two-stream self-tests, and dies later on the pre-existing, unrelated
  `problem/grav_point_mass requires a radial x1`.

So the refusal path is exercised for real, not argued from the code.

**Bitwise effect.**
* dhj gate (c): none.  `inputs/tests/dhj_ck_spherical.athinput:320` states
  `ck_spherical = false` explicitly, and `gates.sh` overrides it on the command line
  anyway.
* `inputs/production/deep_hot_jupiter_cs_prod4.athinput:264` already states `true`: no
  change.
* PRODUCTIONS THAT CHANGE: `bench/sp_mhd_prod3/deep_hot_jupiter.athinput` and
  `bench/cs_mhd_prod3/deep_hot_jupiter.athinput` say nothing about it and both run
  `rt_ck = true` on a curvilinear mesh, so both pick up the spherical thermal two-stream.
* `tst/test_suite/rad/*` run `inputs/mhd/deep_hot_jupiter_rt_eos.athinput`
  (`use_spherical_polar = true`, `rt_ck = true`), so they now exercise the spherical form.

**Inputs touched.** None.  No Cartesian input anywhere in `inputs/` or `tests_*/` names
either flag, so nothing contradicts the refusal.

## 2. `problem/ck_beam_sph`

**Default line.** `src/pgen/deep_hot_jupiter_rt.cpp:462` (was 442).  Same rule, same
`cksph_default` expression, same startup refusal, same guards
(`two_stream_rt.hpp:1492-1519`).  Read by no other problem generator.

**Bitwise effect.** Identical set to switch 1: nothing in the gates (the test input states
it false), `sp_mhd_prod3` and `cs_mhd_prod3` change, `deep_hot_jupiter_cs_prod4` already
stated true.

**Gate number that had to be regenerated.**
`tst/test_suite/rad/test_rad_dhj_ck_cpu.py` asserted, on a NIGHTSIDE column,

```python
assert np.all(qsw == 0.0), "shortwave heating on the nightside"
```

and this FAILS under the new default, correctly: the pseudo-spherical beam lights
twilight columns (mu0 < 0 down to the grazing radius), which is the documented purpose of
the flag.  Old reference: exactly 0.0 in every cell.  New: one cell at
`qsw = 2.73615533e-28` erg/cm^3/s, everything else 0.  The assertion was replaced by a
sign check plus a depth-integrated bound against the insolation,

```python
night_abs = np.sum(qsw[:-1] * np.diff(col[:, 1]))
assert night_abs < 1.0e-6 * ck.SIGMA_SB * hdr["T_irr"] ** 4
```

(measured margin ~1e-28 of the insolation).  The test passes; so does the MPI sibling,
which compares ranks against each other and hard-codes nothing.

## 3. `problem/rt_use_cons`

**Default line.** `src/pgen/deep_hot_jupiter_rt.cpp:545` (was 509), `false` -> `true`.
The file header comment (lines 11-31) was rewritten to match.  `box_convection.cpp:1835`
already defaults true; `red_giant.cpp:1691/1720` keep the mode-3 flag `col3` as the
default and are NOT touched (their header comment was corrected to say so).

**Bitwise effect.**
* `bench/sp_mhd_prod3/deep_hot_jupiter.athinput` NEVER SETS the line — it is the run that
  changes.
* `bench/cs_mhd_prod3/deep_hot_jupiter.athinput:326` already has `rt_use_cons = true` —
  unchanged.
* This is the ONLY switch responsible for the six differing files in gate (b):
  `inputs/tests/dhj_ck_spherical.athinput` is silent about it.

## 4. `problem/rt_col3_ex_iter`

**Default lines.** `src/pgen/box_convection.cpp:1859` and `src/pgen/red_giant.cpp:1739`,
both `false` -> `true`.  `deep_hot_jupiter_rt.cpp` never reads it (it refuses
`rt_implicit_column != 0` outright).

**Where it is consulted.** Exactly one site: `src/utils/two_stream_rt.hpp:3020`,
`c3.ex_iter = rt_col3_ex_iter;`, inside the `rt_implicit_column == 3` branch.  Every other
mention is a validation message.  So it is inert in modes other than 3.

**THE CHECK THE USER ASKED FOR** — run on the UNMODIFIED HEAD binaries, flag stated
explicitly false vs explicitly true (the flag had to be written into a copy of each input
first: AthenaK refuses a command-line override of a parameter the file does not contain):

| run | files compared | result |
| --- | --- | --- |
| `red_giant_column.athinput`, nlim 300 | 4 | **identical** |
| `red_giant_fofc.athinput`, nlim 20, 8x8 | 5 | **identical** |
| `dhj_ck_spherical.athinput`, nlim 20, `ck_spherical=false` | 4 | **identical** |
| `dhj_ck_spherical.athinput`, nlim 20, `ck_spherical=true` | 4 | **identical** |

So the flip is safe and was made.

**Productions.** Both FeCZ productions already state it:
`bench/hestar_fecz/box_w8/he_box_w8.athinput:282` and
`bench/bstar_fecz/prod_w7/fecz_rt_w7.athinput:349`, both `= true`.  Nothing changes for
them, and nothing in `tests_gate_merge/postmerge.sh` had to be regenerated (that script
hard-codes no numbers at all — it `cmp`s against a previously-run `g1_new_m$M` tree).

## 5. `problem/mlt_split_deposit`

**Default line.** `src/pgen/red_giant.cpp:1441` (was 1431), `false` -> `true`.

**box_convection does NOT read this parameter.**  The brief said "box_convection and
red_giant"; grep over `src/` finds it only in `red_giant.cpp`.  Nothing was added there.

**When the flag matters.** Two sites, both requiring `rt_strang_`:
* `red_giant.cpp:4485`, the split deposit: `if (mlt_split_dep_ && rt_strang_ &&
  mlt_alpha_ > 0.0)`, inside `RedGiantRTSplitBody`, which itself returns early unless
  `rt_ck_ || rt_grey_`;
* `red_giant.cpp:3755`, the stage-local deposit, skipped by
  `if (!(mlt_split_dep_ && rt_strang_))`.

So it has an effect only when ALL of `problem/rt_strang`, `problem/mlt_alpha > 0` and the
two-stream (`rt_ck` or `rt_grey`) are on.  That is exactly why the clean-up kept
`rt_strang`: with `rt_strang = false` there is no split half step to deposit into and the
stage-local deposit runs whatever this flag says.

**Productions that change bitwise: NONE of the three named.**
* B-star `bench/bstar_fecz/prod_w7/fecz_rt_w7.athinput` and He
  `bench/hestar_fecz/box_w8/he_box_w8.athinput` are `box_convection` runs, which never
  read the parameter.
* RG `bench/RG_fofc_long3/rg.athinput` is `red_giant` but sets `mlt_alpha = 0.0`
  (line 189) and does not set `rt_strang`, so both guards are false.

Gate (b) confirms it: `b1_rg_col` and `b2_rg_fofc` are byte-identical
(`red_giant_column.athinput` has `mlt_alpha = 0.0` and no two-stream;
`red_giant_fofc.athinput` has `mlt_alpha = 0.0`).

## 6. `problem/rt_split`

**Default line.** `src/pgen/deep_hot_jupiter_rt.cpp:442` (was 433), `false` -> `true`.
The namespace initialiser `two_stream_rt.hpp:233` is left at `false`; every pgen sets it.

**NOT DELETED, and the parameter stays.**  The serial path IS reachable: `rt_ck` (dhj) and
`rt_grey` (box_convection:1818, red_giant:1651/1942) force it true, but
`deep_hot_jupiter_rt` with `rt_ck = false` runs the GREY PICKET FENCE, which sets neither,
so `rt_split = false` still selects the monolithic kernel there.
`inputs/mhd/deep_hot_jupiter_rt_ideal_xe.athinput` is exactly such an input.  Per the
brief's own instruction ("if some pgen path can still run with it false, keep the
parameter default true and report"), the flag and the serial kernel were kept.

**TWO FINDINGS THE USER SHOULD SEE BEFORE ACCEPTING THIS FLIP.**

1. *The split grey sweep sizes its column at compile time.*
   `two_stream_rt.hpp:2097` fatals with `problem/rt_split with grey RT needs RT_NNC >= n1`
   when `n1 > RT_NNC` (default 72), which the monolithic path does not.  Every
   picket-fence input in the repo has `nx1 = 64`, so nothing in-tree breaks, but a
   picket-fence run at nx1 > 72 that relied on the old default will now refuse to start.

2. *The two paths are NOT bitwise equal on the picket fence*, contrary to the comment that
   was there.  MEASURED on `inputs/mhd/deep_hot_jupiter_rt_ideal_xe.athinput`, 20 cycles,
   same binary, `problem/rt_split=false` vs `=true`:
   `dhj.mhd.hst` and the field dump both differ; mass and total energy agree to the
   printed digits, but **1-KE differs by 7.1 %** (2.20858e30 vs 2.05090e30).  A rerun of
   the `false` case is byte-identical to itself, so the run is deterministic and the
   difference is real.  Unexplained; it is confined to the picket-fence path (the grey
   and correlated-k kernels are covered by gates (a)/(b) and are identical).
   Both code comments were corrected to state this.

**Bitwise effect on the gates and productions.** None.  Every dhj input in `inputs/`,
`tests_*/` and `bench/{sp,cs}_mhd_prod3` that states `rt_split = false` also states
`rt_ck = true`, which forces it true regardless, and the box / red-giant pgens force it
true.  Gate (b) shows no change attributable to it.
