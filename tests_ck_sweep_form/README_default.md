# `problem/ck_sweep_form` default 0 -> 1 (tm) under `ck_spherical` — the gate log

2026-09-22, user decision.  Companion to `README.md` (what the two probe-free forms are
and what the probe's lag was worth) and to `tests_ck_sph/README.md` (the spherical form
itself).  **Code only, nothing committed here.**

## What changed

`src/pgen/deep_hot_jupiter_rt.cpp`, the `problem/ck_sweep_form` read: the default is now

    1 (tm)   if  ck_spherical && !ck_implicit && ck_sweep_cache == 2
    0        otherwise

i.e. it is conditioned on **every one of the flag's own refusals** in
`two_stream_rt.hpp` (`!ck_spherical`, `ck_sweep_cache != 2`, `ck_implicit`), so that the
flip cannot turn an input that starts today into a startup fatal.  The `ck_sweep_cache`
condition is not one the decision named; it is included because the refusal exists and
an input carrying `ck_sweep_cache = 1` with `ck_spherical = true` would otherwise die.
`ck_implicit` is read in full further down, so it is peeked at here exactly as
`rt_layer_legacy` already was.  The refusals themselves are untouched.

The resolved value is stated at startup on rank 0, since it is no longer readable off the
input file:

    ### deep_hot_jupiter_rt: problem/ck_sweep_form = 1 (tm), the default here (ck_spherical)
    ### deep_hot_jupiter_rt: problem/ck_sweep_form = 0 (four-pass (probe)), set explicitly

Why: the four-pass form's probe-lagged mixing constant is only `O(beta^2)` — 2.3e-2 of
the face flux at the production area ratio 3.45 — while `tm` and `sd` solve the coupled
column exactly (they agree to 3.5e-17), and `tm` is 1.36x faster on the GPU
(`README.md`).  `src/utils/two_stream_rt.hpp`'s flag note records the same.

`inputs/production/deep_hot_jupiter_cs_prod4.athinput` keeps `ck_sweep_form = 1` set
explicitly and its header note 9 is still correct; it was not touched.

## Gates

Binaries: `build_ckdef` (`-D PROBLEM=deep_hot_jupiter_rt`, Release, CPU) at HEAD
`8b8f246b` = REFERENCE, and the same tree with the change.  All runs 20 cycles,
`inputs/tests/dhj_ck_spherical.athinput` (128 x 8 x 8, the production radial grid) with
`problem/ck_spherical=true` on the command line, `output1/dt=1` (history every cycle) and
`output3/dt=1` (22 `hydro_w` dumps).  "bitwise" = identical `dhj.hydro.hst` (md5) AND
every variable array of all 22 dumps equal under `np.array_equal`.  `ck_sweep_form` had
to be added to a scratch copy of the input where it is set on the command line, because
`parameter_input` refuses a command-line parameter that no block contains.

| # | configuration | reference | new binary | result |
| --- | --- | --- | --- | --- |
| a | `ck_sweep_form = 0` explicit, spherical | form 0 | form 0 | **bitwise** |
| a | `ck_implicit = true` + `rt_use_cons = true`, `ck_sweep_form = 0` explicit | form 0 | form 0 | **bitwise** |
| b | default, spherical | form 1 explicit | form 1 (default) | **bitwise** |
| b | default, `ck_spherical = false` | form 0 | form 0 (fallback) | **bitwise** |
| b | default, `ck_implicit = true` | form 0 | form 0 (fallback), starts | **bitwise** |

(`dhj_ck_implicit.athinput` also needs `problem/rt_use_cons=true`; with the file's own
`false` the implicit solve refuses at the first RT call, which is pre-existing.)

### (c) the regression test

`cd tst && python3 -m pytest test_suite/rad/test_rad_dhj_ck_cpu.py` — **1 passed in
154 s**, with no threshold touched.  Its input
(`inputs/mhd/deep_hot_jupiter_rt_eos.athinput`) is `use_spherical_polar` and the test adds
`rt_ck=true`, so under the new default it really does run `tm`; reproduced by hand with
the same overrides:

    ### deep_hot_jupiter_rt: problem/ck_sweep_form = 1 (tm), the default here (ck_spherical)
    ### two_stream_rt: correlated-k PROBE-FREE sweep ON (problem/ck_sweep_form = 1, tm ...)

So no threshold in that file was tuned tightly enough to see the probe lag, and none was
regenerated.  The MPI sibling `test_rad_dhj_ck_mpicpu.py` was **not run**: no MPI in the
shell's environment on this node (`mpirun` absent).

### (d) the `tests_ck_sph` thermal gates, with the new default

`tests_ck_sph/budget2.py` on a `time/nlim=1` column dump, production radial grid, default
form (= tm), against the `ck_spherical = true` row of `tests_ck_sph/README.md` 0.2/0.3:

| column | max per-cell `|V dep - dPhi|/max|Phi|` | column budget | README (4-pass) |
| --- | --- | --- | --- | --- |
| `m0k5`, mu0 = 0.922 | 2.822e-10 | 3.864e-09 | 2.05e-10 / 2.62e-09 |
| `m1k3`, mu0 = 0.382 | 2.822e-10 | 3.864e-09 | 2.05e-10 / 2.62e-09 |
| `m0k5`, form 0 here | 2.808e-10 | 3.871e-09 | — |

Round-off at `RtF` precision, as before: the identity and the budget are properties of
the telescoping deposit, which both forms share.  The beam block is bit-identical to the
four-pass run (absorbed power 1.31810e+29 / 5.66052e+28, `tau_slant = 1` at 1.24093e+10 /
1.25970e+10, top-face flux ratio 0.99133 / 0.99107) — as it must be, the beam not being
swept by this flag.

Star-off `L(r)`: the last restart of `g_relax_f1` advanced one cycle with the new binary
and the night column (m = 0, k = 2) dumped, then `tests_ck_sph/lprof.py` — **every line
identical to the stored `h_night_f1/col.txt`**, `L(r)/L_int` flat at 992.08 (in units of
the script's normalisation) from the thin region to the top face.

Style: `cpplint` (repo `CPPLINT.cfg`) over the two touched files gives the **same 270
findings as HEAD**, i.e. the change adds none; all added lines are <= 90 columns, no tabs,
no trailing whitespace.

The `build_ckdef` build directory and every run directory made for these gates were
deleted afterwards; the numbers above are the record.
