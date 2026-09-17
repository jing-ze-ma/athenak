# he4_presn on the cubed sphere: the merge, the 1-D gate, and what stopped it
(2026-09-17, viper, branch `he4-presn-global`)

This directory holds the STEP-1 merge check and the STEP-2 1-D gate of the global He-star
run.  **Step 3 (the 3-D smoke run) was NOT reached**: the 1-D gate fails, and per the
brief the work stops at the diagnosis.  All `.bin`/`.rst`/`.cbin` have been deleted and
the two large logs truncated (every number below was read off the full log first); the
`.hst` files, the truncated logs and the submit scripts are kept, and each table's
command is the `sub_*.sh` next to it.

Build: one GPU directory, `../build_gpu_rg` (`-DPROBLEM=red_giant`, HIP,
`Kokkos_ARCH_AMD_GFX942_APU`, MPI on).

    module purge && module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
    cd ../build_gpu_rg && make -j24

---

## STEP 1 -- the merge  -- CLEAN

`git merge --no-edit cs-implicit-transverse` (`a8222bd3`) into `244b1e21`:
**no conflict**, 32 files, +2826/-48.  The two branches are disjoint: the merged branch
owns `src/diffusion/*`, `src/pgen/cs_test.cpp`, `tst/test_suite/rad/`, `tests_adi/` and
`docs/dev/cs_implicit_transverse.md`, while `he4-presn-global` owns `src/pgen/red_giant.cpp`,
`src/utils/two_stream*` and `inputs/hydro/he4_presn_cs.athinput`.

`inputs/hydro/he4_presn_cs.athinput` then switched the transverse radiative operator from
the explicit fallback to the implicit one, exactly as the He box runs it (commit
`e24f9932`): `rad_implicit_ang = true`, `rad_ang_solver = adi`, `rad_adi_scheme = lod2`,
`rad_ang_maxit = 400`, `rad_sts_perplane = true`, `rad_tr_tau_lo/hi = 0`,
`rad_blend_transverse = false`, `rad_gate_rho = 4.97e-10` (= `eos_rad_rho_lo`, the He-star
rescaling of the box's `6e-9`), and `rad_cap_ang` **removed** -- the code refuses both
treatments of the same operator.  The startup line confirms it runs:

    Conduction: transverse solver = ADI (lod2, theta = 1, sweep pairs/stage = 3),
                blocks per line = 1 (x2), 1 (x3)

`problem/rt_bottom_flux = false` was also added to the input, explicitly rather than by
default, so the 1-D relaxation can switch it on from the command line (AthenaK refuses a
command-line override of a parameter the file does not contain).

### Plane-parallel regression -- BITWISE IDENTICAL

`box_chk/`: `bench/hestar_fecz/box_w8/he_box_w8.athinput`, `nx2 = nx3 = 16`,
`meshblock 134x8x8`, 20 cycles, serial CPU, i.e. the same configuration as
`tests_m3/box_new`, run with a binary built from this merged tree **including the
`two_stream_rt` fix below**:

    cmp feczrt.hydro.hst ../../tests_m3/box_new/feczrt.hydro.hst   -> IDENTICAL
    cmp feczrt.user.hst  ../../tests_m3/box_new/feczrt.user.hst    -> IDENTICAL

(md5 pair in `box_chk/box_md5.txt`; the build directory was deleted afterwards for quota.)
So neither the merge nor the fix touches the He box production.

---

## STEP 2 -- the 1-D gate with the base flux injected  -- FAIL

Setup, following `bench/hestar_fecz/relax1d_m3/he_relax1d.athinput`: `inner_bc = wall`
(the box's `bc_mode = 3`), `problem/rt_bottom_flux = true` so that the internal flux
enters on the **two-stream's own lower boundary** rather than on the conduction wall face,
and `<hydro>/rad_flux_inner = L/(4 pi r_in^2) = 2.3066e38/(4 pi (1.18585e11)^2) =
1.305278e15` erg/cm^2/s.  `outer_bc = open` is kept, as are both sponges, mode 3
(`pcr`, `mixed 2`, `warm 1`), `rt_rad_force = true` and `rt_strang = true`.  The angular
grid is collapsed (`nx2 = nx3 = 4`, `meshblock 96x4x4`) and `vpert = 0`.  `sub_relax1d.sh`.
The plumbing is red_giant's own and needed no new code:

    red_giant: F_bot = 1.30528e+15 is carried by the TWO-STREAM's lower boundary
               (rt_bottom_flux); the conduction wall face injects nothing.

### 2a. The first run died at t = 53 s -- and it is a REAL cubed-sphere BUG

`relax1d_prefix/` (log `relax1d_prefix.log`).  Dead at cycle 5, `t = 53.1 s`, `dt` 2e-36,
`rt_use_cons gave a non-positive internal energy in 181 cell read(s)` -- the same death,
to a second, as the pre-merge `tests_1d` D(ii) run (cycle 5, `t = 56.6 s`).  Injecting the
base flux changed nothing, because **the failure is not at the base**.

The history file says what it is at a glance (columns 8/9/10 = 1-KE, 2-KE, 3-KE):

| t [s] | 1-KE (radial) | 2-KE | 3-KE |
| --- | --- | --- | --- |
| 12.93 | 3.27e33 | **2.77e20** | 2.77e20 |
| 25.99 | 2.06e34 | **2.21e30** | 2.20e30 |
| 39.30 | 1.21e35 | **3.11e40** | 3.11e40 |
| 52.83 | 9.49e43 | **3.96e57** | 3.36e57 |

On a spherically symmetric state with `vpert = 0` the transverse kinetic energy must stay
at round-off.  Instead it is amplified by **1e10 per cycle** (`dt` = 13 s) from the first
step, and the run dies when it overtakes the radial energy.  `diag/` (`sub_diag1d.sh`,
six 120-cycle arms) isolates the term; 2-KE at `t = 12.93 s` (one cycle):

| arm | switch | 2-KE at t = 12.93 s | outcome |
| --- | --- | --- | --- |
| `a_base` | -- | 2.77e20 | dead t = 53 s |
| `a_m0` | `rt_implicit_column = 0` | 1.10e20 | **dead t = 52 s** |
| `a_m3_plain` | `rt_impl_mixed = 0 rt_impl_warm = 0` | 1.99e20 | dead t = 54 s |
| `a_noforce` | `rt_rad_force = false` | **2.68e6** | alive at 120 cycles |
| `a_m0_noforce` | mode 0 + `rt_rad_force = false` | **2.18e5** | alive at 120 cycles |
| `a_noang` | `rad_implicit_ang = false` | -- | (arm invalid: `rad_ang_solver = adi` then fatals) |

Mode 0 dies exactly as mode 3 does, and both survive with the force off: **the driver is
`problem/rt_rad_force`, not the column solve.**  This retracts the `tests_m3` section-D
reading that mode 3 was driving the 1-D failure -- mode 0 only looked better there because
`rt_rad_force` was off in the arms that ran 500 cycles.

**The term and the cells.**  `src/utils/two_stream_rt.hpp`, the `Prad grad w` part of the
radiative momentum source, in the apply kernel:

```
f2 = cg*(rhoN(m,k,j+1,i) - rhoN(m,k,j-1,i))/(2.0*size.d_view(m).dx2);
f3 = cg*(rhoN(m,k+1,j,i) - rhoN(m,k-1,j,i))/(2.0*size.d_view(m).dx3);
```

`size.d_view(m).dx2` is the **coordinate** spacing.  On the cubed sphere `x2 = xi` and
`x3 = eta` are the gnomonic angles in `[-1,1]`, so at `nx2 = 4` this is `0.5`,
**dimensionless**, while the physical arc length is `pcoord->dx2(m,k,j,i) ~ r dxi ~ 1e11
cm`.  The transverse components of the radiative force were therefore too large by a
factor **~2e11**, every cycle, in every cell of the EOS taper (`w` strictly between 0 and
1, i.e. `rho` between `4.97e-10` and `1.93e-9`, radial indices i = 88..96, r = 2.33e11 ..
2.40e11) -- the cells the dt-collapse print names, `(m,k,j,i) = (1,5,4,92)`,
`r = 2.364e11`, `v = (4.6e6, 2.18e8, -1.71e8)`, transverse velocity 40x the radial one.
`conduction.cpp:1069` makes exactly this distinction with its `curv` flag; this file did
not.  The radial component of the same term (`X1V(i+1) - X1V(i-1)`) is a length and is
correct, which is why the plane-parallel box, where `size.dx2` IS a length, never saw it.

**Fixed** in commit `4bcdc855`: the transverse derivatives now divide by `DX2`/`DX3`, the
per-cell physical arc lengths on a curvilinear mesh and `size.dx2`/`dx3` on a Cartesian
one, so the box is bitwise (verified above).  No metric inverse enters: what `u0(IM2)`,
`u0(IM3)` hold on the cubed sphere is the COVARIANT angular momentum
(`coordinates.hpp:151-161`), and the covariant component of a gradient along the unit
tangent `e_xi` is the arc-length derivative itself.
**This bug is live for any curvilinear mesh, not just the cubed sphere** -- spherical
polar `dx2 = dtheta` is short by `r` in the same way -- so any `red_giant` /
`deep_hot_jupiter` run with `rt_rad_force = true` on `sp` carries it too.

### 2b. With the fix: 42x further, but STILL FAIL

`relax1d_fix/` (`sub_relax1d_fix.sh`, log `relax1d_fix.log`).  The taper blow-up is gone
and the transverse kinetic energy drops by **fifty orders of magnitude** at the same time
(2-KE at t = 52.8 s: 3.96e57 -> **9.59e7**, i.e. down to the `a_noforce` level of 1.1e8).
The run reaches `t = 2.22e3 s` = **0.47 turnover** instead of 53 s.  The gate table:

| gate | required | measured | verdict |
| --- | --- | --- | --- |
| no fatal / NaN | -- | dt collapse at cycle 221, t = 2220.6 s; never recovers | **FAIL** |
| 5 turnovers (2.35e4 s) | 2.35e4 s | 2.23e3 s (0.095 of it) | **FAIL** |
| `eos_fail` = 0 | 0 | not reached (run never completes) | -- |
| `dt >= 1 s` after 100 cycles | >= 1 | 13 -> 10 s to cycle 200, then 2.6e-2 .. 6e-4 | **FAIL** |
| mass drift | < 0.1 %/turnover | 1.91248151e26 -> 1.91247777e26 in 2121 s = **4.3e-6 %/turnover** | PASS |
| Newton `it_mean` | <= 4 | not printed before the collapse | -- |
| emergent L within 5 % of 2.3066e38 | 0.95-1.05 | `L_rad,out/L` 0.938 at cycle 0, 0.79 at t = 1.9e3, then diverges; `L_rad,cut/L` 0.868 -> **0.664** | **FAIL** |
| column re-stratifies smoothly | -- | see below | **FAIL** |

Column change between `t = 0` and the last healthy profile (`t = 2076 s`, shell means from
`relax1d_fix/rt_profile.bin.gz`, slots 0 and 5):

| layer | r [cm] | rho | T | mean v1 [cm/s] |
| --- | --- | --- | --- | --- |
| FeCZ base | 1.501e11 | +30 % | +7.7 % | -3.4e6 |
| kappa-peak region | 2.101e11 | **+120 %** | -5.6 % | -2.6e6 |
| FeCZ top | 2.292e11 | **-57 %** | -24 % | -3.2e6 |
| top of domain | 2.357e11 | -78 % | -28 % | -7.0e4 |

That is not a relaxation: the column is draining its upper half into the FeCZ at a mean
radial speed of 0.2-0.35 `v_MLT` and the deep radiative luminosity is falling steadily
(`L_rad,cut/L` 0.868 -> 0.664).  The first cell to collapse the timestep is **not** in the
taper any more but **deep in the FeCZ**: `(m,k,j,i) = (1,3,4,53)`, `r = 2.10e11`,
`rho = 1.43e-8`, radial `v1 = 8.8e9 cm/s` with transverse components of 600 cm/s, and a
garbage inverted temperature (the `T = 2.63e14` sentinel of the EOS-inversion failure).

**What this is.**  It is the failure `tests_m3` section D predicted, and injecting the
base flux did not remove it: 30-43 % of `L` is carried convectively through this column,
a 1-D column cannot carry it, so it must restratify toward the (much steeper) radiative
structure.  Giving the bottom face the full `L/(4 pi r_in^2)` bought 42x in time but did
not make the interior a radiative-equilibrium solution -- the FeCZ is heated from below
faster than it can radiate, the upper column falls in, and a cell inverts.  A second,
weaker transverse mode is also still there (2-KE e-folding ~35 s, the same rate with
`rt_rad_force` on or off, so it is NOT the force), but it is 1e-4 of the radial energy
when the run dies and is not what kills it.

**Where to go next (a setup decision, not a code fix -- left to the user).**
1. Relax the 1-D column against a **radiative-equilibrium** target instead of the MESA
   structure: rebuild the ic with `ic_tau_rad` pushed through the whole domain, so the
   1-D state the run is asked to hold is one it CAN hold, then let the 3-D run develop
   the convective flux.
2. Or drop the 1-D gate and go straight to 3-D with `inner_bc = open` (the deep-adiabat
   relaxation, which injects `L` as the entropy of the inflow and lets the convective
   flux exist), gating on `eos_fail`, `dt` and the kinetic-energy doubling instead of on
   a 1-D steady state.  `tests_1d` D(ii) shows the 1-D column has never been the thing
   this setup is good at; GAP_ANALYSIS G.5's "1-D column, 20 turnovers" test may simply
   be the wrong gate for a star whose envelope is convective over most of the domain.
3. Either way the `rt_rad_force` arc-length fix is a prerequisite and is now in.

---

## STEP 3 -- the 3-D smoke run  -- NOT RUN

Step 2's gate is red, so step 3 was not started (the brief: *"If the run fails, diagnose
to the level of tests_m3 (which term, which cells) and stop."*).  Nothing here changes
the 3-D input, which is ready: the smoke grid (`nx1 = 96`, `nx2 = nx3 = 32`,
`meshblock 96x16x16`, 24 blocks, 2 MI300A) and the seed (`vpert = 1e-3`, `v1`, tau window
20..300, `l ~ 2-6`) are already what `inputs/hydro/he4_presn_cs.athinput` carries.

## Files

* `sub_relax1d.sh` / `relax1d_prefix*` -- the gate as first run, dead at t = 53 s.
* `sub_diag1d.sh` / `diag/` / `diag1d.log` -- the six arms that localised `rt_rad_force`.
* `sub_relax1d_fix.sh` / `relax1d_fix*` -- the same gate with the arc-length fix.
* `box_chk/` -- the plane-parallel bitwise regression.
* `mk_ic_from_profile.py` -- writes a `problem/ic_profile` file from the last record of an
  `rt_profile.bin` (the relaxed column -> the 3-D initial state).  Written for step 2's
  hand-off, **unused**: the 1-D run never produced a relaxed column worth handing on.
  It needs no new pgen code -- `rt_profile`'s slots 0 and 6 are already `rho` and `eint`
  on the run's own radial grid -- and it splices the original ic file, rescaled for
  continuity, outside the active range so the radial ghosts are covered.
