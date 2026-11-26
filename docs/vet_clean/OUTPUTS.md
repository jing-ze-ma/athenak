# Outputs of a VET star run: history columns, energy check, warnings

This page describes what a production star run (problem generator `star_vet`, decks
`inputs/stars/star_bsg`, `star_m35mms`, `star_t9l62`) writes as history data and as run-time messages,
and how to check energy and mass conservation from the history file. The input keys are described in `KEYS.md`.

Units: the decks set `<units>` length, mass and time to 1, so code units are cgs (cm, g, s,
erg). The only exception is the temperature: the code temperature is T[K] /
`<rad_vet>/temp_unit_kelvin`. All sums are over the active cells of the wedge, not of the
full sphere. To convert a wedge luminosity to a full-sphere one, multiply by
4 pi / Omega, Omega = (cos x2min - cos x2max)(x3max - x3min).

## 1. History files

`<output1>` of the decks (`file_type = hst`, `dcycle = 1`, `data_format = %24.16e`) writes one
row per accepted cycle into two files in the run directory:

- `<basename>.hydro.hst`: the standard AthenaK hydro sums;
- `<basename>.user.hst`: the columns of the problem generator (`StarHist`,
  `src/pgen/star_vet.cpp`).

Facts that matter when you read them (`src/outputs/history.cpp`, `src/driver/driver.cpp`):

- Line 2 of each file is the header, `#  [1]=time  [2]=dt  [3]=...`. The numbers in the
  header are **1-based**; this page uses **0-based** column indices (what `numpy` uses).
  Labels are cut to 10 characters.
- A row is written after the cycle is accepted and before the next dt is computed, so the
  `dt` column is the dt of the step that **ended** at this row: `dt[i] = time[i] - time[i-1]`.
  A cycle that `dt_retry` takes again writes no row for the failed attempt.
- A restarted run appends to the same files and writes the header again. If the restart file
  is older than the last row, the time column jumps back. Skip that jump (the script in
  section 3 does).
- With `<mesh>/layout_invariant = true` the sums of `hydro.hst` and `user.hst` are done in a
  fixed order (per tile, the tiles' partials exact): the rows do not depend on the number of
  MPI ranks or on the MeshBlock layout made of whole tiles.

### 1.1 `<basename>.hydro.hst`

| col | name | meaning | unit |
|---|---|---|---|
| 0 | time | simulation time at the end of the step | s |
| 1 | dt | step that ended here | s |
| 2 | mass | sum rho V | g |
| 3 | 1-mom | sum rho v_r V (local basis; not a conserved quantity on spherical polar) | g cm/s |
| 4 | 2-mom | sum rho v_theta V (local basis) | g cm/s |
| 5 | 3-mom | sum rho v_phi V (local basis) | g cm/s |
| 6 | tot-E | sum of the conserved gas energy U_E V: internal + kinetic + rho Phi (with `etotgrav = true`). No radiation. | erg |
| 7 | 1-KE | sum (1/2) rho v_r^2 V | erg |
| 8 | 2-KE | sum (1/2) rho v_theta^2 V | erg |
| 9 | 3-KE | sum (1/2) rho v_phi^2 V | erg |

V is the true cell volume (also on curvilinear meshes).

### 1.2 `<basename>.user.hst`

The layout below holds when `<problem>/mlt_flux_frozen = false`, as in all three decks
(`star_vet.cpp`, `nb = 20`). With `mlt_flux_frozen = true` two columns `w_mlt`,
`L_MLT` are inserted at 22-23 and every later column moves up by 2. Read the header to be
safe.

Definitions used in the table:

- `is`, `ie`: the first and last active radial cell; `imid = is + nx1/2`. The decks have one
  MeshBlock along x1, so these are the first, last and middle cells of the mesh.
- `r_int`: the radius where the Rosseland optical depth of the **initial** column, integrated
  from the top of the mesh, reaches `<problem>/star_tau_int` (1.0 in the decks), i.e. the
  initial photosphere (`star_vet.cpp`). If it is never reached, `r_int = r_in`.
  "The L_int face" is the lower x1 face of the cell that contains `r_int`.
- `F`: the radiation flux stored on the x1 face. In the staggered scheme this is the
  face flux variable of the F equation (the radiation momentum times c^2, in flux units). In
  a steady state it equals the energy flux to O(v/c), but it is not the flux that moved E.
- `Phi`: the radiation **energy** flux that the radiation solve applied at a face in the step's
  last solve (the stage-2 solve of an sdirk2 step, or the backward-Euler solve). This is the
  flux that changed E; it includes the advective part (4/3) v E.
- `A`: the area of the x1 face. All L columns are sums over all (k, j) of the wedge.

| col | name | meaning | unit |
|---|---|---|---|
| 0 | time | simulation time | s |
| 1 | dt | step that ended here | s |
| 2 | L_bot | at the x1 face between the first and second active cells (`is+1`), one cell above the inner boundary: with `scheme = clean` the comoving flux sum (Phi - (4/3) v E) A (as L_midC; in earlier versions it was the face flux sum F A); with the old scheme sum F A | erg/s |
| 3 | L_mid | sum F A at the lower face of cell `imid` | erg/s |
| 4 | L_int | sum F A at the L_int face | erg/s |
| 5 | L_top | energy that the step's radiation solves moved **out** through the outer x1 face, divided by dt (positive outward). This is the boundary record of the radiation solve (`BoundaryEnergyAcc`, `src/rad_vet/rad_vet_imex.cpp`). If the record is missing for a step it falls back to sum F A at the top face. | erg/s |
| 6 | L_in | energy that the step's radiation solves moved **in** through the inner x1 face, divided by dt (positive inward). Fallback: `implicit_flux_x1min` x A(r_in) summed over the inner faces. | erg/s |
| 7 | E_rad | sum E V | erg |
| 8 | e_gas | sum e V, e = gas internal energy density | erg |
| 9 | M_int | sum rho V over cells with centre r <= r_int | g |
| 10 | KE_int | sum (1/2) rho v^2 V, r <= r_int | erg |
| 11 | KEr_int | sum (1/2) rho v_r^2 V, r <= r_int | erg |
| 12 | Mr_int | sum rho v_r V, r <= r_int | g cm/s |
| 13 | PV_int | sum P V, r <= r_int | erg |
| 14 | V_int | sum V, r <= r_int | cm^3 |
| 15 | M_tot | sum rho V, whole mesh | g |
| 16 | Mdot_top | outflow part of the hydro mass flux at the outer face, sum max(f_rho, 0) A (flux of the last hydro stage, not the step average) | g/s |
| 17 | Mdot_bot | hydro mass flux at the inner face, sum f_rho A (positive = upward, into the mesh); 0 with the closed wall (`StarWallFlux` sets its flux to 0 before the update) | g/s |
| 18 | Etot | sum (U_E + E) V: gas total energy (internal + kinetic + rho Phi_grav) plus radiation energy. With closed walls its rate is L_in - L_top. | erg |
| 19 | Min_top | inflow part at the outer face, sum min(f_rho, 0) A (<= 0): the mass flux that `StarWallFlux` removed at its call, before `dfloor_flux_limit`. 0 also where a face without mass flux passes L_z and E (see the note below) | g/s |
| 20 | v1sq_wall | mean v_r^2 over the first active radial cells of the whole mesh | cm^2/s^2 |
| 21 | Picard | **Always 0** under `scheme = clean`; the radiation solve reports its passes in the `CLEAN cycle` lines (section 4.2). | - |
| 22 | P_esrc | power that the explicit radiation source term puts into E, sum (chat/c) esrc V, at the state of this row. 0 in the decks (no frozen MLT deposit, no test source). | erg/s |
| 23 | P_sponge | rate of kinetic energy loss of the top sponge at the state of this row. 0 in the decks (`star_sponge_rate` = 0). | erg/s |
| 24 | M_hflx | hydro mass leaving through the physical x1 faces (inner + outer, outward positive) over the step, divided by dt. These are the RK-weighted fluxes the update used (`AmBndFlux`, `src/hydro/hydro_am.cpp:637`), so they close the mass budget exactly. Needs the angular-momentum module (`<hydro>/am_transport`, on in the decks); 0 without it. | g/s |
| 25 | E_hflx | the same for the hydro conserved energy U_E | erg/s |
| 26 | Lz_hflx | the same for the gas L_z | g cm^2/s^2 |
| 27 | Lz_rflx | the radiation's L_z leaving through the physical x1 faces over the step, as a rate: the azimuthal faces' transverse stress torque (D_{r phi} E) at those faces, stage-weighted as the update (LzBook's rbt); a boundary flux, not ledger | g cm^2/s^2 |
| 28 | Lz_tot | L_z of gas + radiation, including the radiation's L_z ledger (`LzTotal`); summed exactly, so dLz_tot + (Lz_hflx + Lz_rflx) dt of a row (the rate of the step it ends) holds to the state's round-off, a small part of eps = 1e-16 times the gross sum V R \|p_phi\| per step: where the net is far below the gross sum, compare with eps x gross, not with eps x Lz_tot | g cm^2/s |
| 29 | L_midE | sum Phi A at the L_mid face: the energy flux that moved E there | erg/s |
| 30 | L_intE | sum Phi A at the L_int face | erg/s |
| 31 | L_midC | L_midE minus the advective part: sum (Phi - (4/3) v E) A at the L_mid face, with v and E the means of the two cells beside the face. This is the comoving radiative flux. | erg/s |
| 32 | L_intC | the same at the L_int face | erg/s |

Notes:

- L_mid, L_int (F) and L_midE, L_intE (Phi) agree in a steady, slowly moving atmosphere and
  differ in a dynamic one; L_midC/L_intC remove the advection of radiation by the gas.
  L_midE..L_intC are 0 on rows where no stored solve flux exists yet (the first row).
- L_mid is a sum of fluxes of both signs over a shell; in a convective interior it can be a
  small difference of large numbers, so small CPU/GPU or layout differences in it are not a
  sign of a problem.
- Mdot_top and Min_top use the last hydro stage's flux; for a mass budget use M_hflx.
- With `hydro/p_ledger = true` ten columns follow (the `am_hist` history of the test
  decks too): `Pr_tot`, gas + radiation radial momentum (sum V m_r plus each cell's halves of
  its x1 faces' V F / c^2), summed exactly, g cm/s; then the step's rates (g cm/s^2):
  `Pr_hflx`, the hydro momentum leaving through the physical faces (outward +); `Pr_ntel`,
  the flux divergence's part that does not telescope; `Pr_src`, the coordinate, gravity,
  am and user sources, less the reference force the radiation withholds from the gas
  (force_reference: the hydro's well-balanced background gives it); `Pr_rad`, the
  radiation update and write-back (boundary radiation stress, geometric terms); `Pr_cut`,
  the cone cut the gas does not take; `Pr_fwb`, the fwb residual's momentum;
  `Pr_flr`, floors
  (dfloor_flux_limit, ConToPrim, floor-gas cap); `Pr_mrg`, the average of the merged
  lateral cells (`merge_lat`); `Pr_gap`, the rest (0). So dPr_tot + Pr_hflx dt - (the
  other rates) dt = 0 to the round-off of the sums (the `END` of `CLEAN PLEDGER`); `ntel`,
  `mrg` (`merge_lat_mom = radial`) and `gap` are round-off. On a
  Cartesian mesh the columns are `Px_*`, `Py_*`, `Pz_*`.
- With `hydro/p_ledger_lz = true` (3-D spherical polar) twelve `LzL_*` columns take the
  place of the `Pr_*` ones: `LzL_tot`, gas + radiation L_z (sum V R m_phi plus each cell's
  halves of its x3 faces' V R F / c^2: `Lz_tot` less the L_z ledger, with the face F), g
  cm^2/s; the step's rates (g cm^2/s^2) of the same
  classes (`LzL_hflx`: the gas's torque through the physical x1 faces, `Lz_hflx`; cut, fwb
  and the reference force stay in `LzL_rad`); `LzL_lzb`, the L_z ledger's growth (`CLEAN
  LZB`: cone cut, update_mcap cap, reference torque, lz_spin's return, what does not
  telescope), and `LzL_rbt`, the radiation's torque through the physical x1 faces (the
  off-diagonal stress D_{r phi} E of the x3 faces' torque form; `Lz_rflx` x dt), booked
  outside the ledger. So `LzL_rad` + `LzL_lzb` + `LzL_rbt` is the write-back's L_z that
  neither holds (round-off).
  `CLEAN LZLEDGER` prints the same per step with `END`.

## 2. Event log `<basename>.log`

`<output4>` (`file_type = log`) writes counters since the previous row (and, for the
`vet_*` columns, totals since the start of the run or restart): `cycle`, `eos_dfloor` (cells
where the density floor acted), `eos_efloor`, `eos_tfloor`, `eos_vceil`, `eos_fail`, `c2p_it`,
`fofc`, `efloor_de` (internal energy created by the energy floor), `eos_tclamp`, `eos_tset`,
`vceil_de` (kinetic energy the velocity ceiling removed), `dfloor_ke` (kinetic energy change of
`dfloor_keep_velocity`), then the radiation totals `vet_bcg_break vet_bcg_fall vet_bcg_keep
vet_newt_fb vet_gn_switch vet_od_fall vet_vlat_fall vet_opn_skip vet_pos_flr vet_pos_un vet_pos_s2
vet_acl_e vet_acl_f` (`src/outputs/eventlog.cpp`). For this scheme the radiation columns
that matter are `vet_pos_flr` (E floor raises) and `vet_pos_un` (energy the E floor created,
x volume); most of the others are not used by this scheme.

## 3. Checking energy and mass conservation from `<basename>.user.hst`

Take a step as the pair of rows i-1 -> i. The rates are those of row i (they describe the
step that ended there) and dt = time[i] - time[i-1] (= the dt column of row i). With the
0-based columns of section 1.2 (t 0, L_top 5, L_in 6, M_tot 15, Etot 18, P_esrc 22,
P_sponge 23, M_hflx 24, E_hflx 25):

    rE_i = (Etot_i - Etot_{i-1})
           - (L_in_i - L_top_i + P_esrc_i - P_sponge_i - E_hflx_i) dt_i
           - dE_floor_i

    rM_i = (M_tot_i - M_tot_{i-1}) + M_hflx_i dt_i - dM_floor_i

Normalise rE by L_in dt (the energy that entered through the base in the step) and rM by
M_tot. A conservative step gives rE/(L_in dt) at round-off (about 1e-9 or less) and rM/M_tot
of about 1e-15.

**The density floor term.** `dE_floor` and `dM_floor` are the energy and mass that the hydro
floors created in the step. They are **not** in the history file.

- The mass that `<hydro>/dfloor` creates adds to M_tot, and each gram adds its potential
  energy Phi_grav = GM (1/r_in - 1/r) >= 0 to Etot (`etotgrav = true`; star_vet sets the
  zero of Phi at r_in, `src/pgen/star_vet.cpp`). `<hydro>/floor_pay = rad_column`
  pays the **internal** energy the floors create from the radiation, so that part does not
  change Etot (what cannot be paid is created and printed as a CLEAN WARNING), but nothing
  pays the potential energy of the new mass.
- So if the density floor is active, rE is not at round-off. Because the floor acts near the
  top, where Phi_grav is largest, the residual can be large compared with round-off: in the
  bsg deck the floor mass gives rE of order 1e-4 L_in dt per cycle at about 1.3e15 erg per
  gram of floor mass (measured), while m35mms stays at about 1e-9.
- A working estimate: on a step with no other mass source, rM (computed without the floor
  term) **is** the floor mass. Then dE_floor is about rM x Phi_grav(r) of the floored cells;
  near the top Phi_grav ~ GM (1/r_in - 1/r_top), with GM = `<problem>/star_gm`. Use the event
  log `eos_dfloor` column to see in which intervals the floor acted.
- The run totals of the floors are printed at the end of the run (section 4.3, line
  `<hydro> EOS floors ...`: mass created x volume and internal energy created x volume).
- Other terms not in the history, which show up in rE when they are non-zero: the kinetic
  energy the velocity ceiling removes (`vceil_de` in the event log), the kinetic energy change
  of `dfloor_keep_velocity` (`dfloor_ke`), unpaid floor energy (CLEAN WARNING), and energy
  created by the radiation E floor (`energy created` in the `vet-positivity` line,
  `vet_pos_un` in the event log).
- `star_wall_zero_flux` (`StarWallFlux`) closes the top faces whose mass flux is inflow or
  exactly 0 when it is called, after the Riemann and FOFC fluxes. `dfloor_flux_limit` runs
  later, in `Hydro::RKUpdate`: a top face that it makes (nearly) massless stays open to the
  diffusive L_z and E of FOFC's flux. That flux is booked in E_hflx and Lz_hflx and conserved,
  but it is not carried by gas. In bsg (100 cycles, Caltech) Min_top is 0, Lz_hflx dt is
  1.96e35 per cycle and the L_z residual 1.57e31 rms per cycle (round-off); this history is
  bitwise that of `star_wall_zero_flux = false`.
- Earlier versions took the removed inflow back after the update, and rE in bsg was a few
  1e-6 L_in dt per cycle once gas fell in through the top. The flux is now 0 before the
  update; rE with inflow was not measured again on bsg.

Python (numpy) script that prints the residuals:

```python
import re, sys
import numpy as np

fn = sys.argv[1]                                   # e.g. bsg3d.user.hst
with open(fn) as f:
    f.readline()
    hdr = f.readline()                             # "#  [1]=time  [2]=dt  [3]=L_bot ..."
col = {k: int(i) - 1 for i, k in re.findall(r'\[(\d+)\]=(\S+)', hdr)}   # 0-based
d = np.loadtxt(fn, ndmin=2)                        # '#' lines (also repeated headers) skipped
g = lambda k: d[:, col[k]] if k in col else np.zeros(len(d))
t = g('time')
dt = np.diff(t)                                    # = the dt column of rows 1..n
ok = dt > 0                                        # skip the jump where a restart appends
rate = g('L_in') - g('L_top') + g('P_esrc') - g('P_sponge') - g('E_hflx')
rE = np.diff(g('Etot')) - rate[1:]*dt              # subtract dE_floor here if you have it
rM = np.diff(g('M_tot')) + g('M_hflx')[1:]*dt      # subtract dM_floor here if you have it
nE = rE/(g('L_in')[1:]*dt)                         # energy residual in units of L_in dt
for i in np.nonzero(ok)[0]:
    print('t %14.6e  rE/(L_in dt) %11.3e  rM/M_tot %11.3e' % (t[i+1], nE[i], rM[i]/g('M_tot')[i+1]))
print('max |rE|/(L_in dt) = %.3e   sum rE / sum(L_in dt) = %.3e' %
      (np.abs(nE[ok]).max(), rE[ok].sum()/np.sum((g('L_in')[1:]*dt)[ok])))
```

Run it as `python3 hst_balance.py bsg3d.user.hst`. It finds the columns by their header
labels, so it also works with `mlt_flux_frozen = true`. A positive rE means energy appeared
that the fluxes do not explain.

## 4. Run-time messages

All messages are printed by MPI rank 0 to standard output.

### 4.1 `CLEAN WARNING` messages

A production run should print **no** CLEAN WARNING. Each one marks an event the code
handled but that the user should know about.

| message (start) | where | meaning | what to do |
|---|---|---|---|
| `CLEAN WARNING cycle N: backward-Euler solve not converged (dt_retry, K so far; last change X)` | `rad_vet_newton.cpp` (`VerdictSolve`) | A backward-Euler radiation solve (first step, the redo after a failed stage, or an auto-BE step) hit `newton_maxit` passes without converging. With dt_retry on (production) the cycle is taken again with a smaller dt; the text says `no dt_retry: accepted` or `no dt_retry: > implicit_resid_fatal` (then the run stops) when dt_retry is off. A `CLEAN WORST` line names the worst cell. | Occasional ones that dt_retry cures are tolerated by the scheme but are not normal; look at the `CLEAN WORST` cell (position, E, T, opacity). Repeated ones mean the step is too large for the solver there: check the keys against the production set (`KEYS.md`, restart rule), then the state at that cell. |
| `CLEAN WARNING cycle N t2s S: stage solve failed, the step is redone with backward Euler (stage fallback K)` | `rad_vet_imex.cpp` (`Time2Restore`) | A stage of the sdirk2 step was not admissible (not converged or not positive); the step is redone with backward Euler (first order for that step). `t2s`: 1 = backward Euler with stored slope, 2 = stage 1, 3 = stage 2. `time2_auto_be` then skips the stage attempt for a few steps. | Rare fallbacks are part of the design; frequent ones lose second-order accuracy in time and cost time. Report them with the deck and the cycle. |
| `CLEAN WARNING cycle N t2s S: stage start not finite in K cells (gas or E; the solve fails and the step is redone)` | `rad_vet_oldvector.cpp` (`OldVector`) | The state entering a stage has NaN/Inf. The stage fails and the step is redone. | Serious. If it repeats, stop the run and find the source (a hydro failure, a floor, a corrupted restart). |
| `CLEAN WARNING cycle N: K cells with a ray J not finite (isotropic)` | `rad_vet_closure_rays.cpp` (`ClosureFromRays`) | The ray sweep gave a non-finite mean intensity in K cells; those cells are given an isotropic closure for this build. | Should not happen. Check the opacities and the temperature in those cells; report it. |
| `CLEAN WARNING cycle N: K ray intensities < 0 set to 0` | `rad_vet_closure_rays.cpp` (`ClosureFromRays`) | Cartesian sc sweep: K ray intensities came out negative and were set to 0. None by construction (the upwind interpolation is bounded, the segment weights are >= 0). | Should not happen; report it with the deck. |
| `CLEAN WARNING clean_init steady: not converged (residual X, tol Y)` | `rad_vet_steady.cpp` (`SteadyStart`) | The steady start (`steady_init = steady`) ended before reaching `steady_tol`. If the last change is above 1e-2 the run stops with a FATAL (diverged start). | Small residuals are usually harmless; otherwise raise `steady_n` or change `steady_dtfac`, or start from a state closer to radiative equilibrium. |
| `CLEAN WARNING cycle N: floor_pay = rad_column left X of the fluid floors' energy unpaid (created)` | `src/hydro/dual_energy.cpp:419` | The radiation of the cell's column could not pay all the internal energy the hydro floors created; X (energy x volume) was created. | Energy is not conserved by X in this cycle (see section 3). If it repeats, the floors act where there is too little radiation: check `dfloor`, `pfloor`, `tfloor`. |
| `CLEAN WARNING cycle N: Ebnd without the step's boundary energy (BoundaryEnergyAcc), its faces take 0` | `rad_vet_ledger.cpp` (`EbndAccept`) | Only with the boundary-energy ledger of the drive/test problem generators (not in production). | None in production. |
| `CLEAN WARNING restart: first dt X != the file's dt Y` | `src/driver/driver.cpp:666` | The first cycle after a restart does not use the dt stored in the restart file, so the restart is not a bitwise continuation. | Check that the restart deck does not change keys that set dt (`cfl_number`, `dt_retry` state, `time2_cfl_post`, `rad_signal_speed`). |
| `CLEAN WARNING restart dump <file>: no progress for N s after ...` | `src/outputs/restart.cpp:739` | Writing a restart file has made no progress for `<job>/restart_warn_s` seconds (default 600). | A file system problem: check the file system and the job; the run may hang. |

### 4.2 Other per-cycle lines of the scheme

With `<rad_vet>/verbose = 0` (default) the code prints the `CLEAN cycle` lines, the
`CLEAN WARNING` lines and the event lines (`CLEAN WORST`, `GNSW`, `KPOS`
without a face slope). `verbose = 1` adds the per-solve diagnostics (`CLEAN MCAP`,
`OLDVEC`, `FWB`, `FCONE`, `LZ`/`LZB`, `WBFCT`, `KPOS`, `KREAL`, `XPAT`, the `XT`
write-back closure, the ray cut). `verbose = 2` adds a `CLEAN PASS cycle N t2s S
pass P it K de D dv V tol T fe F fv G [clos] [opac]` line for every pass of every solve:
its Krylov iterations, its changes against the tolerance the convergence test uses, and the
round-off levels of that test (to see which passes a solve needs and why).
`CLEAN WARNING` lines always print.

| line | meaning | what to look at |
|---|---|---|
| `CLEAN cycle N t2s S it K (sum M) pos 0 mark 0 all 0 neg Q r1 . r2 . first . . guard . solves . cwin .` | One line per radiation solve. `it`: Newton/Picard passes of this solve (**negative = not converged**); `sum`: linear iterations; `neg`: cells with E < 0 after the solve; `cwin`: linear solves that stopped on the stall window. | `it` < 0 and `neg` > 0 should be rare. |
| `CLEAN WORST cycle N ... gid g k j i E E_old T S chi_T` | Printed with a not-converged backward-Euler solve: the cell with the largest last change (`BREAKDOWN` if not finite). | Where and why a solve fails. |
| `CLEAN GNSW cycle N t2s S pass P: stalled, the gas Newton slope off (K so far)` | A stalled solve continues on the bracketed gas root. | Informational; many per cycle point to a hard region. |
| `CLEAN OLDVEC cycle N t2s S eint<0 K (carried) E<0 Q` | The state that starts a stage has K cells with gas internal energy < 0 (carried into the gas row) and Q cells with E < 0. | Occasional values are expected in floor gas; growing counts are a warning sign. |
| `CLEAN FWB cycle N t2s S faces K radiation momentum sum \|dp\| X x1 faces sum dp Y` | `clean_stag_fwb`: the radiation momentum of the well-balancing residuals of this solve. | Diagnostic. |
| `CLEAN FCONE cycle N t2s S faces K momentum cut sum \|dF\| V/c^2 X` | `flux_cone` scaled K faces back to \|F\| <= cE; X is the momentum removed. | Should stay small compared with the radiation momentum. |
| `CLEAN WBFCT cycle N t2s S rounds R left K kept X` | The limiter rounds of the write-back; `left` = cells still limited after the last round. | `left` > 0 means the limiter could not fully protect E > e_floor in those cells. |
| `CLEAN MCAP cycle N t2s S: K face momentum terms at the cap h chi_f E in the last pass` | Face momentum exchange terms at their cap. | Diagnostic. |
| `CLEAN LZB cycle N t2s S ledger X` | 3-D spherical polar: the radiation's L_z ledger (stored in the restart file), without the torque through the physical x1 faces, which is a flux (`Lz_rflx`). | Use with `Lz_tot` in the history. |
| `CLEAN FLOORMOM cycle N: ekin K dp X Y Z \|dp\| X Y Z vceil ... other ...` | The momentum the step's primitive recoveries removed, with their end-of-step weights: `ekin` (`floor_from_ekin`), `vceil`, `other` (density floor, or with a raised E); K cell-conversions, then the signed and absolute sums x volume of (r, theta, phi), Cartesian (x, y, z). Exact sums (any layout / rank count). `initial`: the first conversion; the run end prints the job's totals (not stored in the restart file: the lines of the jobs add up). | Compare with the total momentum; large values mean the floors drive the flow. |
| `CLEAN PLEDGER cycle N dt T hflx .. ntel .. src .. rad .. cut .. fwb .. flr .. mrg .. gap .. P .. END .. gross G` | `hydro/p_ledger`: the step's momentum books by class as rates (one value per component: r on spherical polar / cubed sphere, x y z on Cartesian), P = gas + radiation momentum after the step, END = dP + hflx dt - (the other classes) dt, G = sum V \|p\| (gas + radiation). | END / G is round-off (1e-16); `ntel`, `fwb`, `mrg` and `gap` should be round-off too. |
| `CLEAN KPOS`, `CLEAN KREAL`, `CLEAN XPAT` | Only with those (non-production) keys. | - |
| `DT_RETRY cycle=N t=T dt A -> B reason=R resid=X cut K` | A cycle is taken again with a smaller dt. `reason`: NON-CONVERGED, resid>implicit_resid_fatal, eint<=0, vkick, floor_pay_unpaid, forced. | A few per run are normal for dynamic stars; many in a row mean the run is in trouble. |
| `CLEAN RETRY cycle N t2s S: the redone step's solve converged (last change X, tolerance T)` | A solve of a cycle that dt_retry takes again converged. The redo is held to the tolerance of any step (the time-error budget at the cut dt), not a stricter one: a stalled step cannot be rescued by a cut dt alone otherwise. | - |
| `DT_RETRY recovered cycle=N ...` | The dt cap that dt_retry imposed is lifted. | - |
| `<rad_vet> floorcap cycle N (cumulative, all calls): n_v . n_t . dKE . dEint . (erg, to E_rad) momentum removed \|dp\| . dL_z .` | Every `floorcap_report` cycles: the floor-gas cap totals. The energy is moved to E_rad (conserved); the momentum is removed. | The removed momentum and L_z are a known non-conservation of the cap; keep them small. |

### 4.3 End-of-run summary lines

At the end of the run (`RadiationVet` destructor and `Driver::Finalize`) rank 0 prints a
block of `<rad_vet>`, `<time>` and `<hydro>` lines. The ones that matter for this scheme:

| line | what it tells you |
|---|---|
| `<rad_vet> time_scheme=sdirk2: stage steps=A backward-Euler steps=B stage fallbacks=C ...` | Steps taken with the two stages (A), with backward Euler (B; the first step, auto-BE steps), and the stage fallbacks (C; each one also printed a CLEAN WARNING). A large B or C means the run was often first order in time. The `face_fct` part gives the stage limiter's totals, including `left (Q<=0)` and `emptied` cells. |
| `<rad_vet> vet-positivity: gas cell-solves=K energy moved rad->gas=X \| E floor raises cell-solves=K energy from gas=X energy created=Y ... taken from the column's radiation (implicit_pos_floor_cons)=Z \| hydro EOS floors energy created=W` | The positivity book-keeping of the write-back. **`energy created` (Y) should be 0**: a non-zero value is energy the E floor created (energy x volume). `energy from gas` and `taken from the column's radiation` are moves (conserved). The last field repeats the hydro floors' internal energy. |
| `<hydro> EOS floors (all ranks, accepted steps, RK-weighted): mass created x vol=M energy created x vol=E \| FOFC first-order cell-stages=K` | The run total of the mass and internal energy the hydro floors created (cgs, since the start or restart). M is the density-floor mass of section 3 (its potential energy is not in E). |
| `<time> dt_retry: N checkpoints (...), R retries, cap X` and `<time> dt_retry: K cycles took the backward-Euler redo` | How many cycles were retried, and the dt cap at the end. |
| `<rad_vet> implicit energy solve: cells=... FAILURES=...` | The explicit-transport energy solve; not used by this scheme. |

### 4.4 Start-up messages and FATAL errors

- `### WARNING (scheme = clean): <rad_vet/key> has no effect, the clean solve reads <rad_vet>/other`
  and `### WARNING (scheme = clean): no effect: ...`: the input names a key that the
  solve does not use; use the key named in the message. (`CheckKeys`,
  `src/rad_vet/rad_vet_tasks.cpp`.)
- `### NOTE: <rad_vet>/implicit_maxit = A is not the clean solve's: newton_maxit = B`: that
  key is ignored; name `newton_maxit` to change the limit.
- `### WARNING: N input parameter(s) were never read by the code`: printed after the first
  cycle; check each listed key for a typo (keys that the templates name without effect, `KEYS.md` appendix, may be
  listed).
- `### FATAL ERROR ... <rad_vet/key> belonged to a removed path` / `the key was removed, the code
  runs V (name key=V to continue so)` / `has no counterpart in scheme = clean`: the input names a removed
  key or a removed value; remove it or set it to the accepted value. A restart file of an earlier
  run may carry such a key: override it on the command line (for example
  `rad_vet/time2_stage_local_be=false`).
- `<rad_vet>/scheme = clean currently needs a uniform mesh and chat = c`, `... on spherical
  polar needs ... the rays (closure = sc)`, `force_reference = wb_arad(_md) ... needs
  ... force_reference_work = split`, `clean_stag_fschur cannot take this mesh (nghost >= 3,
  MeshBlock size, uniform mesh)`, `<mesh>/layout_tile N needs ...`:
  set-up checks; the message says which key or mesh property to change.
- `### FATAL ERROR in src/main.cpp: a key that places the cells differs from the restart
  file`: a restart deck changes the grid; restart with the grid the file was written with.
- `clean_init steady: the last pass changed E or F/c by X of max E (diverged: steady_dtfac, _n)`:
  the steady start diverged (last change above 1e-2); change `steady_dtfac` or
  `steady_n`, or the initial state.
- `<rad_vet> clean: a backward-Euler solve did not converge and its last change exceeds
  <rad_vet>/implicit_resid_fatal`: only without dt_retry. Turn dt_retry on (production).
- `### FATAL ERROR: <time> dt_retry: cycle N failed K times (max 8), dt X`: eight cuts did not
  save the cycle; the run stops with the last restart file intact. Look at the CLEAN
  WARNING / CLEAN WORST lines before it.

## 5. Caveats

- The column meanings come from `StarHist` and the functions it reads. Columns can be
  added in later versions; always read the header.
- `Phi` in L_midE/L_intE is the flux of the step's last solve as stored for the next step's
  first slope (`t2pf1`). This follows the code comment and the write path in
  `WriteBack` (`src/rad_vet/rad_vet_update.cpp`); it was not checked against a dump.
- That `E_hflx` contains the potential-energy flux rho v Phi_grav (so that Etot closes with
  etotgrav) follows from the conserved-energy flux the hydro uses; it was not traced line by
  line.
- The size of the density-floor residual quoted in section 3 comes from a measurement for the bsg and m35mms decks; other decks and run phases differ.
