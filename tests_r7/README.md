# he4_presn: THE 0.79 R KILLER IS THE RADIATIVE FORCE, NOT THE FeCZ -- AND THE MLT
# CLOSURE NOW SURVIVES A RESTART

2026-09-17, viper, branch `he4-presn-global`, against `5b71b144` (= the `tests_r6`
commit).  Read `tests_r6/README.md` first: this round executes its section-5 items 1
(diagnose the 0.79 R event) and 2 (put `fmlt1d` in the restart file).

Two binaries.  `tests_r4/athena_v2` = `5b71b144` unchanged, used for every arm of part A
(so the arms are directly comparable with `tests_r6`); `tests_r7/athena_v3` = this
round's commit, used for part B.  All GPU (HIP, `Kokkos_ARCH_AMD_GFX942_APU`, MPI on),
`apudev`.  Every `.bin`, `.rst` and `.cbin` has been deleted after measuring.
Turnover = 4705 s.

---

## 0. THE 1-D COLUMN, AND WHAT IT IS

`one1d.sh`, `bis1d.sh`, `seam1d.sh`, `chk.sh`.  The configuration is arm `Gnr` of
`tests_r6` -- `problem/inner_bc=wall problem/rt_bottom_flux=true
hydro/rad_flux_inner=1.305278e15 problem/mlt_alpha=1.5`, with the input file's own
`dfloor_keep_velocity`, `vceil = 1e8` and `efloor_as_tfloor` -- with the ANGULAR grid
collapsed and the seed off:

    mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0

One 96x4x4 block per panel, 6 blocks, 2304 cells.  It is "1-D" only in the sense that the
state is spherically symmetric: the gnomonic metric is NOT constant over a panel, so the
96 angular cells are not discretely identical, and with `nx2 = 4` and `nghost = 3` EVERY
angular cell is in the panel's outermost row (`js = ks = 3`, `je = ke = 6`).  That matters
in section 2.

## 1. REPRODUCTION: YES, IN 1-D, AT THE SAME TIME  (arm `A1`)

| | 3-D `Gnr`/`Gv5`/`Hnr`/`Dnr` (tests_r6) | 1-D `A1` |
| --- | --- | --- |
| dies at | t = 4879-5074 s = **1.04-1.08 turnover** | t = 5602 s = **1.19 turnover** |
| dt at the collapse | 10.7 -> collapse | 10.5 (t = 4872) -> 3.6 (5525) -> 0.015 (5616) |
| the cell | (11,18,18,34) .. (7,3,4,36), r/R **0.776-0.791** | (1,6,3,98), r/R **1.013** |
| its state | `v1` = `vceil`, T = 1.3-2.3e6 K | `v1` = 9.999980e7 = `vceil`, T = 3.16e6 K |

**The event is neither turbulent nor three-dimensional**: a spherically symmetric column
with no seed at all dies within 10 % of the same time, with the same signature -- ONE cell
with `|v1|` pinned exactly on `vceil` and a temperature of 1-3 MK against an ambient
1.7e5 K.  The radius moves (the 1-D arms mostly die at the top, the 3-D arms at 0.79 R);
section 2 shows why that is not a second event.

### 1.1 THE SHELL MEAN NEVER MOVES -- IN EITHER RUN

`prof.py <rt_profile.bin>` reads the area-weighted shell means.  The 3-D `Gnr` profile and
the 1-D `A1` profile agree to 3-4 digits at every time, and the shell-mean T at 0.79 R is
FLAT in both right up to the collapse:

| t/turn | 0.00 | 0.25 | 0.50 | 0.75 | 1.00 | 1.17 |
| --- | --- | --- | --- | --- | --- | --- |
| `A1`  T(0.791 R) [code] | 1.358e13 | 1.361e13 | 1.361e13 | 1.357e13 | 1.353e13 | 1.380e13 |
| `Gnr` T(0.791 R) [code] | 1.358e13 | 1.361e13 | 1.361e13 | 1.357e13 | 1.353e13 | -- |

(code T / 8.314e7 = T[K]: 1.353e13 = 1.627e5 K, 0.4 % below the IC's 1.634e5 K.)
The hottest shell in the domain is the innermost one at every dump in every arm.

**There is no thermal front and no runaway of any shell-mean quantity.**  The collapse is
a SINGLE CELL, four to five orders of magnitude in T above its own shell.  Two consequences
that answer `tests_r6` section 5 item 1 directly: a shell-mean closure cannot see this
event, and it also cannot be causing it, because nothing the shell mean carries is
changing.

## 2. WHERE: THE PANEL-SEAM ROW, 11 TIMES OUT OF 11

The `(m,k,j,i)` of every `### dt COLLAPSE` print of the four `tests_r6` no-restart arms
(`nx2 = nx3 = 16`, so `js = ks = 3`, `je = ke = 18`, and the outermost angular row is
`k` or `j` in {3, 18} = 60 of 256 cells per panel = 23.4 %):

| arm | (k, j) | on the seam row? |
| --- | --- | --- |
| `Dnr` | (18, 17) | yes (k = ke) |
| `Gnr` | (18, 18), (5, 18) | yes (cube VERTEX), yes |
| `Gv5` | (18, 6), (3, 3), (18, 5), (4, 3) | yes, yes (VERTEX), yes, yes |
| `Hnr` | (3, 4), (15, 18), (18, 16), (3, 5) | yes, yes, yes, yes |

**11 of 11.**  By chance that is 0.234^11 = 1.4e-7.  Two of the eleven are at a cube
vertex.  The radial index is equally concentrated: `i` = 33-36 (r/R 0.77-0.79) or `i` =
92-98 (r/R 0.997-1.013), never in between.

This is what makes the 1-D column the right reproducer AND why it dies at a different
radius: at `nx2 = 4` all four angular cells are seam cells, so the whole 1-D run is a
seam experiment and the event is free to pick whichever radius is marginal first.

## 3. THE BISECTION  (all 1-D, `athena_v2`, `tlim` as stated)

| arm | switch | outcome | t at death | turnover |
| --- | --- | --- | --- | --- |
| `A1` / `x_base` | -- | dies, (1,6,3,98) r/R 1.013 | 5602 | 1.19 |
| `b_mlt0` | `problem/mlt_alpha=0` | dies, (2,6,6,92) r/R 0.997 | 5728 | **1.22** |
| `g_percol` | `problem/mlt_mean=false` | dies, (2,3,4,98) r/R 1.013 | 4789 | **1.02, WORSE** |
| `x_nostrang` | `problem/rt_strang=false` | dies, (0,3,3,95) r/R 1.005 | 6651 | **1.41, delayed** |
| `x_noforce` | `problem/rt_rad_force=false` | **ALIVE**, no collapse at all | > 11366 | **> 2.42** |
| `x_base` (1 rank) | -- | dies, (5,5,6,36) **r/R 0.791** | 5541 | 1.18 |
| `e2_nouform` | `hydro/rad_x1_uform=false` | **BITWISE `x_base`** -- a true null | 5541 | 1.18 |
| `s_nofofc` | `hydro/fofc=false` | dies, six cells at `i` = 92 | 5556 | 1.18 |
| `s_dc` | `hydro/reconstruct=dc` | dies, (2,6,4,97) r/R 1.011 | 5597 | 1.19 |
| `r64` | `nx1 = 64` (same stretch map) | dies, `i` 26 (r/R 0.804) then 64/66 | 5890 | **1.25** |
| `e_noimplx1` | `hydro/rad_implicit_x1=false` | REFUSED at startup, `conduction.cpp:254` | -- | -- |
| `f_noang` | `hydro/rad_implicit_ang=false` | REFUSED at startup, `conduction.cpp:488` | -- | -- |
| `s_novfill` | `mesh/cs_vertex_fill=false` | REFUSED, `parameter_input.cpp:390` | -- | -- |
| `s_wenoz` | `hydro/reconstruct=wenoz` | REFUSED (FOFC ghost count), `hydro.cpp:305` | -- | -- |
| `a8` | `nx2 = nx3 = 8` | no collapse by t = 5586 when its wall ran out | -- | inconclusive |

### 3.0 THE RANK COUNT PICKS THE RADIUS

`A1` and `x_base` are THE SAME CONFIGURATION, run on 2 MPI ranks and on 1.  `A1` dies at
t = 5602 s in cell (1,6,3,98), r/R = 1.013; `x_base` dies at t = 5541 s in cell
(5,5,6,36), r/R = **0.791**, the 3-D arms' own radius.  The only difference between them is
the order of the shell-mean MPI reductions, i.e. a round-off-level perturbation.

So the two radii of section 2 are ONE event, and which cell and which radius it takes is
decided at the level of 1e-16.  That is the signature of a marginal numerical instability
seeded by round-off, and it is incompatible with a physical front, which would have a
radius of its own.  `e2_nouform` (`rad_x1_uform = false`) came out BITWISE identical to
`x_base` -- same cycle, same time to seven digits, same cell -- which also confirms the
input file's claim that the radial conduction operator is inert under the whole-column
two-stream blend (`rad_tau_lo = 1e5`), including for the MLT closure's deficit.

`r144` did not get a slot.  The four refusals are configuration guards, not results; the
`nx1 = 64` arm reuses the SAME stretch polynomial (it is a map on the normalised index, so
the two grids are geometrically similar and nested, which is what a resolution test wants),
and its 1.25 turnover against 1.19 says **the death time barely moves with the radial
grid** -- a physical front would move with the grid or converge, and a marginal explicit
operator would not.

### 3.1 THE ONE NULL: `problem/rt_rad_force`

`x_noforce` ran to **t = 11366 s = 2.42 turnovers with no `dt COLLAPSE` print at all**,
against 1.19 for the baseline (it stopped on its own 3-minute wall, not on anything
physical).  It is the ONLY switch in the table that removes the event rather than moving
it by a tenth of a turnover.  This reproduces `tests_1d`'s section-D(ii) result
(`d_noforce` was one of two arms that survived the cycle-5 death of the un-floored 1-D
column) in the FLOORED, walled, `mlt_alpha = 1.5` configuration and at 20x the duration.

### 3.2 THE MLT CLOSURE IS INNOCENT, AND THE PER-COLUMN REMEDY IS WORSE

`tests_r6` section 5 item 1 asked whether the shell-mean closure's blindness to a single
column is the mechanism, and proposed `mlt_mean = false` as the remedy.  Both halves are
answered NO:

* with the closure switched off entirely (`mlt_alpha = 0`) the star still dies, 0.03
  turnover LATER, at a seam cell with `v1 = 1.000000e8 = vceil` and T at the table edge;
* with the per-column closure (`mlt_mean = false`) it dies 0.17 turnover EARLIER, which is
  the `mlt_mean` declaration's own documented failure mode (a per-column flux responds to
  the cell's own gradient and deposits what it carries into a few cells).

## 4. THE OPACITY: THE CELL IS ON THE FALLING SIDE OF THE Fe BUMP, SO THE FEEDBACK IS
## STABILISING

`tests_r6` section 5 item 1's third candidate was a physical thermal instability of the
Fe bump (heat -> more opacity -> absorb more -> heat more).  The table says the opposite.
The OPLIB He X = 0, Z = 0.02 Rosseland table at the 0.79 R cell's own state (T = 1.64e5 K,
rho = 2.34e-9, from `A1/column_A1.txt`, r = 1.8766e11):

| T [K] | 1.4e5 | 1.6e5 | **1.64e5** | 1.8e5 | 2.0e5 | 2.5e5 | 3.0e5 | 1e6 | 3.16e6 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| kappa at rho = 2.34e-9 | 0.713 | 0.802 | **0.796** | 0.718 | 0.562 | 0.313 | 0.231 | 0.204 | 0.199 |
| d ln kappa / d ln T | -- | +0.70 | **-0.94** | -1.75 | -2.53 | -2.19 | -0.62 | -0.01 | -- |

At fixed GAS PRESSURE (rho T = const), which is the relevant path for a perturbation that
stays in pressure balance:

| T [K] | 1.4e5 | **1.64e5** | 2.0e5 | 3.0e5 | 1e6 |
| --- | --- | --- | --- | --- | --- |
| rho | 2.74e-9 | **2.34e-9** | 1.92e-9 | 1.28e-9 | 3.84e-10 |
| kappa | 0.729 | **0.796** | 0.537 | 0.223 | 0.204 |

**The cell sits essentially ON the bump maximum** (0.802 at T = 1.6e5 against its own
0.796) with `d ln kappa / d ln T = -0.94`, and heating it by 20 % drops kappa by a third,
by a factor 1.8 at 2e5 K and by a factor 3.6 at 3e5 K.  A cell that heats LOSES opacity,
absorbs LESS and leaks more: the feedback is stabilising by a factor of four over the
runaway's own temperature range.  **The kappa-mechanism hypothesis is refuted for this
layer**, and with it the reading that the FeCZ is thermally unstable to a single-cell
perturbation that the closure cannot suppress.

## 5. THE VERDICT

**A NUMERICAL defect, and the operator is the RADIATIVE FORCE (`problem/rt_rad_force`).**
The evidence, in the order that matters:

1. It reproduces in a spherically symmetric, seedless 1-D column at the same time (1.19
   against 1.04-1.08 turnover), so it is not convection, not turbulence and not the seed.
2. No shell-mean quantity moves at all, in 1-D or in 3-D: T at 0.79 R is flat to 0.4 %
   from t = 0 to the collapse.  There is no front to be physical about.
3. It is a single cell, and that cell is in the panel's outermost angular row in 11 of 11
   3-D collapses (P = 1.4e-7 by chance).  A physical instability of a spherical layer has
   no reason to know where the cube's seams are.
4. Its terminal state is `|v1|` EXACTLY on `vceil` and T at or just below the EOS table's
   own upper edge, `eos_logt_max = 6.5` -> 3.1623e6 K -> 2.6293e14 in code units, which is
   the number printed in `A1`, `b_mlt0` and five of the six `s_nofofc` collapses.  The cell
   runs off the top of the table; it does not settle anywhere physical.
5. The opacity feedback at that radius is stabilising by a factor 4 (section 4).
6. `rt_rad_force = false` is the ONLY switch that removes it (> 2.42 turnovers clean);
   `mlt_alpha`, `mlt_mean`, `rt_strang`, `fofc`, `reconstruct` and the radial grid all move
   the death by 0.03-0.22 turnover and none of them removes it.

### 5.1 THE PROPOSED REMEDY (designed, NOT yet tested)

`rt_rad_force` applies `(1-w) rho kappa F/c + P_rad grad w` with `w = w(rho)` -- the EOS
radiation taper, `eos_rad_rho_hi = 1.93e-9` to `eos_rad_rho_lo = 4.97e-10`, read off the
column at tau = 3.0 and tau = 0.30.  **`w` is a function of DENSITY ALONE.**  The taper was
sized so that it "sits ENTIRELY ABOVE the FeCZ and straddles the photosphere" (the input
file's own comment), and that is true of the INITIAL column -- but it is not a property the
run preserves.  The ambient density at 0.79 R is 2.34e-9, a factor 4.7 above
`eos_rad_rho_hi`; `Gnr`'s collapsing cell there had `rho = 5.098e-10`, i.e. it had been
evacuated by 4.6x into the taper window, 2 % above `eos_rad_rho_lo`, where `w` is ~0.01 and
the cell receives **the full optically-thin surface radiative force at kappa F/(c g) = 1.23
in the middle of the convection zone**, plus the `(1-w) a T^4` the taper hands back as
internal energy.  Any cell that a seam-row flux error empties by half a decade is given the
surface treatment, and the force then empties it further: that is the loop, and it closes
on one cell because `w` is local.

So the two things to try, cheapest first:

1. **Gate `w` on optical depth or radius as well as density**, so that an interior cell
   cannot be handed the surface treatment because its own density dropped.  The natural
   form is `w = max(w(rho), w_tau(tau))` with `w_tau = 1` below the existing 20/300 or
   tau ~ 3 handover -- the taper is meant to mark WHERE the two-stream takes the radiation
   over from the gas, and that is an optical-depth statement, not a density one.
2. **Limit the force to the local g** (a per-cell cap at, say, 2 GM/r^2), which is a floor
   on the momentum error rather than a fix, but is one line and would tell immediately
   whether the acceleration is the whole story.

Both are guesses until measured.  The measurement that decides them is the arm above,
`x_noforce`, run to 5 turnovers: if the star settles with the force off, the force is the
whole story; if it dies of something else at 3 turnovers, this diagnosis covers only the
first killer.

### 5.2 WHAT WAS NOT DONE

The per-cycle energy budget of the collapsing cell (two-stream deposit, conduction, MLT
divergence, PdV) was NOT instrumented.  It would have needed a targeted debug print, a
rebuild and a rerun; the two null arms (`x_noforce` alive past 2.4 turnovers, everything
else dead within 0.22 turnover of the baseline) answer the same question -- which operator
-- without it, and section 4 removes the one candidate a budget was needed to rule out.
`a8` (`nx2 = nx3 = 8`) and `r144` never got a clean slot, so the ANGULAR resolution
dependence of the seam concentration is untested.  Both belong in the next round.

---

# PART B: THE MLT CLOSURE'S RELAXED PROFILE NOW GOES INTO THE RESTART FILE

`tests_r6` section 3 measured a restart kick: the radial kinetic energy jumped 2.5x (arm
G) / 4.6x (arm H) in the one history interval containing a restart, while the horizontal
kinetic energy was perfectly continuous, because `fmlt1d` -- the relaxed 1-D subgrid flux
profile -- was pgen state and not restart state, so the first source call of every process
re-seeded it to the instantaneous target (`red_giant.cpp:3262,3409`).  That is now fixed.

## B.1 THE MECHANISM, AND WHY THE FILE FORMAT DID NOT HAVE TO CHANGE

There was no hook for pgen state in the restart file at all: `restart.cpp` STEP 3 stores
the Z4c output time, the puncture trackers and the turbulence RNG and nothing else, and
`rot_potential` and `bcc0` -- the two previous restart bugs of this pgen family -- were
both fixed by REBUILDING the quantity, which a relaxed profile cannot be.  So this adds
one, as a **marked, optional block**:

| file | change |
| --- | --- |
| `src/pgen/pgen.hpp` | `PgenRestartStateFnPtr`, the 8-byte marker `kPgenRstMagic` = `"PGENST01"`, and two members: `pgen_rst_write_func` (enrolled by the pgen) and `pgen_rststate` (what a restart read back) |
| `src/outputs/restart.cpp` | STEP 3 writes marker + length + payload, and `step3size` accounts for it -- **only when the pgen enrolled a writer and returned a non-empty block** |
| `src/pgen/pgen.cpp` | the restart constructor peeks eight bytes where the variable data size used to be.  Marker -> read the length, the payload and then the real data size.  Anything else -> those eight bytes ARE the data size, and nothing has been lost |
| `src/pgen/red_giant.cpp` | `RedGiantMltRestartState()` / `RedGiantMltRestartRestore()`: `int32 nface, int32 seeded, nface Reals`, enrolled and consumed inside the existing `if (mlt_alpha_ > 0.0)` allocation block |

Why the peek is safe: the marker is exactly `sizeof(IOWrapperSizeT)` long (a
`static_assert` enforces it) and its eight bytes as an integer are 3.5e18, while the value
that occupies that position in the old layout is the per-MeshBlock variable data size,
O(1e6) on this grid and O(1e9) on any grid that fits in memory.  **So no switch is needed
in the input file and no version number in the file**: every restart file this code has
ever written takes the fall-back branch, and a `mlt_alpha = 0` run -- which is every box
run, every `red_giant_cs` run and every other problem generator -- writes no block, so its
file is byte-for-byte what it was.  Restarting a `mlt_alpha = 0` file with
`mlt_alpha = 1.5` prints a warning and re-seeds (the old behaviour); a stored profile of
the wrong length is a fatal error rather than a silent misread, because a face-by-face
profile has no meaning on a different radial grid.

`cpplint` over the four changed files: **34 errors before, 34 after** -- every one
pre-existing.

## B.2 THE VERIFICATION  (`rchk.sh`, `cmphst.py`, `firstline.py`, `rchk_saved.log`)

The chain is 1-D, `mlt_alpha = 1.5`, 0.5 turnover + restart + 0.5 against 1.0 turnover in
one piece, run with the OLD binary (`tests_r4/athena_v2`) and with the NEW one
(`tests_r7/athena_v3`), **with `problem/rt_rad_force = false`**: that is the one
configuration in which `tests_1d` section D(iv) established that a chained restart is
reproducible to 1 ulp, so any larger difference is the closure re-seed and nothing else.
The first resumed history line against the continuous run's line at the same time:

| column | OLD (no closure state) | NEW (closure state stored) |
| --- | --- | --- |
| `dt` | 0 (bitwise) | 0 (bitwise) |
| `mass` | **1.48e-08** | **0 (bitwise)** |
| `1-mom` (radial) | **8.82e-03** | **1.22e-16** (1 ulp) |
| `tot-E` | **2.29e-05** | **0 (bitwise)** |
| `1-KE` (radial) | **1.65e-01** | **9.82e-16** (4 ulp) |
| `2-mom` / `3-mom` | 1.07e-01 / 4.59e-03 | 8.32e-02 / 5.48e-03 |
| `2-KE` / `3-KE` | 1.39e-03 / 6.32e-04 | 1.38e-03 / 2.91e-04 |

and the kick keeps growing in the old binary (`1-KE` 16.5 % -> 25.4 % -> 31.4 % over the
next two history lines) while in the new one it stays at 3e-15.  The number of history
times the two legs still agree on rose from **2 to 35**.

The restart prints what it did:

    ### red_giant: MLT closure state read from the restart file (103 faces, seeded = 1)

103 = `nx1 + 2*nghost + 1` = the face count the closure allocates.

**Not bitwise, and here is the honest reason.**  Two residuals survive.

1. `1-mom` and `1-KE` differ by 1-4 ulp on the first resumed line.  That is the
   reduction-order difference `tests_1d` D(iv) already documented for this code (the
   restarted leg's shell sums and history sums are formed over the same data in a
   different order), and `dt`, `mass` and `tot-E` are bitwise, so it is not a state
   difference.
2. The ANGULAR columns `2-mom`, `3-mom`, `2-KE`, `3-KE` differ by 1e-3 to 8e-2 in BOTH
   binaries, old and new alike.  They are cancellation sums: `2-mom` = 4e15 against
   `1-mom` = 3.7e31, i.e. 1e-16 of it, and `2-KE` = 1e10 against `1-KE` = 1.2e37, i.e.
   1e-27 of it.  A relative difference of 8e-2 on a quantity that is 1e-16 of its own
   radial partner is round-off in the sum, not motion; it is present in the unmodified
   binary at the same size and is therefore not something this change introduced.

So the claim this round makes is the one that matters for production: **the restart kick is
gone** -- a factor 1.65e-1 in the radial kinetic energy became 9.8e-16 -- and the residual
is the pre-existing 1-ulp reduction-order noise.

## B.2a THE FORMAT GATE  (`fmt.sh`, `fmt2.sh`, `fmt.log`, `fmt2.log`)

Five checks that the file format is untouched where it must be:

| check | result |
| --- | --- |
| `mlt_alpha = 0`, old binary vs new: file SIZE | **2968901 bytes both** -- no block written |
| `mlt_alpha = 0`, old vs new: bytes | 21 bytes differ in 31017..31052 -- **and the OLD binary run TWICE differs in the same 21 bytes in the same range** (`fmt2.sh`), so it is not this change |
| new binary reads the OLD binary's `mlt_alpha = 0` file | runs, no message, `nlim = 8` reached |
| new binary reads an OLD file with `mlt_alpha = 1.5` | prints `this restart file carries no MLT closure state; the relaxed profile is re-seeded` and continues -- the old behaviour, announced |
| OLD binary reads a NEW `mlt_alpha = 1.5` file | **refused loudly**: `CC data size read from restart file not equal to size of Hydro, MHD, Rad, and/or Z4c arrays, restart file is broken` (`pgen.cpp:313`).  Backward compatibility is one-way by construction; it fails on the size check rather than reading garbage |

**A pre-existing defect found by the second row, NOT fixed here.**  `restart.cpp` STEP 1
writes `mesh_indcs` and `mb_indcs` as raw `RegionIndcs` structs, and that struct's
`cnx1..cke` coarse-cell fields (`mesh.hpp:41-42`) are only ever assigned when the mesh is
multilevel.  On a uniform grid they are never initialised, so **every restart file this
code writes contains ~21 bytes of indeterminate memory and two runs of the same binary
produce different files.**  It is harmless to a restart (the reader takes the coarse
indices from the input, not the file) but it makes byte-for-byte restart regression tests
impossible and it writes uninitialised memory to disk.  One-line fix: zero-initialise
`RegionIndcs`.  Left alone in this round because it is unrelated and touches the mesh
header.

## B.3 WHAT THIS UNBLOCKS AND WHAT IT DOES NOT

`tests_r6` section 5 item 2 is closed: a `mlt_alpha > 0` run can now be chained across
24-hour slots.  It does NOT make the star live longer -- arms G and H of `tests_r6` died at
0.6 turnover because of the kick, but `Gnr` died at 1.08 without one, and part A shows why.

## Files

* `README.md` -- this file.
* `one1d.sh` (arm `A1`), `bis1d.sh` (`b_mlt0` .. `g_percol`), `seam1d.sh`
  (`s_novfill` .. `r144`), `chk.sh` (`x_noforce`/`x_nostrang`/`x_base` + the first restart
  chain), `rchk.sh` (the four-chain old-vs-new restart comparison).
* `prof.py` -- the `rt_profile.bin` reader (shell means, hottest cell per record);
  `cmphst.py`, `firstline.py` -- the history comparators.
* `A1/`, `b_mlt0/`, ..., `rc_*/` -- `he4.hydro.hst`, the column and `mltfaces` dumps.
* `A1.log`, `bis.log`, `seam.log`, `chk.log`, `rchk_saved.log` -- the job logs.
  `rchk.log` was overwritten by a duplicate submission; `rchk_saved.log` is the copy that
  carries the `oldNF`/`newNF` comparison quoted above.
* `athena_v3` -- this round's binary (`build_gpu_rg`, MI300A, `PROBLEM=red_giant`).  NOT
  committed (46 MB); rebuild it from `52d16d6b`.  Part A's arms all use
  `tests_r4/athena_v2` = `5b71b144`.
