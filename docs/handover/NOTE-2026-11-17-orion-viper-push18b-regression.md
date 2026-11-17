# NOTE 2026-11-17 Orion -> viper (USER): push 18b CPU regression of 42f94c81b vs push 18 b53e35565

This answers TASK-2026-11-17-viper-orion-push18b-regression.md.

**A: bitwise except the reflecting-x1 FV mirror class, plus round-off.** No CLEAN WARNING / DT_RETRY / FATAL in any A run.
The movers are cs feed1, the RT feed decks (feed, feedx1, feedmr, csCnr), pulse0, flz_q/flz_Q and dif/lat. Every one of
them has a reflecting x1 face with FV rays, either set as fv1 or taken as the default for cell + vet_gd.
- The bisect puts cs feed1, pulse0, flz_q and lat in 9fea91e22..c33ff8dd2. That range is 12ea9dcbb (= your 37686d20a) plus
  81be1beb0 (the tie rule), 8b716f18e (the zero-step restart) and c33ff8dd2 (ParameterDump).
- The RT feed decks and csCnr are confirmed only at the c33ff8dd2 end. flz_Q is not bisected.
- dif/lat move at 1e-13 and rhse at 7e-18.
- Two movers are not on your list:
  - **pulse0** (O(1) in F, dE stays at round-off).
  - **flz_q**, a chaotic deck. It now reaches t 86.0 instead of 32.9 in 40 steps, with 6966 instead of 7909 iterations, and
    CLEAN cycle 0 already differs.
- Conservation is unchanged on all of them.
**Restart: bitwise.** That holds in all 9 A restarts, the 3 extrap restarts, and the forced retry with views = clean and
with views = all + poison.
**Extrap arm: unchanged from push 18.** rw F/X are byte-identical to push 18. rw with `clean_rays = fv1` and no vet_gd keys
still segfaults (pre-existing).
**clean_stage_floor: FINDING.** In a closed box whose equilibrium lies below tfloor, the default run (dt_retry on) does not
conserve energy by converging. It stalls: dt is cut down to 1e-16 and t freezes (Zeno), so E stays at round-off only because
the box stops evolving. With dt_retry off, the floor energy (+11..34 %) is still there, and the run turns into consecutive
backward-Euler steps: all 40 of 40 steps are BE, with 6 fallbacks and 33 time2_auto_be steps.
**Cell layout and fref_wacc: not booked.** The cell write-back pays h rho a_ref.v at the velocity after the kick, in stage
and BE steps alike. That is the booking our fcfc0929e / dccd3e031 replaced in the staggered layout.
**Timing (Raven, dedicated node):** cs_closed +3.4 % (all 3 rounds), clean stage +3.4 %. flz_q is not comparable, because
its trajectory changed: -14 % per cycle, -2 % per iteration.

## Binaries (Orion, git archive trees, gcc/13 + openmpi/4.1, Release, MPI+OpenMP, -ffp-contract=off, NMAX 4096, built_in_pgens, kokkos d8e9af031)
| head | md5 |
|---|---|
| 42f94c81b (push-18b candidate) | `85e0c3273c5f8f03840009b7592efc2e` |
| b53e35565 (push 18, rebuilt in this dir) | `747da680df7eaed54e5cc8bf28d45108` |
| bisect 9fea91e22 / c33ff8dd2 / 90428c592 | `a1653869` / `71b779c6` / `a48f0a82` |

- The b53e35565 md5 differs from push 18's 851a7eb5 only because the source path is embedded in the binary. Its A-arm data
  are identical to push 18's A18 (checked on the decks the bisect reran).
- Runs used hold 231122 (srun --overlap, explicit -c, OMPI_MCA_hwloc_base_binding_policy=none).

## A: 42f94c81b (A19) vs b53e35565 (A18), no keys
Every run has rawdiff > 0 because of the ParameterDump widths (c33ff8dd2), so the bin headers differ. The table compares
data blocks only.

| run | bindata | maxrel | first differing file | bisect |
|---|---|---|---|---|
| cc ladder x5, thin_p/r/f, s88, s81, cs uni/closed/thin1/thin1mb/thin1mbmr, rw x4 | same | 0 | none | |
| RT hcyl, sp3mb, sp3mbmr, thin1, rot1s (incl. restarts) | same | 0 | none | |
| WB fv1 / fv / lld / cart | hst same | vcs_max 1.21364e-09 / 2.96263e-09 / 4.75321e-09 / 7.38032e-09 | | |
| rhse | hst | 7.1e-18 | rh.user.hst | (your rhse_sp2 round-off) |
| dif / lat | 1 of 2 | 5.7e-14 / 1.1e-13 | sp.m1.00001.bin (cycle 241) | lat: 9fea91e22 = A18, c33ff8dd2 = A19 (dif: c33ff8dd2 = A19) |
| cs feed1 | 2 of 4 | 9.6e-01 | cs.hydro_u.00001.bin (cycle 10), mom1 | 9fea91e22 = A18, c33ff8dd2 = A19: **12ea9dcbb** (37686d20a) |
| RT feed / feedx1 / feedmr | 8 of 12 | 2.3e-02 | f/bin/cs.hydro_u.00001.bin (cycle 16), mom1 | c33ff8dd2 = A19 (the 9fea91e22 end was stopped; its CLEAN cycle 0 equals both) |
| RT csCnr (no clean_rays key: fv1 default) | 8 of 12 | 2.5e-02 | f/bin/cs.hydro_u.00001.bin (cycle 16), mom1 | c33ff8dd2 = A19 (9fea91e22 end stopped) |
| **pulse0** (reflect x1, fv1 default) | 19 of 20 | 2.5e+00 at the last dump, 3.0e-02 at the first | sp.m1.00001.bin (t 0.215), F2 | 9fea91e22 = A18, c33ff8dd2 = A19: **12ea9dcbb** |
| **flz_q** (fv1, reflect x1, chaotic) | 2 of 2 | 1.1e+00 | sp.hydro_u.00001.bin (cycle 40), mom1; CLEAN cycle 0 already differs | CLEAN cycle 0: 9fea91e22 = A18, c33ff8dd2 = A19 |
| flz_Q | 2 of 2 | 9.1e-05 | sp.hydro_u.00001.bin (cycle 40), mom2 | not bisected (stopped; node load 800+) |

## B table (42f94c81b vs b53e35565)
| deck | metric b53e35565 | metric 42f94c81b |
|---|---|---|
| cc hcyl cbg cb shd nltek1 | 1.2061e-02 1.1003 0.45999 1.1805e-02 9.0716e-05 | same |
| thin_p / r / f, dif, lat | L1(E) 1.289e-01 / 1.179e-01 / 2.942e-02, 1.536e-02, 2.915e-01 | same |
| pulse0 | dE +8.882e-16, it 58 | dE +4.441e-16, it 68 |
| s88 / s81 | 2.578e-07 | same |
| rhse | vcs_max 1.85058e-01, it 806 | same |
| cs_uni / cs_closed | 0 / sum E V rel 0.00e+00 | same |
| cs_thin1 / thin1mb | L1(E) 3.891e-02 | same |
| cs_feed1 | dE rel 0.00e+00, it 257, E_rad 7.0052e+01 | dE rel -1.11e-16, it 255, E_rad 7.0040e+01 |
| flz_q | t 32.88 in 40 steps, dLz 5.87e-16, it 7909 | t 86.04, dLz 7.83e-16, it 6966 |
| flz_Q | t 80.39, dLz 1.25e-15, it 711 | t 80.39, dLz 1.10e-15, it 711 |

- 0 CLEAN WARNING, 0 DT_RETRY and 0 FATAL in every A run of both builds.

## Extrap arm, restarts, forced retry (42f94c81b)
- **Extrap:** sp3, thin1 and hcyl restarts are bitwise, and so is the hcyl fv1 base. cs thin1 X: L1(E) 3.891e-02, 98
  one-pass / 0 swept. cs feed1 X: dE -1.11e-16, 5 one-pass / 3 swept. flz_q X: 31 one-pass / 20 swept, it 6901.
  0 W / DT / F throughout.
- **rw fv1 + vet_gd keys (F19) and with extrap (X19):** the tabs are byte-identical to push 18's F18 / X18, at both N and
  both signs. So extrap's effect on rw is unchanged: F1 3.6e-1 at N32 and 3.4e-1 at N64.
- **Forced retry** (`dt_retry_force_cycle = 4`, fac 1, 16 cycles; 1 `reason=forced` each): sp3, thin1 and hcyl are bitwise
  to the straight extrap runs with `dt_retry_views = clean`, and also with `dt_retry_views = all` + `dt_retry_poison = true`.
  Their own restarts are bitwise. One-pass count 14 = 14 (push 18 had 15 vs 14).
- **A restarts:** 9/9 give `hydro_u same m1 same` (hcyl, sp3mb, sp3mbmr, feed, feedx1, feedmr, thin1, csCnr, rot1s).

## clean_stage_floor (default true): closed boxes with tfloor
- **cs feed-box (cubed sphere), gas T 1 -> T_eq ~0.77, tfloor 0.9:** the hydro's tfloor never touches the cons on this
  mesh: T stays at 0.62 < 0.9 in the output, and a kappa = 0 box at T 0.5 stays at 0.5.
  - So b53e35565 conserves E to 0 here.
  - 42f94c81b **FATALs at cycle 0** (`dt_retry: cycle 0 failed 9 times`, reason=eint<=0, the first-step BE write-back). This
    happens with clean_stage_floor on and off, so the cause is the write-back check 0e8a3a00b + f198987d0, not the stage
    floor.
  - With rho_lat_amp 0 it does not FATAL, but it stalls: 32-52 eint<=0 retries in 40 steps, last dt 1e-11..1e-17, and t
    reaches 2.6e-2 / 2.5e-3 / 2.5e-4 (kappa 1 / 10 / 100) instead of 1.15.
- **Cartesian periodic box** (rad_m1_beam equil, 32x4, rho 1, eint 1.5 (T 1), E 1e-3, arad 1, tfloor 0.9, c 10, 40 cycles):
  here the hydro floor acts, and T_eq ~0.77 < tfloor. Columns: `p18` = b53e35565; `on` / `off` = 42f94c81b with
  clean_stage_floor true / false; `nr` = dt_retry false.

| kappa | arm | t after 40 cycles | E(gas+rad) rel | DT_RETRY (eint<=0), last dt | stage / BE steps, fallbacks, auto_be |
|---|---|---|---|---|---|
| 1 | p18 | 0.616 | +3.34e-01 | 0 | 39 / 1, 0, 0 |
| 1 | **on (default)** | **0.0239** | +8.2e-15 | **50, dt 1.1e-16** | 39 / 1, 0, 0 |
| 1 | off | 0.523 | +3.32e-01 | 1 | 39 / 1, 0, 0 |
| 1 | on, nr | 0.616 | +3.32e-01 | 0 | **0 / 40, 6, 33** |
| 1 | off, nr | 0.616 | +3.34e-01 | 0 | 37 / 3, 1, 1 |
| 0.1 | p18 / on / off | 0.616 / **0.210** / 0.508 | +1.11e-01 / -2.0e-14 / +8.6e-02 | 0 / **40, dt 2.9e-14** / 1 | 39 / 1 |
| 0.1 | on, nr / off, nr | 0.616 | +1.09e-01 / +1.10e-01 | 0 | 12 / 28, 5, 22 / 37 / 3, 1, 1 |
| 10 | p18 / on / off | 0.616 / **0.0025** / 0.385 | +3.37e-01 / +8.9e-16 / +3.37e-01 | 0 / **50, dt 2.8e-17** / 3 | 39 / 1 |
| 10 | on, nr / off, nr | 0.616 | +3.37e-01 | 0 | 0 / 40, 6, 33 / 39 / 1 |

- **Reading:**
  - With dt_retry on, every step that would take the gas through tfloor is rejected and dt is halved. The gas creeps up to
    the floor and dt collapses: 1e-16 after 40 cycles, at t 0.024 instead of 0.62.
  - E(gas + rad) is at round-off only because the box stops evolving. Please check t (or the dt column of the hst) in your
    +9.7 % -> round-off box: "15 retries" may be the same stall.
  - When the equilibrium is below tfloor, no dt satisfies the check. The floor's energy is physical for that setup, so
    rejecting cannot converge.
  - My suggestion: do not retry or redo when the cell's own exchange root (gas + E of the cell alone at h -> infinity) is
    below tfloor. Count that floor energy instead (a CLEAN line or a WARNING once). Keep the rejection for the
    undershoot-only case: root above the floor, iterate below it.
- **Decks where the stage floor does fire turn into BE runs.**
  - With dt_retry off, `on` gives 5-6 stage fallbacks in 40 steps. The cell layout's default time2_auto_be = 16 then makes
    22-33 of the next steps BE, and at kappa 1 and 10 all 40 of 40 steps are BE.
  - On our staggered R136 restart (r136.00009), BE solves in a row locked the face F of a semi-transparent moving layer into
    an odd-even oscillation: v ~2.5e8 cm/s, hydro dt 5.4 -> 2.1 s. That is why staggered now has time2_auto_be = 0
    (stag2 32b045052).
  - The cell layout keeps 16, so clean_stage_floor can meet that instability in a moving, semi-transparent layer. I could
    not reproduce it in these uniform boxes: there is no transport there, and no instability showed.
  - Also, the "off" arm with dt_retry on shows that once a cell sits at the floor (eint = ef, not > ef), neither check sees
    it again.
- **fref_wacc in the BE path (wb_arad, force_reference_work = split):**
  - **The cell layout does not book fref_wacc.** rad_m1_clean_cell.cpp CleanCellWriteBack (around line 2023) pays
    wr = h rho a_ref.(p + dm)/rho at the velocity after the kick, in stage and BE steps alike. No rad_m1_clean*.cpp file
    reads FrefWacc().
  - The problem generator (he_star_m1.cpp hs_grav_fws_x around line 3788) gives the gas bdt rho a_ref v at each hydro
    stage's start velocity. In a hesdirk2 BE step (FrefWaccOn: first step, a redo, time2_auto_be) it also accumulates
    fref_wacc, which nothing in the clean scheme consumes.
  - On the staggered layout that same mismatch left +0.13 L_top dt per BE step on r136.00009. fcfc0929e fixed it: in a BE
    step E pays h/dt fref_wacc.
  - The stage version left 1.6-3.1e-5 L_top dt per stage step. dccd3e031 fixed it: each stage pays the work of its hydro
    stage (gam0 = 0), and no slope carries a payment.
  - clean_stage_floor adds BE steps, so it adds steps that hit the BE half of this leak. Porting both fixes to the cell
    write-back (wr -> fra*fref_wacc where FrefWaccOn, and the slope rule of dccd3e031) would close it.
  - I have not measured this on a cell-layout wb_arad deck.

## Timing (Raven job 32466704, one dedicated node ravc3286, Xeon 8360Y, `general --exclusive`; 3 rounds a: P18, 18b; b: 18b, P18; c: P18, 18b)
| build | flz_q s/cycle (it) | cs_closed s/cycle | cs_closed clean stage steps (s) |
|---|---|---|---|
| b53e35565 | 2.665 (2.685 2.665 2.645), it 7909 | 2.201 (2.184 2.221 2.198) | 3.449 (3.436 3.466 3.446) |
| 42f94c81b | 2.297 (2.284 2.307 2.299), it 6966: **-14 %/cycle, -2 %/iteration** | 2.277 (2.259 2.292 2.278) **+3.4 %** | 3.567 (3.536 3.595 3.569) **+3.4 %** |

- cs_closed: the data are identical (A) and every 42f round is slower than every b53 round. So this +3.4 % is a real CPU
  cost, outside the ±1-2 % spread.
  - Candidates: ed3e404b3, the BiCGStab true residual (one more operator application per converged solve; closed has
    122 its).
  - The cln_old parallel_reduce of c90273fd5 (it now reduces even where clean_stage_floor cannot act).
  - Not separated.
- flz_q's trajectory changed (larger dt, fewer iterations), so its s/cycle is not a like-for-like timing.
- Binary md5s on Raven: ab076454 (b53e35565), e6db2ad8 (42f94c81b).

## Files (Orion: /orion/ptmp/jinma/clean_1115/push18b_reg)
- Source trees and builds:
  - src_p18 and src_p18b, with builds build_p18 and build_p18b.
  - The bisect trees src_9fea91e22, src_c33ff8dd2 and src_90428c592, with their build_* dirs.
  - The binaries athena_*.
- Harness and comparisons:
  - reg/ is the push-18 harness copy, without its arm outputs.
  - The arms are A18, A19, X19 / R19 / B19 / F19, A19_RTA* (views = all + poison), and the bisect arms K1/K2/K3 (the slow 9fea91e22 RT and flz runs were stopped by PID on the loaded node)
    (9fea91e22 / c33ff8dd2 / 90428c592).
  - The comparisons are cmp_A19_A18.txt and ana_A18_A19.txt.
- Stage-floor boxes: sfl/ (cubed sphere, fl_*.athinput) and sflc/ (Cartesian eq_k*, plus radwave a*_t* decks where only the
  hydro crosses the floor: all arms equal).
- Scripts: runA.sh, xarm19.sh, bis.sh, sfl/go.sh and tools/{cmpx.py, pd.py, rwcmp.py}.
- Raven timing: /raven/ptmp/jinma/stag2/p18btim (job18b.sbatch, node18b.log). Its runs are in
  /raven/ptmp/jinma/stag2/p16bis/reg (*P8B*).
- All .rst dumps were deleted (Orion and Raven).

-- Orion
