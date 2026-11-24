# NOTE 2026-11-24zm orion -> Caltech: dfloor_flux_limit fixed by construction (bsg c777); arms H1 / T1 / H2 / H3 when a node frees

Branch pushed: `caltech-dflfix-5b6159ad0` (rc4 c21bc2c33 + 1 commit). Test arms only, no production switch yet. Run them on a node that frees up, or beside a production run if it fits. Never stop production for them. The H200 arms are this fix's GPU gate (MI300A is still queued). Report per the table. If H1 / H2 / H3 are clean, I'll send the production step.


## Code
- Branch dflfix-rc4, tip 5b6159ad0 = rc4 (ia123-rc4 c21bc2c33) + one commit. COORD / MERGE123 push it as caltech-dflfix-5b6159ad0 (or carry it in the next rc).
- Build as usual (CUDA 12.9 + g++13), PROBLEM=he_star_vet.
- Restarting an rc3 rst, or using the old gate deck, needs `rad_vet/vet_gd_ls=0`: rc4's 26edb6abf FATALs when a product set and vet_gd_ls are both named.
- The fix has no key. It is part of `hydro/dfloor_flux_limit`, which every star deck names.
  - A cell whose outflow the limiter scaled by th < 1 gets its momentum scaled by the same th, except the L_z part (phi_hat, as in floor_from_ekin). E is held fixed, so the kinetic energy becomes heat.
  - This is applied after all of the stage's momentum sources.
  - The held floor gas then moves only as fast as the mass flux it really had. Gravity can no longer push it to vceil while its potential-energy release (carried by the throttled mass flux under etotgrav) never arrives. That is the c777 mechanism BSGINF found.
  - The end-of-run line `<hydro> dfloor_flux_limit: outflow-limited cell-stages (all attempts)=N, their x1 momentum x vol removed=P, kinetic energy x vol made heat=K` books what it took.
  - A run in which no cell is limited is bitwise the same as rc4.

## Orion result (CPU, reduced bsg 256 x 64 x 64, BSGINF's rc3 deck as ran, from its rst at t = 4.32e4)
dt at about the same t (the arms' 50-cycle lines):
| t | fix 5b6159ad0 | base rc3 af836f5a3 | dfl off (rc3) |
|---|---|---|---|
| 6.2e4 | 302 | 302 | 302 |
| 7.6e4 | 304 | 294 | 304 |
| 9.0e4 | 304 | 239 | 304 |
| 1.07e5 | 306 | 157 | 304 (stopped at 1.07e5) |
| 1.22e5 | 303 | 122 | |
| 1.37e5 | 302 | 64 | |
| 1.53e5 | 304 | died at 1.50e5 (dt 1e-3, closure inf) | |
| 1.68e5 | 306 | | |
- A FRESH fix arm from t = 0 (rc3 deck + vet_gd_ls=0, to t 1.7e5) prints the same 50-cycle dt lines as the restart arm, to every digit, and it is also clean: 0 floored / kedom / e_c<=0 cells, 0 CW. An rc4 c21bc2c33 base arm from the same rst reproduces rc3's failure (first floors at c215, r 4.011e12; dt 294 at c250).
- Fix arm, to tlim 1.7e5: floored cells 0 in every cycle (the base had its first at c215, r 4.011e12); e_c <= 0 0; kedom 0; EOS floor energy created 0; 0 CLEAN WARNING; 0 DT_RETRY.
- dfl line: 3e8 limited cell-stages; x1 momentum x vol removed -6e31 (inward, i.e. the infall it stopped); heat 2e37 erg over 1.25e5 s. For comparison, the base's floors paid 2e39 erg per cycle before it died.
- Against the dfl-off arm at matching t (6e4 - 1.05e5): L_top within 1e-4; KE_int within 1e-3 (the floor gas KE); M_tot -1.8e-8 (dfl off creates mass at the c2p floor; the fix conserves it). Against the base before its failure: M_tot equal to 1e-15, L_top within 1e-4.

## Arms (1 node, 4 H200; bsg FRESH from t = 0 with the rc3/rc4 bsg deck; nlim = 1000, ~11 min each)
- Common args: `time/nlim=1000 output2/dt=-1 output3/dt=-1 output5/dt=-1 output6/dt=-1 output7/dt=-1` (hst + log only), with `hydro/floor_audit=1` kept.
- Report for each arm:
  - the first CW cycle, or none;
  - dt at cycles 200 / 400 / 600 / 777 / 1000;
  - the eos_efloor column of bsg3d.log;
  - the end-of-run dfl line;
  - the c770-790 floor_audit lines.

| arm | binary | args (on top of the common ones) | expected |
|---|---|---|---|
| H1 | 5b6159ad0 | (none: dfl on, as the deck) | clean past c1000, dt steady ~300 like G1 (dfl off) |
| H2 | 5b6159ad0 | `hydro/fofc=true hydro/fofc_retest=true` (production form) | clean, dt as H1 (FOFC first-order cell-stages reported) |
| H3 | 5b6159ad0 | production restart: newest bsg rst (FOFC on) + `rad_vet/vet_gd_ls=0` if it is an rc3 rst, 1000 cycles | clean, no new floored cells at r ~ 4e12 |
| T1 | 5b6159ad0 vs the binary t9l62 runs now | t9l62 restart from its newest rst, 200 cycles, both binaries | the dfl line of both; if N = 0, hst bitwise (the fix does not act); else L_top / L_in and the floor_audit lines of both |

Order: H1, then T1 (it answers whether t9l62's WD0 floor energy is this effect), then H2, H3.

## t9l62 (WD0 floor energy): orion's reading
- The WD0 floored cells have rho 1e-8 .. 1e-6 at r ~ 2.44e12, and dfloor is 3.6e-19. The limiter acts only where a stage would take rho below dfloor, so those cells cannot be dfl-limited cells.
- WD0's energy is therefore not this mechanism. T1 checks it directly: dfl line N and hst bitwise.

-- orion (coordinator)
