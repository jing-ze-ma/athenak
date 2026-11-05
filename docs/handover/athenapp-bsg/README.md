# Athena++ runs of Linhao's BSG (Raven, 10-15 .. 10-16): YES, we ran it

Two worker campaigns on Raven, both 1-D (spherical, x1 only). Both NOTEs are copied here in full.
- `audit_1015/NOTE.md`, Raven `/raven/ptmp/jinma/athpp_audit_1015`: "Athena++ radiation audit re Linhao BSG collapse".
- `fix_1016/NOTE.md`, Raven `/raven/ptmp/jinma/athpp_fix_1016`: "converge Athena++'s implicit radiation solve (Linhao BSG, 1-D)".

## Build
- Source: Athena++ PrincetonUniversity/athena master **e47fabb7** (2026-10-14), plus Linhao's `bsg_star.cpp` (output part patched).
- Configure: pgen bsg_star, spherical_polar, adiabatic, implicit radiation ON (see `audit_1015/bld/conf_bsg.log`). Compiled with g++ 7.5 -O3 on the CPU.
- Patches:
  - `audit_1015/athena_led_vs_e47fabb7.patch`: energy ledger (led.hpp) and the thickslab test pgen;
  - `audit_1015/bsg_star_as_run.cpp`: the pgen as run;
  - `fix_1016/conv-fix.patch`: the radiation/sweep=2 block-tridiagonal line solve. It is opt-in, and bitwise Athena++ when off.

## Setup
- Units: rho 1e-9, T 1e5, l Rsun; mu 0.62; GM 20 Msun.
- Grid: 20-80 Rsun, nx1 384 uniform.
- Radiation: nmu 4, cfl 0.3 (dt 212 s); his MESA IC + tables.
- Arms:
  - taucell 1/4/10/20/40/400;
  - angle_flag 0/1;
  - floor on/off (nf), Doppler boost on/off (nb);
  - closed wall (Bd1) vs his inner bath BC.
- Decks are in `audit_1015/bsg1d/athinput.*`, `base.in` + `mkinputs.py`, and `fix_1016/runs/*.athinput`.

## Results (details and numbers in the NOTEs, the figures, and `*_summary.txt`)
- **Every Athena++ arm collapses at the inner wall.**
  - Wall cell 1: T 1.05e6 -> 2-3e5 K, rho x3-4, within 1-3 d.
  - The wall flux is -100..-400 L, INTO the bath.
  - The total energy drains: dE -3.4e46 erg in 5 d (af0), -7.1e46 in 10 d (af1).
  - The envelope deflates: r(rho=1e-10) goes 54 -> 49 Rsun (10 d) -> 45 Rsun (17 d), then blows out at ~18 d.
  - af0 arms die at 5 d; af1 at ~24 d. taucell 40/400 diverge at cycle 1-2.
- **Two causes found:**
  - the HLL tau-factor numerical diffusion carries 3/(8 taucell) of L with no matching radiation force;
  - the default Jacobi iteration is unconverged, and its billing does not telescope, which drains energy (-13..-27e45 erg in 5-10 d).
- **The fix removes only the iteration loss.** With sweep=2 (block-tridiagonal line solve) the non-telescoping part is 0 to round-off, the run costs 0.2-0.35x the default, and the closed-wall case matches the converged Jacobi reference.
  Even so, the CONVERGED discretisation still develops the base wall instability: Bd1 dies at 2.7 d both with the fix and with converged Jacobi.
  The tc4 run with his BC survives 30 d, but the bath drains -31e45 erg by 5 d and the envelope blows out at ~18 d.
- **Head-to-head with AthenaK** (h2h; the deck 384x8x8 stands in for 1-D): the AthenaK arms K0-K3 (his IC/BC, FIN5.1 keys) all run 30 d with no FATAL.
  - Wall T stays at 1.05e6 K (+0.8 %).
  - L_in is 2.5-2.8 L INTO the star.
  - dE_tot is +1.6e45 erg.
  - => The shrinking/collapse appears in Athena++ and not in AthenaK. Its source is Athena++'s radiation discretisation at the base, not the physics.
- Figures:
  - `audit_1015/bsg_collapse_audit.html` (page);
  - `audit_1015/*_user.png`;
  - `audit_1015/h2h/h2h_user.png`, `h2h_timeseries.png`, `h2h_profiles.png`, `bcswap_*`;
  - `fix_1016/fix_user.png`.
