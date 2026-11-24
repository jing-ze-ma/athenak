# NOTE 2026-11-24zk orion -> Caltech: floorcap_lz_spin (restores m35mms dt with exact L_z): build, test FS1 / FS1r / FS2, then switch m35mms

Branch: `caltech-ia123rc3-fcapspin-6f2973dc8` = rc3 af836f5a3 + one commit (rad_vet/floorcap_lz_spin, read only when named; key off bitwise to rc3). Build incrementally (rad_vet_fluid.cpp / rad_vet.cpp change), tiny gate with the seam keys (3.686e-11).

**Plan:**
1. At m35mms's next rst stop, run FS1 + FS1r below on its node (~2 min), then FS2 when M35ZAMS is at an rst.
2. If FS1 is clean (CW 0, FATAL 0, restart bitwise) and dt returns to ~10-11 (the pre-floorcap_lz value), continue **m35mms production on this binary with `rad_vet/floorcap_lz_spin=true`** (keep the pre-switch rst). Under the 24r watch rules. This keeps L_z exact (floorcap_lz) AND keeps floor gas under the cap, so it doesn't bind dt.
3. M35ZAMS: expect no dt gain (floor gas at the cap still binds there); report FS2 anyway.


Problem (Caltech M35DT_LZ, dtrep_dtreg.txt): with `rad_vet/floorcap_lz = true` the floor-gas cap keeps v_phi, so floor gas (rho = dfloor 1e-17) keeps |v_phi| 3.8 -> 6.7e7 > floorcap_v 3e7 and binds dt (dir = 3): m35mms dt 11.2 -> 5.9. With floorcap_lz = false, L_z is not exact (ledger: the cap removed 4.7e30, Lz_tot -3.6e30).

Fix (key `rad_vet/floorcap_lz_spin`, read only when named; key off bitwise): the cap scales the WHOLE velocity of floor cells (|v| <= floorcap_v, as with floorcap_lz off), and the L_z X it removes goes to the gas at or above floorcap_rho as one rigid rotation dOmega = X/I, I = sum rho V R^2 (fixed-order sum, layout-invariant in tile mode). E_rad of each receiving cell pays the rotation's kinetic energy. L_z and E exact with no ledger. The floorcap report's `dL_z` field is now the NET (removed - returned), so it should read ~0 (round-off of the atomic counter sum).

## Build (rc3 + 1 commit, CUDA 12.9 + g++ 13, the same recipe as athena_ia123rc3_af836f5a3)
- Source: af836f5a3 + `git am stag2/fcapspin/patch_rc3/0001-*.patch` = 6f2973dc8 (branch fcapspin-rc3, local). rc4 form: 96ba7f353 + patch_rc4 = bc0ccc884 (branch fcapspin).
- Key off = bitwise vs af836f5a3 (Orion CPU, M35ZAMS cut, hst + dumps). Tiny BSG gate pass 1 must be 3.686e-11.

## Arm FS1: m35mms (one node, 4 H200, scratch dir, no bin / rst outputs)
Same restart and args as M35DT_LZ, plus the spin keys:
```
-r /resnick/groups/carnegie_poc/jingze/prod_rc3_1124/m35mms/rst_keep/m35mms.00006.rst.pre-lz
time/cfl_number=0.3 hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4
rad_vet/floorcap_lz=true rad_vet/floorcap_lz_spin=true rad_vet/floorcap_mesh_vol=true
rad_vet/floorcap_v=3e7 rad_vet/floorcap_rho=1e-14 rad_vet/floorcap_report=10
hydro/dt_report=10 time/nlim=<rst cycle + 50> job/final_restart=false
output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30
```
(50 cycles; M35DT_LZ ran to nlim 6674 and its DTREP lines start at cycle 6630. floorcap_v / floorcap_rho / mesh_vol named for the rst-stored-key rule; new keys are added on the restart command line with a "### NOTE ... added" line.)

Return:
- dt per cycle (hst col 2) vs M35DT_LZ's 11.2 -> 5.9; the DTREP / DTREG lines: the binder must no longer be rho 1e-17 with |v| > 3e7 (floor gas now |v| <= 3e7);
- floorcap lines (n_v, dKE, dL_z net);
- hst Lz_tot (col 28) and Etot (col 19) per cycle; M35DT_LZ's values over the same cycles for comparison;
- s/cyc, CLEAN solves / BiCGStab its per cycle; FATAL / CLEAN WARNING counts.

## Arm FS1r (if a slot is free): restart bitwise with the key
- FS1 to rst cycle + 25 with a restart dump there (kept), then restart -> rst cycle + 50. hst md5 of the last 25 rows vs FS1. Delete the rst afterwards.

## Arm FS2: M35ZAMS (one node, 4 H200)
- Restart from the latest m35zams rst of the floorcap_lz run (Z3 / the PROD_RC3_1124 chunk), with the SAME layout and command line as that run, plus
  `rad_vet/floorcap_lz=true rad_vet/floorcap_lz_spin=true rad_vet/floorcap_mesh_vol=true rad_vet/floorcap_report=10 hydro/dt_report=10`, 50 cycles, no bin / rst outputs.
- Return the same items as FS1. Expectation: floor gas is capped at |v| <= 3e7 again, so dt returns to the floorcap_lz = false value (ZAMSFAST: 2.49 s = 0.3 h / (3e7 + c_s), floor gas at the cap). Floor gas at the cap can still bind M35ZAMS: ZAMSFAST's `hydro/dt_floor_courant` is the lever for that part, independent of this key.

## Decision after FS1 / FS2
If dt is back to the non-floor binder (m35mms) and Lz_tot / Etot hold to round-off, name `rad_vet/floorcap_lz_spin = true` wherever `floorcap_lz` is named (uniform settings; a USER decision for the production decks).

-- orion (coordinator)
