# NOTE 2026-11-01 viper -> Caltech: divergence-free MHD remap (C32 -> C64/C128, optional new radial grid) is on fork/fin51-r-fix c62542349

Follow-up to NOTE-2026-11-01-viper-caltech-dhj-sample-seidel2026 (last bullet). Use it when you want to take an MHD run
(e.g. WASP-121b after the hydro spin-up and the MHD switch-on at C32) to higher resolution.

## What is on fin51-r-fix (fast-forward c8aac483c..c62542349, on top of your roche commits)
- **547652da3** merges dhj-remap-mhd: the pgen reader `problem/remap_bfile`.
  - It is read only when named, so decks without it embed exactly the same input as before.
  - docs/handover/scripts/dhj_remap.py already had all the features (+ --hfilter-deep, radial head); it rewrites the c64deep remap files bitwise.
- **c62542349** adds docs/handover/scripts/dhj_remap_HOWTO.md: hydro + MHD commands, keys, the seed, the transient caveat and the gate numbers.
- **Seed:** with remap_bfile set, the pgen dipole (bbot / bbot_gauss) is now SKIPPED, with a log line. B is the remapped field only. remap.athinput carries bbot_gauss over from the restart, so no hand edit is needed.

## Method
- Hydro: conservative MC-limited prolongation per parent; mass, e_int and L_z (restored per shell) exact.
- B: radial-gauge edge integrals of a vector potential built from the restart face fluxes, interpolated horizontally (linear). The code sets b = curl / area with its own areas. Coarse face fluxes are kept exactly, so div B stays at round-off.

## Gates (viper)
- **Restart bitwise, no remap keys:** WASP-121b MHD rot-320 rst, 300 cycles on MI300A. The whole restart is identical to c8aac483c (and to 3e0bd80a9 before the rebase).
- **CPU C8 -> C16:**
  - budgets: mass 0, e_int 0, L_z -1e-16;
  - div B 8e-16 at t0, 6e-15 after 40 cycles;
  - code mass +1e-10, from the Alfven-speed density floor in the top cells only.
- **GPU C32 -> C64 at rot 320:** div B 1.3e-14 at t0, 1.4e-14 after 30 cycles. C32 x 112 -> C128 takes 9 s, 2.2 GB RSS, and writes 0.6 GB.
- **CPU identity remap:** fresh start = source to 5e-14 (rho, E) / 5e-13 (B). After 10 cycles the deep cells match the -r path to 1e-11; the top cells differ, attributed to the fresh start resetting the ck caches (not isolated further).
- **Tests:** restart_bitwise, rad dhj ck/srclim and all mhd *_cpu tests pass. No new style errors.

## Commands
```bash
P=<repo>/docs/handover/scripts/dhj_remap.py
python3 $P <run>/rst/dhj.NNNNN.rst <new>/c64 --grid "mesh/nx2=64 mesh/nx3=64 meshblock/nx2=32 meshblock/nx3=32"
python3 $P <run>/rst/dhj.NNNNN.rst <new>/id        # identity control on the same path
# check: mass / e_int / L_z budgets < 1e-13, "new grid div B" ~ 1e-16
srun ... $BIN -i <new>/c64/remap.athinput -d <rundir> time/tlim=<t0+...>   # tlim MUST be overridden (embedded tlim = t0)
```
- For C128: `mesh/nx2=128 mesh/nx3=128`. For a new radial grid, add nx1, meshblock/nx1 and the stretch keys to --grid.
- Every command-line key on `-i` must already exist in remap.athinput. Later links are ordinary `-r` restarts.
- Smoke 10-30 cycles first.

## Caveat
- The remap raises the magnetic energy, by linear prolongation of grid-scale B_t: C32 -> C64 at rot 320 gave +2.5 % total and +19 % in the deep equatorial shells.
- That state still carried the pre-fix fossil noise; a clean state should show much less.
- Discard the first 0.5 rot (1.5 rot for deep magnetic budgets). Judge resolution against the identity control, not the old chain.
