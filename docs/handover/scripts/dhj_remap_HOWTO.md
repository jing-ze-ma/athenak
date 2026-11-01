# HOWTO: remap a WASP-121b (deep_hot_jupiter_rt, cubed sphere) restart to a finer grid

Tool: `docs/handover/scripts/dhj_remap.py` (python3 + numpy + scipy). Code side: the pgen
keys `problem/remap_file` (hydro state) and `problem/remap_bfile` (MHD face fields) of
`src/pgen/deep_hot_jupiter_rt.cpp`. Both keys are absent by default; a run without them is
bitwise the same as before (gate below).

What can change: the cubed sphere C(n) -> C(f n) (integer f, e.g. C32 -> C64 / C128), the
MeshBlock layout (any meshblock/nx2, nx3 dividing the new nx2, nx3; one MeshBlock in r),
and the radial grid (nx1, the polynomial stretch keys, a raised x1max; x1min is fixed).

## 1. Binary

Any build of this branch with `-D PROBLEM=deep_hot_jupiter_rt` (GPU for the run).

## 2. Remap file and input (login node; C32 x 112 -> C128 took 9 s, 2.2 GB RSS, 0.6 GB out)

```bash
RST=<run>/rst/dhj.NNNNN.rst                  # read in place, never written
P=<repo>/docs/handover/scripts/dhj_remap.py
# horizontal C32 -> C64 (MeshBlocks 112 x 32 x 32), same radial grid
python3 $P $RST <new>/c64 --grid "mesh/nx2=64 mesh/nx3=64 meshblock/nx2=32 meshblock/nx3=32" \
  | tee <new>/c64.log
# identity control: same grid, same fresh-start path (run it next to the remapped run)
python3 $P $RST <new>/id | tee <new>/id.log
# radial grid too: add the planet_setup.py override string (nx1, meshblock/nx1, stretch
# c1..c8, bumps; G_PROD in grid_w121*.env) to --grid
```

- Output: `remap.dat` (hydro: rho, v, e; E - e_mag for MHD), `remap_b.dat` (MHD only) and
  `remap.athinput` = the restart's embedded input with the new grid keys,
  `time/start_time` = t(RST), `problem/remap_file` and `problem/remap_bfile` set. Output
  counters and every `last_time` carry over; nothing else changes.
- Read the printout:
  - budget table: mass, e_int and L_z must change by < 1e-13 (round-off). KE rises by
    ~1e-3 (sub-cell velocity with momentum conserved), sum rho Phi changes by ~5e-4
    (quadrature), E_tot by ~1e-5. Expected.
  - `B: new grid div B ... ~1e-16`.
  - `B gate (c)`: magnetic energy old -> new per pressure level (see the caveat below).
- MHD options: `--bh-interp linear` (default; each child gets 1/f^2 of the parent radial
  flux, no overshoot) or `cubic` (adds sub-cell structure: C32 -> C128 ME +2.4 %, |B|max
  +58 %; not recommended). `--no-bfield` writes the hydro file only (continue with the
  analytic seed instead).
- Restrictions: cubed sphere only, restarts with RT warm-start or M1 slabs are refused,
  x1min fixed.

## 3. Run: a FRESH start (`-i`, not `-r`) in a NEW run directory

```bash
export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1
srun -n <2 x nodes> -c 24 $BIN -i <new>/c64/remap.athinput -d <rundir> [keys] [-t hh:mm:ss]
```

- Every command-line key given to a fresh `-i` start must exist in `remap.athinput`
  (add it there as `key = value` first, otherwise the run stops).
- Seed: with `problem/remap_bfile` set the dipole seed (`problem/bbot` / `bbot_gauss`)
  is NOT applied; the pgen prints one line saying so. B is the remapped field only.
- The first outputs are written at t0; rst/bin numbering continues from the restart,
  hst/log start fresh. ck caches, warm starts, flux-history accumulators and floor
  bookkeeping start fresh, as on any fresh start. Later links are ordinary `-r`
  restarts of the new run (the remap keys are ignored on a restart).
- Smoke first (apudev, ~10-30 cycles, same binary, input and keys).

## 4. Transient caveat (MHD): discard 0.5-1.5 rotations, compare with the identity control

- The remap keeps every coarse face flux exactly, but the magnetic ENERGY is not
  conserved: the code's cell energy uses face-averaged B, and the CT-consistent linear
  prolongation raises it by (L-R)^2/16 per child pair where opposite faces of a coarse
  cell differ. C32 -> C64 at rot 320 (c64deep_1031): total +2.5 %, the deep grid-scale
  field (shells 2-10) +19 % at the equator. The excess decays: +2.5 % at 0.75 rot,
  +0.4 % at 1.5 rot against the C32 run.
- So discard the first 0.5 rot at least (1.5 rot for the deep magnetic budget), and
  judge resolution effects against the identity-remap control (same fresh-start path,
  same reset caches), not against the old restart chain.
- Hydro: KE +1e-3 at t0 from the sub-cell velocity profile; nothing else to discard.

## 5. Gates of the merge onto fin51-r-fix (remapmhd_1101, 11-01)

- Tool: the merged dhj_remap.py writes remap.dat / remap_b.dat bitwise equal to the
  c64deep_1031 files (rot-320 C32 restart, identity and C64).
- Off = bitwise: no remap keys, 300 cycles from the rot-320 MHD restart on MI300A: the
  whole restart after `<par_end>` is identical to the parent (vs 3e0bd80a9, apudev job
  12526590; after the rebase onto c8aac483c, job 12527882).
- CPU, C8 x 112 MHD (analytic start + 20 cycles, then remap):
  - identity remap, fresh start vs the source restart at t0: rho, E 5e-14, B 5e-13 of
    |b|max. After 10 more cycles vs the `-r` path: 1e-11 in the deep cells (i < 24), up
    to 1e-3..1e-1 in the upper, RT-dominated cells, from the reset ck caches of a fresh
    start (not the remap).
  - C8 -> C16: tool budgets mass 0, e_int 0, L_z -1.1e-16; div B 8e-16 at t0, 6e-15
    after 40 cycles. Code mass at t0 +1.1e-10: the Alfven-speed density floor
    (valfven_max_cgs) raises rho in the top 26 of 112 cells; below, rho equals remap.dat
    exactly.
- GPU, C32 -> C64 at rot 320: div B 1.3e-14 at t0 and 1.4e-14 after 30 cycles.
