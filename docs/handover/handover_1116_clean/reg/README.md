# clean-cc pre-push regression harness

This is the CPU matrix that Orion ran before every clean-cc push (pushes 7-14). It needs
python3 with numpy and scipy, and it takes about 2 min on one node with 72 or more cores
(29 runs, 30 at a time).

## Run

```
cd reg                                   # work in a copy; outputs go under it
./run.sh NEW /path/to/new/athena         # arm NEW (job.sbatch: ARM=NEW BIN=... ; ~2 min)
./run.sh OLD /path/to/baseline/athena    # same for the last pushed head
./ana.sh NEW OLD                         # scores, iterations, warnings, conservation
./cmp.sh NEW OLD                         # files that differ, total CLEAN WARNING of NEW
./wb.sh /path/to/new/athena              # well-balancing levels (~1 min)
grep '^RT' <log of run.sh>               # restart subset
```

Build: Release, `-ffp-contract=off`, MPI+OpenMP, `-D Athena_M1_VGD_NMAX=4096`. All runs are
OpenMP only (no mpirun), so the bitwise checks also hold across machines with the same
compiler. Set ATHENA_VIS to a vis/python with bin_convert.py; the default is cc/ana.
`ref_62d71accc.txt` holds the ana/wb/RT output of push 13 (62d71accc, bitwise to
0e2194db2). Compare it with `./ana.sh NEW`; cpu times are not reference values.

## Gates (cmds.txt)

| gate | deck | checks |
|---|---|---|
| hcyl cbg cb shd nltek1 | cc/decks, run_ladder.sh x0.5 (shd x0.25, nltek1 x1) | Cartesian L1 vs exact (score7.py); shd umbra, cbg FWHM via cc/ana/score.py |
| thin_p thin_r thin_f | gtg.py thin (sp wedge, gd closure; p: periodic theta, f: full sphere) | L1(E), L1(F) vs the thin point source |
| dif lat | gtg.py dif / lat | diffusion and lateral-structure L1(E) |
| pulse0 | gtg.py pulse0, reflecting x1max | closed-box energy (energy.py), round-off |
| s88 s81 | clean_init = steady, sph | steady start converges, 0 warnings |
| rhse | rhse/rhse_sp.athinput | steady-start passes and E, max v/c_s (no WB) |
| cs uni closed thin1 thin1mb feed1 | cs/*/cs.athinput (cubed sphere) | uni exact; closed sum E V; thin1 L1 (also 2 radial blocks); feed1 E_rad+E_gas and L_z round-off |
| flz q / Q | fb/fb.py flz, chi=1 rotating box, fv1 | L_z(gas+rad) and E drift round-off, iterations |
| rw x4 | rw/*/rw.athinput, amp +-1e-3, N32 / N64 | odd (linear) and even parts vs baseline (bitwise unless physics changed) |
| RT hcyl sp3mb feed thin1 | rst/rt.sh: 16 cycles vs 8 + restart + 8 | `hydro_u same m1 same` (bitwise restart; the predictor-increment warning text is stale) |
| wb fv1 fv lld cart | wb.sh: rhse at rest, force_reference = clean + hydro WB, t = 4 | max v/c_s at solver tolerance: ~8e-8, 7e-8, 6e-8, 1.9e-7 |

## Pass criteria

- rc = 0 everywhere, and **zero CLEAN WARNING** in every run (`cmp.sh` total, the `warn`
  column of `ana.sh`, wb.sh).
- RT: all four `hydro_u same m1 same`.
- Pure refactors, scheduling, GPU or initialiser fixes: `cmp.sh NEW OLD` gives `differ 0`
  for every run directory.
- Physics changes: state which decks change and why. v = 0 decks must stay bitwise unless
  the change touches them. Conservation lines stay at round-off (pulse0 energy ~1e-15,
  closed/feed1 rel <= 2e-16, flz L_z drift ~1e-15). rw changes shrink with resolution
  (N64 vs N32) when they are truncation-level.
- WB levels stay at the reference (a jump to 1e-6..1e-3 is a regression; push 9 had one).

## Files

run.sh (per-arm deck copies, then cmds.txt with $REG $BIN $ARM), job.sbatch (set partition,
modules), ana.sh, cmp.sh, wb.sh; gtg.py + geo/gt.py (sp wedge), fb/ (flz), energy.py,
cs_tests/ct.py, cc/ (run_ladder.sh, score7.py, decks, ana: exact solutions at 128 cells,
others by `python3 hcyl.py exact` / shd.py's cache), rst/ (rt.sh, cmpbin.py, decks), wb/.
