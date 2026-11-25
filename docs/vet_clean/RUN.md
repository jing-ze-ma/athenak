# Running the clean VET scheme

## 1. A first run (CPU, one minute)

The test decks in `tests/vet_clean/decks/` are complete, small inputs. The closed rotating
scattering shell on a 3-D spherical-polar wedge is a good first run:

```bash
mkdir -p run/shell && cd run/shell
OMPI_MCA_hwloc_base_binding_policy=none OMP_NUM_THREADS=4 \
  ../../build_cpu/src/athena -i ../../tests/vet_clean/decks/shell_sp.athinput
```

It runs 40 cycles in about 15 s and writes:
- `sp.hydro.hst` (gas totals) and `sp.user.hst` (the energy and angular-momentum budget,
  see [OUTPUTS.md](OUTPUTS.md));
- `bin/sp.vet.*.bin` (radiation E and F) and `bin/sp.hydro_u.*.bin` (gas conserved
  variables), readable with `vis/python/bin_convert.py` or `tests/vet_clean/vett.py`.

Any key can be changed on the command line, `block/key=value`, for example
`time/nlim=100 mesh/nx1=64 meshblock/nx1=64`.

The log ends with a summary of the radiation solver (solves, iterations, fallbacks, floor
clips); [OUTPUTS.md](OUTPUTS.md) explains which lines matter. A healthy run prints no
`CLEAN WARNING` line.

## 2. Production star runs

Templates: `inputs/stars/star_bsg.athinput`, `star_m35mms.athinput`,
`star_t9l62.athinput`. They hold the full production key set (tile layout mode). They need:

1. a binary built with `-D PROBLEM=star_vet` (see [BUILD.md](BUILD.md));
2. the opacity tables and initial-condition profiles named by the `star_opac_table`,
   `star_planck_table`, `star_ic_file` (and for the BSG `star_ic_profile_file`) keys. Replace
   `TABLES/` with the directory that holds them. The tables are not in the repository.

Typical launch, one MPI rank per GPU:

```bash
srun -N 2 -n 8 --gpus-per-node=4 env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 \
  ./athena -i star_t9l62.athinput -d run_t9l62 time/tlim=...
```

The decks use MeshBlocks of 128 x 128 lateral cells spanning the whole radial range (one
MeshBlock along r, which tile mode requires): the BSG deck is 256^3 cells in 4 MeshBlocks, the
m35mms deck 336 x 256 x 256 in 4, the t9l62 deck 512 x 256 x 512 in 8. To use more GPUs, make the lateral
MeshBlock smaller and set `<mesh>/layout_tile` to the new lateral MeshBlock size. A run
with tiles of N x N cells gives the same answer, bit for bit, for every layout of whole-tile
MeshBlocks. Different N give answers that differ at the solver tolerance.

On a phi-periodic wedge the width must be a multiple of the ray set's azimuthal period (GL
7 x 18: 20 deg; LQ10: 90 deg; the default icosahedral set: 180 deg), else the run stops with a
FATAL; a narrower wedge names `vet_gd_wall_sym = true` with `vet_gd_wall_interp = true`.

The templates name `<rad_vet>/clean_verbose = 1`, which prints the per-solve `CLEAN`
diagnostic lines (useful context for any `CLEAN WARNING`); your own decks can leave it at the
default 0 ([OUTPUTS.md](OUTPUTS.md), section 4.2).

On the first cycle the code relaxes the radiation field to its steady state in the initial
gas (`clean_init = steady`). The log line `clean_init steady pass 1 ...` is a convenient
fingerprint: a new binary, machine or compiler must reproduce it to a few digits.

### Restarts

A restart (`-r file.rst`) reads the whole input deck that is stored in the restart file.

> **Rule: any key you do not give on the command line of a restart takes the value stored
> in the restart file, not the code default.**

This bites when a deck or the code default changed since the restart file was written. For
example, an old restart file stored `time2_k1_stiff = 1.0`, and a restart that relied on the
new default 0 silently ran with 1.0, which changed the photospheric luminosity by about
50 %. So, when you restart from a file written with other keys (another code version,
another run), name every key that matters on the command line:

```bash
./athena -r run/hegiant.00205.rst rad_vet/time2_k1_stiff=0 rad_vet/clean_stag_fwb_stage=false \
  rad_vet/clean_stag_fwb_t2=false rad_vet/clean_w_fse_tau=true rad_vet/time2_gtime_mom=true \
  rad_vet/clean_pred_stage=false rad_vet/clean_stag_central=true rad_vet/clean_stag_kpat=false \
  rad_vet/clean_stag_eicarry=true rad_vet/clean_gas_secant=max_stage2 rad_vet/clean_pred_xch=true \
  rad_vet/clean_stag_fcone_gas=false rad_vet/clean_wb_exch_local=true hydro/floor_pay=rad_column \
  hydro/floor_pay_skip=true hydro/dfloor_flux_limit=true hydro/rad_signal_speed_diff=true
```

and check the start of the log: it prints the value in effect for every key (grep for the
keys you named). A restart of your own run with an unchanged deck needs no extra keys.

Restarts are bitwise: N cycles straight and N/2 + restart + N/2 give identical output (the
`restart_*` tests check this). The restart file holds every piece of solver state that
carries from step to step (the stage slopes, the ray-closure history, the ledgers).

A key the restart file does not hold takes the fresh default when it corrects the scheme
(the clean-scheme keys, `time2_*`, `merge_lat_mom`, `fofc_flag_exchange`, `clean_cs_lat_wb`
...), with a `### NOTE: ... the fix applies`; naming the old value keeps the old behaviour. A
few keys (the angular set and its rotation and averaging, the closure history, the problem
set-up) keep the restart default; the start-up log echoes it with the reason.

**Restarting a file of an older scheme.** Files written before a key existed store neither it
nor the keys it goes with, so the run changes at the restart: the NOTEs list what applies and
how to keep the old value. A stored value always wins, also an old one that a newer key needs
differently; the run then stops with a FATAL that names both keys and what to name. Do not
restart such a file with nothing named: name every key of the production deck (the list of
[KEYS.md](KEYS.md), "Keys that must always be named"), read the NOTEs and check the start of
the log before you let it run on.

### Time step and retries

`time/dt_retry = true` (production) lets a step that fails its checks (no convergence, a
negative energy) be repeated with a smaller dt (`time/dt_retry_fac*`), and dt grows back
after `dt_retry_ngood` good steps. The retry is bitwise: a retry with factor 1 repeats the
same step exactly (the `dt_retry` test). Each retry prints a `DT_RETRY cycle=...` line.

## 3. Keys

[KEYS.md](KEYS.md) lists the keys a user sets, with meaning, default and production value,
the restart rule, and the keys that must always be named.

## 4. Checking a run

- Energy: see "Checking energy conservation" in [OUTPUTS.md](OUTPUTS.md).
- `CLEAN WARNING` lines and the end-of-run summary: [OUTPUTS.md](OUTPUTS.md).
- The test suite: `tests/vet_clean/run_tests.sh --exe <binary> --quick`.

## 5. Why is my run slow? (the time-step report)

```bash
mpirun -np 4 ./athena -r rst/star.00010.rst -d . time/nlim=<cycle + 100> hydro/dt_report=10
grep -E "^DTREP|^DTREG" log | tail -20
```

`<hydro>/dt_report = N` (read only when named, no effect on the state) prints every N
cycles which limit sets dt: a `DTREP` line with the cell of the global minimum (rank, MeshBlock,
k j i, r, rho, T, sound and flow speeds, the radiation terms) and, on a spherical-polar mesh,
`DTREG polar` / `DTREG rest` for the minima near the poles and elsewhere. If the minimum sits
at the poles or in hot floor gas, change the mesh or the floors there rather than the
radiation keys. `<mhd>/dt_report` does the same for MHD runs
(radiation with MHD is not tested: SCHEME.md).

## 6. GPU memory on spherical-polar wedges

If an A100 run with `vet_gd` on a spherical-polar wedge stops with an out-of-memory error
(often at a ray-cut change, `vet_vgdgb_keep`), name

```
<rad_vet>
vet_gd_band_ring = true
```

It lends the band intensities per shell from a few wide slabs instead of keeping every
shell's band: the same results bit for bit, less memory, a few per cent slower. Unnamed, it
is on where it applies except on the raw band of spherical-polar meshes.
