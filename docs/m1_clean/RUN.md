# Running the clean M1 scheme

## 1. A first run (CPU, one minute)

The test decks in `tests/m1_clean/decks/` are complete, small inputs. The closed rotating
scattering shell on a 3-D spherical-polar wedge is a good first run:

```bash
mkdir -p run/shell && cd run/shell
OMP_NUM_THREADS=4 ../../build_cpu/src/athena -i ../../tests/m1_clean/decks/shell_sp.athinput
```

It runs 40 cycles in about 15 s and writes:
- `sp.hydro.hst` (gas totals) and `sp.user.hst` (the energy and angular-momentum budget,
  see [OUTPUTS.md](OUTPUTS.md));
- `bin/sp.m1.*.bin` (radiation E and F) and `bin/sp.hydro_u.*.bin` (gas conserved
  variables), readable with `vis/python/bin_convert.py` or `tests/m1_clean/m1t.py`.

Any key can be changed on the command line, `block/key=value`, for example
`time/nlim=100 mesh/nx1=64 meshblock/nx1=64`.

The log ends with a summary of the radiation solver (solves, iterations, fallbacks, floor
clips); [OUTPUTS.md](OUTPUTS.md) explains which lines matter. A healthy run prints no
`CLEAN WARNING` line.

## 2. Production star runs

Templates: `inputs/rad_m1_clean/star_bsg.athinput`, `star_m35mms.athinput`,
`star_t9l62.athinput`. They hold the full production key set (tile layout mode). They need:

1. a binary built with `-D PROBLEM=he_star_m1` (see [BUILD.md](BUILD.md));
2. the opacity tables and initial-condition profiles named by the `he_opac_table`,
   `he_planck_table`, `he_ic_file` (and for the BSG `he_ic_profile_file`) keys. Replace
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

The templates name `<rad_m1>/clean_verbose = 1`, which prints the per-solve `CLEAN`
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
./athena -r run/hegiant.00205.rst rad_m1/time2_k1_stiff=0 rad_m1/clean_stag_fwb_stage=false \
  rad_m1/clean_stag_fwb_t2=false rad_m1/clean_w_fse_tau=true rad_m1/time2_gtime_mom=true \
  rad_m1/clean_pred_stage=false rad_m1/clean_stag_central=true rad_m1/clean_stag_kpat=false \
  rad_m1/clean_stag_eicarry=true
```

and check the start of the log: it prints the value in effect for every key (grep for the
keys you named). A restart of your own run with an unchanged deck needs no extra keys.

Restarts are bitwise: N cycles straight and N/2 + restart + N/2 give identical output (the
`restart_*` tests check this). The restart file holds every piece of solver state that
carries from step to step (the stage slopes, the ray-closure history, the ledgers).

Some keys have one default for a fresh run and another for a restart of a file that does
not name them, so that old runs continue unchanged. The start-up log echoes such
defaults with the reason, e.g. "(default with vet_gd on a fresh start)".

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
- The test suite: `tests/m1_clean/run_tests.sh --exe <binary> --quick`.
