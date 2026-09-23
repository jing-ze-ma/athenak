# runs_4h_bcadv: implicit_bc_advect becomes the default (branch m1-bcadv)

Base 5c9e4432 (rt-integration). Runs, binaries, inputs and scripts:
/viper/ptmp2/jinma/bcadv_0923 (`measure.sh`, `redo.sh`, `gates.sh`, `run.sh`,
`run_t7.sh`, `analyze.py`, `cmp_payload.py`). CPU, gcc + openmpi, login node, nice 10.

## What changed

`<rad_m1>/implicit_bc_advect` (tests_m1/runs_4f_drift: a Marshak x1 end face also
carries the upwinded enthalpy flux A E) now defaults to **true**.
- Always read (GetOrAdd) and echoed into restarts, as before.
- On a restart (`global_variable::restart_run`) a file that lacks the key (written
  before 5c9e4432) keeps the old value false. Files written by 5c9e4432 echo false.
- Explicit input values override. The value is printed in the implicit-transport banner
  when an x1 end is Marshak.

## Measurement: He slab / box, false vs true (HEAD 5c9e4432 binary, key in the input)

Seeded 2-D He slab (`slab2d_nd` of runs_4b_defaults, vpert = 1e-2; top x1 = marshak with
E_bath = 0, bottom = flux). Relative differences true vs false at the end unless noted.

| quantity | Edd 200 s (1 rank) | vet_sc 200 s (2 ranks) | Edd 1000 s (1 rank) | 3-D box, 60 cycles (4 ranks) |
|---|---|---|---|---|
| top-cell <E> ratio - 1 | -1.5e-6 | -2.0e-6 | -1.5e-6 | -1.1e-6 |
| top-cell <v1> (2nd half) | 1.95e4 cm/s (v/c 6.5e-7), unchanged | 2.01e4, unchanged | 1.88e4, unchanged | n/a (no profile dump) |
| F1top/Fin mean, 2nd half | 1.000002 -> 1.000001 | 1.000004 -> 1.000003 | 1.000001 -> 1.000000 | 1.000008 -> 1.000007 |
| KE1 / KE2 | 2.5e-6 / 1.5e-6 | 3.9e-6 / 1.8e-6 | 1.4e-5 / 2.9e-6 | 1.3e-5 / 2.2e-9 |
| total energy (hydro) / rad Etot | 1.3e-9 / 1.2e-8 | 8.5e-10 / 1.9e-8 | 3.0e-10 / 1.2e-8 | 9.3e-10 / 8.1e-9 |
| T(z) max rel diff, top 10 cells (last profile) | 3.9e-7 | 3.2e-7 | 4.1e-7 | n/a |
| NON-CONVERGED / pos. fb / line_jacobi fb | 0/0/0 both | 0/0/0 both | 0/0/0 both | 0/0/0 both |
| Picard mean false -> true | 2.0483 -> 2.0492 | 2.6232 -> 2.6240 | 2.3751 -> 2.3754 | 2.2333 -> 2.2333 |

The top cell holds c q E / F1top = 1.0000 in both arms; the O(v/c) = 6e-7 shift of
E_top is the expected size. The change is at the 1e-6 level: default switched on.

## Gates (new binary vs HEAD 5c9e4432 binary)

| test | result |
|---|---|
| explicit false vs HEAD (no key): slab Edd 1 rank, vet_sc 2 ranks, box 4 ranks | bitwise (hst, bin payload from `<par_end>`) |
| explicit false vs HEAD: T7 M5 N512 be (runs_4f_drift input, names false) | bitwise (tab) |
| new default (no key) vs HEAD + true: slab Edd, vet_sc, box | bitwise (hst, bin payload) |
| new default (no key) vs HEAD + true: T7 M5 N512 be | bitwise (tab) |
| d0c59f7c restart file (no key), t = 20, 2 ranks: HEAD vs new binary | bitwise (hst, bin, rst); new echoes `implicit_bc_advect = 0` |
| HEAD restart file (key = false) continued by the new binary vs HEAD full run | bitwise (bin, rst, 20 hst rows) |
| new default full run to t = 40 vs restart at t = 20 (2 ranks) | bitwise (bin, rst, 20 hst rows); echoes `= 1` |
| cpplint rad_m1_implicit.cpp | only the pre-existing build/include_subdir errors |
