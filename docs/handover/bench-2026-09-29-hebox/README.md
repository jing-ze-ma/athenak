# bench-2026-09-29-hebox: 3-D He box, implicit M1 + vet_sc, cross-cluster timing

> **HOLD LIFTED (user 09-29):** run it; see `docs/handover/TASK-2026-09-29-hebox-timing.md`.

Timing only. The same fresh-start run of the 3-D He-star FeCZ box with implicit grey M1 radiation
and the short-characteristics closure (`closure = vet_sc`, full tensor) is timed on viper
(MI300A, HIP), Caltech (H200, CUDA) and DeltaAI (GH200, CUDA). **One setting only: hesdirk2 at
cfl 0.3** (the accurate-phase setting; the cfl 0.9 arm was dropped, user 09-29). Full
instructions for the other machines: `docs/handover/TASK-2026-09-29-hebox-timing.md`.
The structure follows `docs/handover/bench-2026-09-28/` (WASP-121b).

| file | what it is |
|---|---|
| `inputs/hebox_bench.athinput` | the benchmark input (data paths = placeholder `HE_BOX_DATA`) |
| `run_bench_hebox.sh` | one fresh-start run: `run_bench_hebox.sh <binary> <1\|2\|4> <new run dir>` |
| `ana_bench_hebox.py` | reads the logs, prints one line per run |
| `bench_viper_hebox.sub` | the viper job script (reference for the other machines) |
| `RESULTS_viper.md` | the viper rows |

Quick use (inside a GPU allocation; the run dir must be new):

    MACHINE=viper|caltech HE_BOX_DATA=<.../he_box> ./run_bench_hebox.sh /path/athena 2 /scratch/hb/n2_r1
    python3 ana_bench_hebox.py /scratch/hb/*/bench.log

## Pinned commit

**401875f0** on fork/rt-integration (`401875f01beb8cd69c215eb096bbebd020e0886b`), PROBLEM =
box_convection, Kokkos = the in-tree submodule. It contains the M1 default flips of 09-28
(implicit_op_team_red / implicit_vimp_fold / implicit_halo_ovl_faces off, implicit_one_pass 0,
time2_lin_tol_fac 1, sp precond rbgs_fwd). The docs commit that adds this package is later; only
`docs/` differs, so the binary is the same.

## Grid and GPU counts

84 x 104 x 104 cells (x1 vertical, x2/x3 periodic), **4 MeshBlocks** of 84 x 52 x 52 (1 x 2 x 2).
Legal GPU counts with 1 rank per GPU: **1, 2, 4** (4 blocks, 2 or 1 per rank; 4 is the maximum).
vet_sc on several blocks/ranks (VetMBInit, `src/rad_m1/rad_m1_vet.cpp`) needs a uniform mesh,
periodic x2/x3 and a per-layer ray reach (w2 x w3 cells) no wider than a block (52 x 52): here
dx1/dx2 = 0.25, so a ray crosses 1-2 cells per layer; the ray set (64 rays, nmu 4 x nphi 8) is not
split over ranks (vet_mb_agroup = 1, the default, needs only ranks % 1 == 0). No vet_mb_* key is
set; all are at their defaults. Runs clean at 1, 2 and 4 ranks on viper (RESULTS_viper.md).

## The input: what changed, and why

Source: `docs/handover/caltech-2026-09-26/inputs/he_box/box.athinput` (= the 09-26 viper production
input `hebox_cfl_0926/inp/box.athinput` with the data paths replaced by `HE_BOX_DATA`), reconciled
with `/viper/ptmp2/jinma/hebox_cfl2_0927/inp/box.athinput` (the 09-27 cfl sweep with the ke-dt fix
keys). Every changed line is tagged `bench-2026-09-29-hebox` in the input.

| key | Caltech 09-26 | hebox_cfl2_0927 | bench | reason |
|---|---|---|---|---|
| `<rad_m1>/implicit_opac_update` | false | true | **true** | ke-dt fix (09-27); also the default with force_reference = wb_arad |
| `<rad_m1>/implicit_one_pass` | 8 | 0 | **0** | ke-dt fix; 0 is the default since 09-28 |
| `<rad_m1>/force_reference_work` | (absent) | split | **split** | ke-dt fix; the default auto resolves to split for box_convection |
| `<rad_m1>/implicit_precond` | mg_gf | mg_gf | **mg** | the current default on a single-level multi-D Cartesian mesh (plain mg fastest on H200 2 GPUs, handover 09-26 sect. 6) |
| `<rad_m1>/implicit_gf_modes` | 2 | 2 | removed | read only under mg_gf |
| `<rad_m1>/implicit_mg_fuse` | true | true | removed | no longer read by the code (obsolete key) |
| `<time>/nlim` | -1 | -1 | **1500** | fixed-length timing run |
| `<time>/tlim` | 1000 | 1000 | 1e6 | never reached |
| `<time>/cfl_number` | 0.3 | 0.3 (arms set it) | **0.3** | the timed setting |
| `<problem>/column_dump`, `rt_profile_dt/_file` | set | set | removed | diagnostics off (rt_profile default 0 = off) |
| `<output1>` hst dt | 1.0 | 1.0 | 10.0 | hst only |
| `<output2>`, `<output5>` (bin, dt 1e9), `<output9>` (rst) | present | output2/5 | removed | no bin/rst output during timing |
| data paths | HE_BOX_DATA | viper absolute paths | HE_BOX_DATA | the viper files have the same md5 as the tarball he_box/ files |

Unchanged (both sources agree): hesdirk2, closure vet_sc + vet_tensor full, implicit transport
with bicgstab, implicit_tol 1e-8, lhllc + plm, polytropic well-balancing, etotgrav, general EOS
table, fofc, gas Newton + EOS cache, the seed (vpert 1e-2 eint, seed 1234), the Marshak top and
imposed bottom flux. Every other `<rad_m1>` key written in the input equals the current default
at 401875f0 (checked key by key against the GetOrAdd defaults).

**Fresh start, not the production restart.** The viper production arms restart from an evolved
box (t = 3800, bench/m1_vet3dcfl_0923/C15), which cannot be shipped. The benchmark starts from
the static IC column plus the seed, so the flow is nearly at rest in the timed window; the
cycle cost of the implicit solve is close to the evolved box but not identical (compare the
viper row with hebox_cfl2_0927 R3: 92.7 ms/cycle, 5.06 Picard passes/solve, 1 GPU, ROCm 6.3).

## Timing window and outputs of ana_bench_hebox.py

- A run is 1200 cycles (`time/nlim`), diag line every 10 cycles. On viper 2 GPUs that is ~2 min
  of stepping (102 ms/cycle in the smoke).
- **Window: cycles 400-1200** (`ana_bench_hebox.py` default, `--c0/--c1` to change). dt is set by
  the hydro sound speed at the bottom: 0.16142 -> 0.16121 over the smoke's 200 cycles. ms/cycle per
  20 cycles stays within 95-107 ms from cycle ~20 on (smoke, 2 GPUs); the first step is a
  backward-Euler step (hesdirk2 start). Starting
  the window at 400 leaves a wide margin (RESULTS_viper.md compares it with 100-1200).
- Printed per run: ms/cycle (median of the 10-cycle windows, and the whole-window mean),
  wall s per simulated s (same two ways), mean dt, and from the `<rad_m1>` summary lines at the
  end of the run (whole run, cycles 0-1200): Picard passes per implicit transport solve, BiCGStab
  inner iterations per transport solve and per linear solve, NON-CONVERGED solves, the sum of the
  fallbacks (positivity, vimp positivity, hesdirk2 stage, line_jacobi), FATAL and nan counts.
  hesdirk2 does 2 implicit solves per step (399 solves in 200 cycles).

## Smoke

Viper job 12018403 (apudev, 2 GPUs, ROCm 7.2 build md5 fe495394fa25b6ecc89f67c3e9384044):
200 cycles, ndiag 1, `rad_m1/implicit_picard_log=20`, rc 0. FATAL 0, nan 0, NON-CONVERGED 0,
positivity / vimp / stage / line_jacobi fallbacks 0, bicgstab breakdowns 0, vet clips 0,
vet_tensor guard 0; Picard 6.76 passes per solve (max 10), 6.75 BiCGStab iterations per linear
solve; 102.1 ms/cycle (cycles 50-200), dt 0.1613.
Command for the other machines: `EXTRA="time/nlim=200 time/ndiag=1 rad_m1/implicit_picard_log=20"`,
then `python3 ana_bench_hebox.py <run>/bench.log --c0 50 --c1 200`.
