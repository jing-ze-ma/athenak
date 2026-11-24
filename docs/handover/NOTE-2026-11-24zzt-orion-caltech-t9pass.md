# NOTE 2026-11-24zzt orion -> Caltech: t9l62 per-cycle profile (T9PASS): zero-cost print now, profile arms later

**Now (zero cost):** section 0 below, at t9l62's next restart (together with 24zzp). **Later, only on a node freed by bsg's pause / end, AFTER the Compton arms (24zzs) and only if it doesn't delay bsg's resume or CSHERO:** PE, then PL, then PO. Never on production nodes.

## 0. Zero cost, production T2 t9l62
- At its next restart (with the halved rst cadence of 24zzp), also name `rad_vet/clean_verbose=2`.
  - This prints one `CLEAN PASS` line per pass (it, de, dv, tol, fe, fv, clos / opac). It is a host printf of host scalars:
    hst bitwise (Orion gate f676a91b3, 11-24 02:40). The release binary already has it.
  - It costs ~10-14 lines per cycle in the log. Keep the log, or grep `^CLEAN PASS` into a side file each chunk.
- From then on, every phase (outburst, D, E/F) gives its pass profile without a separate arm.

## 1. Arms (one free node, 4 H200; never on the production nodes while production runs)
Binary for all: the T2 binary `athena_ia123rc4v_262affd6f` (md5 7d7c4a71). The release a23a78794 has the same src.
Tool: `profc.sh` from H200PROF (docs/handover/caltech_runs_1121/h200prof/), NCU_ON=0 (nsys phase only).

| arm | rst | keys | N0 -> N0+NC | purpose |
|---|---|---|---|---|
| PE | T2 t9l62's latest rst (early, t < 4.5e5) | T2 restart args + `rad_vet/clean_verbose=2` | +25 (nsys window cycles 2..) | early cycle: kernel families, GPU idle, per-pass lines |
| PL | prod_fresh_1121b/t9l62/rst/t9l62.00031.rst (t ~1.34e6, J4 chunk 0, cfl 0.3 / conv 1e-5) | **every** `<time>` / `<hydro>` / `<rad_vet>` key of the T2 deck as block/key=value (awk line of 24zg), + `mesh/layout_invariant_mode=tile mesh/layout_tile=128 time/cfl_number=0.3 rad_vet/clean_verbose=2` | first a plain +40 (no nsys) to leave the scheme change behind, write one rst at its end (job/final_restart=true); then PL from that rst +25 under nsys | late cycle (dt ~96): same tables |
| PO | T2 t9l62's first rst with 4.6e5 < t < 5.2e5, once it exists | as PE | +25 | outburst cycle (dt 12-40) |

- PE / PO / PL each take ~10 min. PL's warm-up takes ~5 min.
- Check run.log for the `added from the command line` NOTEs (PL) and grep `cfl` / `clean_conv_max`.
- No rst is kept (except PL's warm-up rst, deleted after PL).

## 2. Return per arm
- `res.txt`, `n/kname_r0..3.txt`, `n/gaps_r0..3.txt`, `n/stats_r0.txt` (text only).
- `grep -E "^CLEAN (PASS|cycle|GNSW|WARNING)" n/run.log` (all lines of the window).
- s/cyc over the window, plain (no nsys) s/cyc of PL's warm-up over its last 20 cycles.
- dt min / median / last; CW / fallbacks / DT_RETRY counts.

## How Orion reads it
- Families with stag2/h200prof/fam.py: closure (rays), Picard RHS + face fluxes, F solve, gas Newton, E pc, BiCGStab, hydro,
  halo, FCT, rest; fixed per cycle vs per pass (kernel counts / passes).
- Per-pass lines: which of de / dv holds the late stages open, at which pass, and against which tol (budget vs strict, and the
  opacity passes 1-2 / closure pass 1 / xdonor freeze from pass 4 / gas Newton switch from pass 6).
- Orion runs the same tables on a reduced CPU deck (512x32x32), to test cuts before they come to Caltech.

-- orion (coordinator)
