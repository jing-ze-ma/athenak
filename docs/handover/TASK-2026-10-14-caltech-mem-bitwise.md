# TASK raven -> Caltech: fin4-mem (MEMORY reductions, no method change) BITWISE vs FIN4 on IDORT dirs — SHORT, debug QoS

From Raven worker MEM, via the coordinator (USER: move short CPU checks to Caltech). Branch fin4-mem only changes allocation /
layout / buffer sizes; every result must be BITWISE equal to FIN4. CPU-only check (the GPU checks stay on Raven/Viper).
Raven runs the same 20 dirs (Raven job memidort 31149292, queued); whichever finishes first counts. Production NOT affected.

## 1. Binaries
- MEM: branch **`fin4-mem` @ `1f4e1724`** on github.com/jing-ze-ma/athenak (pushed). No kokkos change. Plain CPU build with
  EXACTLY the flags of your FIN4 / ENTHBE / mbband builds (gcc 13.2.0 + openmpi 5.0.1, Kokkos OpenMP, MPI on, Release,
  built_in_pgens). Report the md5.
- FIN4: your plain CPU build of fin4-1014 @ d505d21f from TASK-2026-10-14-caltech-fin4-suite (report its md5).

## 2. Keys: NONE. Decks AS IS + args.txt (the runner adds only time/nlim and an rst output if the deck has none).

## 3. Run (one debug job, 32 cores, ~10-15 min): `bash run_mem.sh <MEM> <FIN4> 4` (bundle docs/handover/open-1013/caltech_mem/)
- list_mem.txt = 20 dirs (10 Cartesian vet_sc, 10 sp vet_gd) with your list thread counts; dirs looked up in MB (mbband tree)
  then EB (ENTHBE tree); edit MB / EB if your trees moved. Writes <dir>/mem/ (O = FIN4, N = MEM, 10 cycles each), final rst
  data compared byte-for-byte after <par_end> (cmprst.py) -> <dir>/mem/cmp.txt; 3 dirs (RS) also MEM 5 -> -r -> 10 vs N.
- Expected: 20/20 BITWISE, 3/3 restart BITWISE.
- Optional if the debug job has time left (second debug job is fine): `LIST=list_b145.txt` (the 145 part-B dirs) via
  `bash run_mem.sh <MEM> <FIN4> 8 $PWD/list_b145.txt` (no restart arms beyond RS). Expected all BITWISE except the 2 known p2 FAIL.

## 4. Report: docs/handover/NOTE-2026-10-14-caltech-mem-bitwise.md
- md5s, job id, `COUNTS:` line, every non-BITWISE line (full), push mem_cmp.txt under docs/handover/caltech_mem/.

## Coordinator notes
- Priority: short (debug QoS). Do not delay running production.
