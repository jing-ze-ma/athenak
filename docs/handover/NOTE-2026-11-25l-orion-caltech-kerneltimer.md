# NOTE 2026-11-25l orion -> Caltech: kernel timer instead of the remaining HG ABBA (replaces the rest of 25k)

The HG arms scatter by about 2 % between identical binaries (NEW1 / NEW2 of 6d20f1337 differ by 5 %), so 30-cycle ABBA cannot decide 1-4 %. Please stop the HG ABBA after the arm now running (keep its table; post what is done) and use a Kokkos kernel timer instead. The sp arms (t9l62, m35mms) may finish.
- Tool: docs/handover/tools/kt.cpp (55 lines, no paths). Build: `g++ -O2 -shared -fPIC -o libkt.so kt.cpp` (host compiler that matches the MPI/Kokkos ABI). Run: `KOKKOS_TOOLS_LIBS=<dir>/libkt.so KT_OUT=<outdir>/kt` in the environment of every rank; one file per rank at finalize. It fences after each kernel: s/cyc is not comparable with the other runs, hst must still be printed-bitwise with the base of the same length. No nsys.
- **HG, 02-03 + 02-09, 8 ranks, hegiant.00205, HGN's command:** binaries 79c1cabef (9bebf6a3), 6d20f1337, 7b4141e1b, e8fca9eb8, e77765501, abcfc1e4a and 482a396df (69a4e783), each in two lengths (6 and 26 cycles); base repeated once at each length (so 16 runs).
- **sp: t9l62 and m35mms on one node each (02-16 and whichever is free), 79c1cabef and 7b4141e1b, 6 and 26 cycles.**
- Any run that hangs or fails: note it, go on with the next (timeout 10 min each). Return the kt files (text, gz) under docs/handover/caltech_results/KT/<arm>/ and one line per arm (rc, hst bitwise, CW). Orion evaluates; no tables are needed from you.
- Then continue 25c-25i's order; stop all by 08:00 PST.

-- orion (coordinator)
