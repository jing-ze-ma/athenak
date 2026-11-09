# NOTE 2026-11-09 Orion -> viper: status of TASKs 5-7 (raycheap2, hgblend, raycheap3)

**raycheap2: runs done. The full NOTE (NOTE-2026-11-09-orion-raycheap2.md) follows from the agent. Headline numbers:**
- R1 restart (FIX 4e9f602c8): bitwise m1 4/4, m1_lat 12/12.
- R2 (BASE f34ada236): differs, max |diff| 1.191e-3. This is the bug.
- R3 FIX vs BASE: bitwise on H and P.
- R4 FIX 4 vs 1 rank: bitwise.
- -> The rstabe-1109 fix closes the G3 restart gap.
- C1 (Z_Q64, as is):
  - d 8.680e-5, 3.029e-5, 7.497e-6; p 1.52, 2.01.
  - It reproduces viper's Raven CPU numbers (8.7e-5, 3.0e-5, 7.5e-6) across machines.
- C2 (Z_Q128T8 with implicit_tol 1e-13, implicit_lin_tol 1e-14):
  - d 6.975e-6, 5.135e-6, 2.850e-6; p 0.44, 0.85.
  - This is essentially unchanged from the 1e-10 run (7.2e-6, 5.1e-6, 2.9e-6; p 0.49, 0.84).
  - -> The d ~ 3e-6 floor is NOT the solver tolerance.
- All runs rc=0, 0 FATAL. Fallbacks occur only at the coarsest level: C1 2 and C2 3 stage fallbacks, with 6 and 8 BE steps.

**hgblend (1f314023a): running.**
- REF reuses the hgpos b798f1a30 build.
- All 7 arms run at once (O1 in its own 120G allocation).
- NOTE-2026-11-09-orion-hgblend.md in about 20-30 min.

**raycheap3: running.**
- It reuses the 9fbd2d3c6 CPU build: 9b9951606 has an identical src/, CMakeLists and kokkos (git diff --stat is empty).
- Deviation for speed (USER: "make tests fast"): all 5 arms x 4 levels were submitted at once, including c0.01875 for every arm. conv.py is reported on 3 and on 4 levels.
- C2/C3 hit the rbgs / krylov_fuse stop and were resubmitted with the fuse workaround. Those two are queued behind Orion's 6-job user limit and start when hgblend ends.
- NOTE-2026-11-09-orion-raycheap3.md in about 1-1.5 h.
