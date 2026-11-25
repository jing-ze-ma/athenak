# NOTE 2026-11-25z orion -> Caltech: the merged next line e8d3e64d9 (memory + slab + speed, trimmed decks), CUDA gates

Branches: caltech-next2-e8d3e64d9 = 79c1cabef + trimmed star / suite decks + cc_split (HIP) + the memory series (slab included) + the speed series 1d10fc373; caltech-skip-ab9994907 = e8d3e64d9 + one commit (fill passes skipped if no hole). CPU suites green (55/55 at bf007230a; e8d3e64d9's own suite and the poison arms are running on Orion and HIP arms are bitwise). On free nodes, in this order, after 25y:
1. Build star_vet and the suite of e8d3e64d9 (recipe as before, tiny gate 3.686e-11, .prov); **full CUDA suite** (55 tests, 4 lists, the two long tests alone; new tests star_cs_ring, ray_cut_slab). Values against 79c1cabef's run where they share a test. First FAIL at once.
2. **Restart chain m35zams** as for 79c1cabef: +1, 25 + rst + 25 vs 50, vs 79c1cabef 50 cycles (bitwise payloads expected).
3. **HG, 8 H200, T6, ABBA 26 cycles 79c1cabef vs e8d3e64d9** (hst bitwise) with nvidia-smi peak per GPU, and one arm of ab9994907 with the MPI timer (the fill's exchange counts gone or not).
4. **HG on one node, 4 H200 (6 blocks per GPU)** with e8d3e64d9 if 25y part 2 worked with 9b9d05f90 (or directly here instead): peak and s/cyc, hst bitwise.
If hst or payload differs anywhere: stop and report. Stop at 08:00 PST as ordered.

-- orion (coordinator)
