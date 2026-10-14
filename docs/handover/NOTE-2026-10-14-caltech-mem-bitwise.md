# NOTE Caltech -> raven: fin4-mem BITWISE vs FIN4 (TASK-2026-10-14-caltech-mem-bitwise) — PASS

## Binaries (plain CPU: gcc 13.2.0 + openmpi 5.0.1, Kokkos 4.6.2 OpenMP, MPI on, Release, built_in_pgens; identical cmake flags)
- MEM  fin4-mem @ 1f4e1724: athena_cpu_omp_1f4e1724  md5 **1c15546fba4e55589bbdbed2b5e987de** (build job 4490171, /resnick/groups/carnegie_poc/jingze/mem_1014)
- FIN4 fin4-1014 @ d505d21f: athena_cpu_omp_d505d21f md5 **bba193085d4b6590eda1c708c7f51dad**

## list_mem.txt (20 dirs), debug job 4490172 (3:36), `bash run_mem.sh <MEM> <FIN4> 4`, MB/EB defaults unchanged
`COUNTS:      20 BITWISE`
- restart arms: sp/rsw_T/HFfZ_1841da_1, cart/cb/FIN3_r4, sp/core_L16_u/HFfZ_1841da_128 -> restart: BITWISE (3/3)
- non-BITWISE lines: none. Full list: caltech_mem/mem_cmp.txt

## Optional list_b145.txt (145 dirs), parallel 8, debug jobs 4490189 (TIMEOUT at 97/145, cached cmp.txt reused) + 4492370 (10:13, finished)
`COUNTS:     145 BITWISE`
- non-BITWISE lines: none (the 2 expected p2 FAILs did not occur here). Full list: caltech_mem/mem_cmp_b145.txt
