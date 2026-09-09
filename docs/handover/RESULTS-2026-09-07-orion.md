# Orion results, 2026-09-07 (answers to HANDOVER-2026-09-07.md)

Everything below was done on orion at 5e9b1273 with a CPU build; orion has no GPU.
Run output is under `/orion/ptmp/jinma/Athenak/`, which viper cannot see, so the numbers
that matter are transcribed here.

## 1. The two radiation regression tests PASS on orion

`test_rad_dhj_ck_cpu` (45 s) and `test_rad_cs_raddiff_cpu` (36 s), both with a plain
`cmake -D PROBLEM=... -D CMAKE_BUILD_TYPE=Release` build under gcc 13.

Trap: orion's default `python3` is 3.6 and `tst/run_test_suite.py` dies in it before
building anything, with `Popen() got an unexpected keyword argument 'text'`. Load
`anaconda/3/2023.03` (python 3.10, pytest 7.1.2). Toolchain that works:

    module purge; module load gcc/13 cmake/3.28 anaconda/3/2023.03
    export CC=gcc CXX=g++

## 2. Determinism (handover item 4): the CPU restart path is BIT-IDENTICAL

`/orion/ptmp/jinma/Athenak/cs_determ`, MPI + OpenMP build, one node. A cubed-sphere MHD
run on the 09-07 defaults (WB polytropic + rot_potential + radiative blend + ck table,
cs 128 x 32 x 32, refit poly stretch, point-mass gravity, bbot 3 G, max_eta 1e13, STS
off) ran 800 cycles from scratch; three runs then restarted from that one restart file
and advanced 200 cycles:

| arm | layout | final rst | final bin | history |
| --- | --- | --- | --- | --- |
| B1 vs B2 | both 16 ranks x 7 threads | identical | identical | identical |
| B1 vs B3 | 16x7 against 28x4 | identical | identical | identical |

Identical means equal md5. **Even changing the rank x thread layout, and with it the
block decomposition and the reduction order, reproduces the run bit for bit.** So the
restart state is exact and the cs seam exchange, the EOS resistivity, the WB cache and
the radiative blend are all deterministic on the CPU path.

The 5e-6 divergence measured on viper is therefore GPU-specific -- atomics or reduction
order on the MI300A -- or a difference between the two binaries used there. It is NOT
the restart machinery. The remaining test is the same B1/B2 pair on one GPU binary.

Caveat: viper's production input is not in git, so this input was reconstructed from the
memory notes. It exercises every code path the handover named, but it is not byte-for-byte
viper's.

## 3. The EOS inversion table is regenerated (handover item 5)

`/orion/ptmp/jinma/Athenak/eos_table/eos_table.txt` with a `README.md` giving the exact
command. Production composition, 281 x 451 nodes. Validated against the run's own
`problem/ck_dump_file` column: **max relative error 7.6e-4 in T and 7.0e-4 in p**
(median ~7e-5) over T = 3318-6523 K and p = 8e-4 to 247 bar. The residual is the
`scripts/dhjcs.py` reader's bilinear interpolation against the run's Hermite evaluation;
regenerating four times finer does not change it.

**TRAP, and it cost most of that task: the ck column dump's `i`, `j` AND `k` are all
ghost-inclusive.** The header prints `k = 2, j = 5 ... (is = 2, ie = 65)`, which flags
the radial index but not the other two. With nghost = 2 that column is
`[0, k-2, j-2, :]` of the binary dump. Using j and k literally picks a neighbouring
column and yields a smooth, monotonic, entirely plausible profile that is wrong by 1.5 %
in p and 0.1 % in T -- it reads exactly like table interpolation error and is not one.
It survived a 4x grid refinement and two composition A/Bs before an all-columns scan
found the real one. Subtract nghost from all three.

## Still open

Handover item 3 (a GPU build carrying the NaN guard d3d74f2b) needs viper. Items 1, 2
and 5 of "Open questions" are waiting on viper's own runs after the 09-12 maintenance.
