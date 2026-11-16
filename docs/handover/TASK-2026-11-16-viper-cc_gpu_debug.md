# TASK 2026-11-16 Orion -> viper (USER): debug the sp/cs clean-solve divergence on GPU

Re NOTE-2026-11-16-viper-orion-cc-gpu-g14-curvilinear-fail (050a79319). The user asks viper to run the debugging
on the GPUs. Cartesian is exact on both GPUs; sp and cs differ from the CPU with different symptoms on MI300A and
A100, which points to uninitialised device memory, a race, or a missing host->device copy in a path only the
curvilinear cell layout uses. Orion is auditing the code in parallel and will push any fix as a new origin/clean-cc.

Code: origin/clean-cc **094488c1f** (push 11; sp/cs paths unchanged w.r.t. the bug, use it for everything below).
Decks: task_1116_cc_gpu1/decks (g3_spthin, g4_cs_thin1, g3_rhse_sp2), 1 GPU each; nlim 2 is enough unless noted.

## D1. CUDA compute-sanitizer on A100 (most informative)
Build a Release+lineinfo binary (`-D CMAKE_CXX_FLAGS="-ffp-contract=off -lineinfo"`, otherwise the A100 recipe),
then for each of the three decks with nlim 2 (and g3_rhse_sp2 with nlim 3, since its A100 failure starts at cycle 1):
```
compute-sanitizer --tool initcheck --track-unused-memory no ./athena -i <deck> time/nlim=2
compute-sanitizer --tool racecheck ./athena -i <deck> time/nlim=1
compute-sanitizer --tool memcheck ./athena -i <deck> time/nlim=1
```
Send the first ~30 distinct reports of each tool (kernel name, file:line from lineinfo).

## D2. Kokkos bounds checking
Same commit, `-D Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON` (Release otherwise), on A100 and on MI300A, the same three
decks, nlim 2. Send any "View bounds error" messages (they name the View label).

## D3. First differing cell vs CPU
For g3_spthin and g4_cs_thin1 on MI300A (where they differ already in the first solve): run nlim 1 with a bin output
every cycle (`output1/dt` tiny or `output1/dcycle=1`, variable m1, bin_double = true) on the MI300A and send the two
dump files of cycle 0 and 1; Orion compares them with the CPU and localises the first differing cells.

## D4. MI300A with XNACK off
Rerun g3_spthin and g4_cs_thin1 on MI300A with `HSA_XNACK=0` (and the HIP build unchanged). If a host dereference of
a device View exists, it now faults instead of silently reading host memory; send the error.

Reply: NOTE-...-viper-orion-cc-gpu-debug.md with the outputs (or paths) of D1-D4. D1 and D4 first if time is short.
